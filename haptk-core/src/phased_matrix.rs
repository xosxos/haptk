use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::collections::HashMap;
use std::ops::Bound;
use std::ops::RangeBounds;
use std::ops::RangeInclusive;
use std::path::Path;
use std::path::PathBuf;
use std::sync::Arc;

use ndarray::iter::AxisIter;
use ndarray::s;
use ndarray::Array2;
use ndarray::ArrayView1;
use ndarray::ArrayView2;
use ndarray::Axis;
use serde::{Deserialize, Serialize};

use crate::error;
use crate::ploidy::Ploidy;
use crate::variant::Coord;
use crate::variant::HapVariant;

#[derive(Debug, Default, Clone, PartialEq, Serialize, Deserialize)]
pub struct ReadMetadata {
    pub file_path: PathBuf,
    pub fetch_range: (u64, u64),
    pub lookups: Vec<[bool; 2]>,
    pub indexes: Vec<usize>,
    pub contig_len: Option<u64>,
    pub sharded: bool,
    pub remove_no_alt: bool,
    pub include_indels: bool,
    pub is_genome_wide: bool,
}

#[derive(Debug, Default, Clone)]
pub struct Matrix {
    pub data: Array2<u8>,
}

impl Matrix {
    fn new(data: Array2<u8>) -> Self {
        Self { data }
    }

    pub fn haplotype(&self, sample_idx: usize, range: RangeInclusive<usize>) -> Vec<u8> {
        self.data.slice(s![sample_idx, range]).to_vec()
    }

    pub fn genotype(&self, sample_idx: usize, var_idx: usize) -> u8 {
        self.data[[sample_idx, var_idx]]
    }

    pub fn genotypes(
        &self,
        variant_idx: usize,
    ) -> ndarray::ArrayBase<ndarray::ViewRepr<&u8>, ndarray::Dim<[usize; 1]>> {
        self.data.index_axis(Axis(1), variant_idx)
    }

    pub fn nrows(&self) -> usize {
        self.data.nrows()
    }

    pub fn ncols(&self) -> usize {
        self.data.ncols()
    }

    pub fn nhaplotypes(&self) -> usize {
        self.data.nrows()
    }

    pub fn nvariants(&self) -> usize {
        self.data.ncols()
    }
}

#[derive(Debug, Default, Clone)]
pub struct PhasedMatrix {
    pub start_coord: Coord,
    pub variant_idx: usize,
    pub matrix: BTreeMap<Arc<Coord>, Matrix>,
    samples: Vec<String>,
    coords: BTreeSet<Coord>,
    pub indexer: HashMap<Coord, (Arc<Coord>, usize)>,
    pub ploidy: Ploidy,
    pub metadata: ReadMetadata,
}

impl PhasedMatrix {
    pub fn new(
        variant_idx: usize,
        start_coord: Coord,
        matrix: Array2<u8>,
        samples: Vec<String>,
        coords: BTreeSet<Coord>,
        ploidy: Ploidy,
        metadata: ReadMetadata,
    ) -> Self {
        let start = Arc::new(coords.first().unwrap().clone());

        let indexer: HashMap<Coord, (Arc<Coord>, usize)> = coords
            .iter()
            .enumerate()
            .map(|(i, v)| (v.clone(), (start.clone(), i)))
            .collect();

        let mut map = BTreeMap::new();
        map.insert(start.clone(), Matrix::new(matrix));

        Self {
            variant_idx,
            start_coord,
            matrix: map,
            samples,
            coords,
            ploidy,
            indexer,
            metadata,
        }
    }

    pub fn nsamples(&self) -> usize {
        self.samples.len()
    }

    pub fn nhaplotypes(&self) -> usize {
        self.samples.len() * *self.ploidy
    }
    pub fn ncoords(&self) -> usize {
        self.coords.len()
    }

    pub fn samples(&self) -> &Vec<String> {
        &self.samples
    }

    pub fn is_genome_wide(&self) -> bool {
        self.metadata.is_genome_wide
    }

    pub fn variant_idx(&self) -> usize {
        self.variant_idx
    }

    pub fn set_variant_idx(&mut self, variant_idx: usize) {
        let new_start_coord = self.get_coord(variant_idx).clone();
        self.set_start_coord(new_start_coord);

        self.variant_idx = variant_idx;
    }

    pub fn variant_idx_pos(&self) -> u64 {
        if !self.coords.is_empty() {
            self.start_coord.pos
        } else {
            0
        }
    }

    pub fn get_pos(&self, idx: usize) -> u64 {
        self.coords.iter().nth(idx).unwrap().pos
    }

    pub fn get_coord(&self, idx: usize) -> &Coord {
        self.coords.iter().nth(idx).unwrap()
    }

    pub fn coords(&self) -> &BTreeSet<Coord> {
        &self.coords
    }

    pub fn coords_mut(&mut self) -> &mut BTreeSet<Coord> {
        &mut self.coords
    }

    pub fn set_coords(&mut self, coords: BTreeSet<Coord>) {
        self.coords = coords;
    }

    pub fn start_coord(&self) -> &Coord {
        &self.start_coord
    }

    pub fn get_coord_idx(&self, coord: &Coord) -> usize {
        self.coords.range(..=coord).count() - 1
    }

    pub fn set_start_coord(&mut self, coord: Coord) {
        self.start_coord = coord;
    }

    pub fn get_contig(&self) -> &String {
        &self.coords.first().unwrap().contig
    }
}

// General helpers
impl PhasedMatrix {
    pub fn get_nearest_coord_by_pos(&self, pos: u64) -> &Coord {
        if let Some(idx) = self.coords.iter().position(|c| c.pos >= pos) {
            let coord = self.coords.iter().nth(idx).unwrap();

            if idx == 0 || coord.pos == pos {
                return coord;
            }

            let before = self.coords.iter().nth(idx - 1).unwrap();

            match pos - before.pos >= coord.pos - pos {
                true => coord,
                false => before,
            }
        } else {
            // return last element if no element is larger
            self.coords.last().unwrap()
        }
    }

    pub fn get_sample_name(&self, index: usize) -> String {
        self.samples.get(index / *self.ploidy).unwrap().clone()
    }

    pub fn get_sample_names(&self, indexes: &[usize]) -> Vec<String> {
        indexes.iter().map(|i| self.get_sample_name(*i)).collect()
    }

    pub fn get_idxs_for_samples(
        &self,
        samples: &[String],
    ) -> std::result::Result<Vec<usize>, error::Error> {
        let idxs: Vec<_> = self
            .samples
            .iter()
            .enumerate()
            .filter(|(_, s)| samples.contains(s))
            .flat_map(|(i, _)| {
                ((i * *self.ploidy)..(i * *self.ploidy) + *self.ploidy).collect::<Vec<usize>>()
            })
            .collect();

        if idxs.is_empty() {
            return Err(error::Error::SamplesNotFound);
        }

        Ok(idxs)
    }

    // Inclusive ranges not supported so remember to add + 1 to stop_idx
    pub fn find_haplotype_for_sample<R: RangeBounds<Coord> + std::fmt::Debug>(
        &self,
        range: R,
        sample: usize,
    ) -> Vec<HapVariant> {
        self.coords()
            .range(range)
            .map(|coord| {
                let gt = self.matrix_point_coord(sample, coord);

                HapVariant {
                    contig: coord.contig.to_string(),
                    pos: coord.pos,
                    alt: coord.alt.clone(),
                    reference: coord.reference.clone(),
                    gt,
                }
            })
            .collect()
    }

    // Inclusive ranges not supported so remember to add + 1 to stop_idx
    pub fn find_u8_haplotype_for_sample<R: RangeBounds<Coord>>(
        &self,
        col_range: R,
        sample: usize,
    ) -> Vec<u8> {
        let (col_first_idx, col_last_idx) = match (col_range.start_bound(), col_range.end_bound()) {
            (Bound::Included(first), Bound::Included(last)) => (first, last),
            _ => unreachable!(
                "Range for find_u8_haplotype_for_sample should always be inclusive x..=y"
            ),
        };
        let (matrix_key_first, index_first) = &self.indexer[col_first_idx];
        let (matrix_key_last, index_last) = &self.indexer[col_last_idx];

        if matrix_key_first == matrix_key_last {
            let matrix = self.matrix.get(matrix_key_first).unwrap();
            return matrix.haplotype(sample, *index_first..=*index_last);
        }

        self.matrix
            .range(matrix_key_first.clone()..matrix_key_last.clone())
            .flat_map(|(key, next_matrix)| {
                let index_stop = match key == matrix_key_last {
                    true => *index_last,
                    false => next_matrix.ncols() - 1,
                };

                let index_start = match key == matrix_key_first {
                    true => *index_first,
                    false => 0,
                };

                next_matrix.haplotype(sample, index_start..=index_stop)
            })
            .collect()
    }
}

// Inner matrix methods
impl PhasedMatrix {
    pub fn matrix_nrows(&self) -> usize {
        self.matrix.values().map(|m| m.nrows()).sum()
    }

    pub fn matrix_ncols(&self) -> usize {
        self.matrix.values().map(|m| m.ncols()).sum()
    }

    pub fn matrix_point_coord(&self, x: usize, coord: &Coord) -> u8 {
        match self.indexer.get(coord) {
            Some((matrix_key, index)) => {
                let matrix = self.matrix.get(matrix_key).unwrap();
                matrix.data[[x, *index]]
            }
            None => {
                let nearest_coord = self.get_nearest_coord_by_pos(coord.pos);
                tracing::warn!(
                    "Could not find coord {}, using {} instead",
                    coord,
                    nearest_coord
                );
                let (matrix_key, index) = self.indexer.get(nearest_coord).unwrap();
                let matrix = self.matrix.get(matrix_key).unwrap();
                matrix.data[[x, *index]]
            }
        }
    }

    pub fn matrix_column(&self, coord: &'_ Coord) -> ArrayView1<'_, u8> {
        let (matrix, index) = self
            .indexer
            .get(coord)
            .unwrap_or_else(|| panic!("Could not find coord {coord} from the indexer"));

        let matrix = self.matrix.get(matrix).unwrap();

        matrix.genotypes(*index)
    }

    pub fn insert_matrix(&mut self, start_coord: Coord, matrix: Array2<u8>) {
        self.matrix
            .insert(Arc::new(start_coord), Matrix::new(matrix));
    }

    // Sort inside select_rows to avoid bugs down the line
    pub fn select_rows(&mut self, mut to_keep: Vec<usize>) {
        to_keep.sort();
        self.matrix_select(0, &to_keep);

        self.samples = to_keep
            .iter()
            .map(|index| self.get_sample_name(*index))
            .collect();
    }

    pub fn matrix_select(&mut self, axis: usize, to_keep: &[usize]) {
        if axis == 0 {
            for (_, matrix) in self.matrix.iter_mut() {
                *matrix = Matrix::new(matrix.data.select(Axis(0), to_keep))
            }
        } else {
            unreachable!("matrix select on columns is not yet implemented")
        }
    }

    // NOTE: Legacy compare-to-haplotype methods, will be deprecated when index based access to matrices is completely removed
    pub fn matrix_axis_iter(&self, axis: usize) -> AxisIter<'_, u8, ndarray::Dim<[usize; 1]>> {
        let (_, matrix) = self.matrix.iter().nth(0).unwrap();
        matrix.data.axis_iter(Axis(axis))
    }

    pub fn set_matrix(&mut self, start_coord: Coord, matrix: Array2<u8>) {
        let mut map = BTreeMap::new();

        map.insert(Arc::new(start_coord), Matrix::new(matrix));
        self.matrix = map;
    }

    pub fn slice_cols<R: RangeBounds<usize>>(&'_ self, col_range: R) -> ArrayView2<'_, u8> {
        let (col_first_idx, col_last_idx) = match (col_range.start_bound(), col_range.end_bound()) {
            (Bound::Included(first), Bound::Excluded(last)) => (*first, last.saturating_sub(1)),
            (Bound::Included(first), Bound::Included(last)) => (*first, *last),
            _ => unreachable!(
                "Range for find_u8_haplotype_for_sample should always be x..=y or x..y"
            ),
        };

        let (_, matrix) = self.matrix.iter().nth(0).unwrap();

        matrix.data.slice(s![.., col_first_idx..=col_last_idx])
    }

    pub fn write_npy(&self, path: &Path) -> std::result::Result<(), error::Error> {
        let (_, matrix) = self.matrix.iter().nth(0).unwrap();

        let mut writer = std::io::BufWriter::new(std::fs::File::create(path)?);
        npy_v1::write(&mut writer, &matrix.data)?;
        Ok(())
    }
}

mod npy_v1 {
    use byteorder::{BigEndian, ByteOrder, LittleEndian, NativeEndian, WriteBytesExt};
    use ndarray::prelude::*;
    use ndarray::Data;
    use std::io;

    static MAGIC_VALUE: &[u8] = b"\x93NUMPY";
    /// npy Version 1.0
    static NPY_VERSION: &[u8] = b"\x01\x00";

    /// Types that can be serialized using this crate.
    pub trait DType<B> {
        fn dtype() -> &'static str;
        fn write_bytes<W: std::io::Write>(self, w: &mut W) -> io::Result<()>;
    }

    impl<B: ByteOrder> DType<B> for u8 {
        fn dtype() -> &'static str {
            "u1"
        }

        fn write_bytes<W: std::io::Write>(self, w: &mut W) -> io::Result<()> {
            w.write_u8(self)
        }
    }

    trait NumpyEndian {
        fn endian_symbol() -> &'static str;
    }

    impl NumpyEndian for LittleEndian {
        fn endian_symbol() -> &'static str {
            "<"
        }
    }

    impl NumpyEndian for BigEndian {
        fn endian_symbol() -> &'static str {
            ">"
        }
    }

    fn get_header<A, B>(shape: &[usize]) -> String
    where
        A: DType<B>,
        B: NumpyEndian,
    {
        use std::fmt::Write;
        let mut shape_str = String::new();

        for (i, s) in shape.iter().enumerate() {
            if i > 0 {
                shape_str.push(',');
            }
            write!(&mut shape_str, "{}", s).unwrap();
        }

        format!(
            "{{'descr': '{endian}{dtype}','fortran_order': False,'shape': ({shape})}}\n",
            endian = B::endian_symbol(),
            dtype = A::dtype(),
            shape = shape_str
        )
    }

    /// Write an ndarray to a writer in the numpy format.
    ///
    /// Can be saved with file extension `npy` and loaded using `numpy.load`.
    pub fn write<A, S, D, T>(w: &mut T, array: &ArrayBase<S, D>) -> io::Result<()>
    where
        S: Data<Elem = A>,
        D: Dimension,
        A: DType<NativeEndian> + Copy,
        T: std::io::Write,
    {
        let header = get_header::<A, NativeEndian>(array.shape());
        let padding = (16 - (MAGIC_VALUE.len() + 4 + header.len()) % 16) % 16;
        let header_len = header.len() + padding;
        // the following value must be divisible by 16
        assert_eq!(
            (MAGIC_VALUE.len() + 4 + header_len) % 16,
            0,
            "Invalid alignment of the npy header"
        );
        assert!(
            header_len <= u16::MAX as usize,
            "Length of the npy header overflowed."
        );

        w.write_all(MAGIC_VALUE)?;
        w.write_all(NPY_VERSION)?;
        w.write_u16::<LittleEndian>(header_len as u16)?;
        w.write_all(header.as_bytes())?;

        // padding
        for _ in 0..padding {
            w.write_u8(b' ')?;
        }

        // actual data
        for x in array.iter() {
            x.write_bytes(w)?;
        }

        Ok(())
    }
}
