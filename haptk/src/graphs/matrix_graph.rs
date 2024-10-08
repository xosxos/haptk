// use color_eyre::Result;
use svg::node::element::path::Data;
use svg::node::element::Element;
use svg::node::element::Path;
use svg::node::Text;
use svg::Document;
use svg::Node;

use crate::args::GraphArgs;
use crate::structs::MatrixSlice;
use crate::structs::PhasedMatrix;

pub fn determine_line_color(
    vcf: &PhasedMatrix,
    decoy_samples: &Option<Vec<String>>,
    longest_alleles: &Option<Vec<usize>>,
    row_idx: usize,
    variant_idx: usize,
) -> image::Rgba<u8> {
    // Color based on if only longest or decoy samples coloring is wanted
    let sample = vcf.get_sample_name(row_idx);
    if variant_idx == vcf.variant_idx() {
        return image::Rgba([0, 0, 0, 255]);
    }
    match (decoy_samples, longest_alleles) {
        (Some(vec), None) => match vec.contains(&sample) {
            true => image::Rgba([0, 255, 0, 255]),
            false => image::Rgba([255, 255, 255, 255]),
        },
        (None, Some(vec)) => match vec.contains(&row_idx) {
            true => image::Rgba([255, 255, 255, 255]),
            false => image::Rgba([0, 255, 255, 255]),
        },
        (Some(vec), Some(_)) => match vec.contains(&sample) {
            true => image::Rgba([0, 255, 0, 255]),
            false => image::Rgba([255, 255, 255, 255]),
        },
        (None, None) => image::Rgba([255, 255, 255, 255]),
    }
}

#[derive(Clone, Debug)]
pub struct Point {
    gt: u8,
    y: f32,
    x: f32,
}

impl std::fmt::Display for Point {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let line = format!("Point {{ x: {}, y: {} }}", self.x, self.y);
        write!(f, "{line}")
    }
}

#[derive(Debug)]
pub struct MatrixGraph<'a> {
    pub document: Document,
    pub s: GraphArgs,
    pub vcf: &'a PhasedMatrix,
    marker_width: f32,
    row_height: f32,
    btm_padding: f32,
    decoy_samples: Option<Vec<String>>,
    longest_alleles: Option<&'a Vec<usize>>,
}

impl<'a> MatrixGraph<'a> {
    pub fn new(
        vcf: &'a PhasedMatrix,
        s: GraphArgs,
        decoy_samples: Option<Vec<String>>,
        mark_shorter_alleles: bool,
        only_longest: Option<&'a Vec<usize>>,
    ) -> Self {
        let document = Document::new()
            .set("viewBox", (0, 0, s.width, s.height))
            .set("style", format!("background-color:{}", s.background_color));

        let longest_alleles = match mark_shorter_alleles {
            true => only_longest,
            false => None,
        };

        Self {
            document,
            vcf,
            marker_width: 0.0,
            row_height: 0.0,
            btm_padding: s.width * 0.015,
            decoy_samples,
            longest_alleles,
            s,
        }
    }

    pub fn draw_graph(&mut self, order: &[usize]) {
        let nrows = self.vcf.matrix_nrows();
        let ncols = self.vcf.matrix_ncols();
        self.marker_width = self.s.width / ncols as f32;
        self.row_height = (self.s.height - (self.s.height) * 0.015) / nrows as f32;

        for (y, row_idx) in order.iter().enumerate() {
            let row = self
                .vcf
                .matrix_slice(MatrixSlice::Point(*row_idx), MatrixSlice::All);
            for (x, gt) in row.iter().enumerate() {
                self.create_box(
                    Point {
                        gt: *gt,
                        x: x as f32,
                        y: y as f32,
                    },
                    *row_idx,
                    x,
                );
            }
        }

        self.draw_bottom_margin();
    }

    pub fn draw_bottom_margin(&mut self) {
        let y = self.s.height - self.btm_padding / 3.;

        let start = self.vcf.coords().first().unwrap();
        let stop = self.vcf.coords().last().unwrap();

        let mut element = Element::new("text");
        element.assign("x", self.s.width * 0.02);
        element.assign("y", y);
        element.assign("fill", "black");
        element.assign("font-size", format!("{}px", self.s.font_size));
        element.append(Text::new(format!("start: {}", start.pos)));
        self.document.append(element);

        let mut element = Element::new("text");
        element.assign("x", self.s.width * 0.25);
        element.assign("y", y);
        element.assign("fill", "black");
        element.assign("font-size", format!("{}px", self.s.font_size));
        element.append(Text::new(start.contig.clone()));
        self.document.append(element);

        let mut element = Element::new("text");
        element.assign("x", self.s.width * 0.45);
        element.assign("y", y);
        element.assign("fill", "black");
        element.assign("font-size", format!("{}px", self.s.font_size));
        element.append(Text::new(format!("markers: {}", self.vcf.ncoords())));
        self.document.append(element);

        let mut element = Element::new("text");
        element.assign("x", self.s.width * 0.65);
        element.assign("y", y);
        element.assign("fill", "black");
        element.assign("font-size", format!("{}px", self.s.font_size));
        element.append(Text::new(format!("haplotypes: {}", self.vcf.nhaplotypes())));
        self.document.append(element);

        let second_place = self.s.width * 0.86;
        let mut element = Element::new("text");
        element.assign("x", second_place);
        element.assign("y", y);
        element.assign("fill", "black");
        element.assign("font-size", format!("{}px", self.s.font_size));
        element.append(Text::new(format!("stop: {}", stop.pos)));
        self.document.append(element);
    }

    pub fn create_box(&mut self, p: Point, row_idx: usize, var_idx: usize) {
        let color = self.determine_line_color(&p, row_idx, var_idx);

        let x = (p.x + 1.0) * self.marker_width;
        let y = (p.y + 1.0) * self.row_height;

        let vertical_line = Data::new()
            .move_to((x, y))
            .line_to((x, y + self.row_height / 1.5));

        let vertical_line = Path::new()
            .set("fill", "none")
            .set("stroke", color)
            .set("opacity", 1)
            .set("stroke-width", self.s.stroke_width)
            .set("d", vertical_line);

        self.document.append(vertical_line);
    }

    pub fn determine_line_color(&self, p: &Point, row_idx: usize, var_idx: usize) -> &str {
        // Color based on if only longest or decoy samples coloring is wanted
        let sample = self.vcf.get_sample_name(row_idx);

        if var_idx == self.vcf.variant_idx() && p.gt == 1 {
            return "#000";
        }

        if row_idx == self.vcf.nhaplotypes() / 2 {
            return "#000";
        }

        match p.gt {
            0 => "#ff0087",
            1 => match (&self.decoy_samples, &self.longest_alleles) {
                (Some(vec), None) => match vec.contains(&sample) {
                    true => "#00ff00",
                    false => "#fff",
                },
                (None, Some(vec)) => match vec.contains(&row_idx) {
                    true => "#fff",
                    false => "#00ffff",
                },
                (Some(vec), Some(_)) => match vec.contains(&sample) {
                    true => "#00ff00",
                    false => "#fff",
                },
                (None, None) => "#fff",
            },
            _ => panic!(),
        }
    }
}

// pub fn matrix_graph_png(
//     vcf: &mut PhasedMatrix,
//     s: GraphArgs,
//     mark_shorter_alleles: bool,
//     decoy_samples: Option<Vec<String>>,
//     order: &[usize],
// ) -> Result<image::RgbaImage> {
//     let longest_alleles: Option<Vec<usize>> = match mark_shorter_alleles {
//         true => Some(vcf.only_longest_indexes()),
//         false => None,
//     };

//     let width = vcf.matrix_ncols() * 5;
//     let height = vcf.matrix_nrows() * 7;

//     let mut imgbuf: image::RgbaImage = image::ImageBuffer::new(width as u32, height as u32);

//     // for y in 0..vcf.matrix.nrows() {
//     for (y, row_idx) in order.iter().enumerate() {
//         let row = vcf.matrix_slice(MatrixSlice::Point(*row_idx), MatrixSlice::All);

//         for (x, gt) in row.iter().enumerate() {
//             let ny = y * 7;
//             let nx = x * 5;
//             for yy in ny..ny + 7 {
//                 for xx in nx..nx + 5 {
//                     let pixel = imgbuf.get_pixel_mut(xx as u32, yy as u32);
//                     if (yy == ny || yy == ny + 7) || (xx == nx || xx == nx + 5) {
//                         *pixel = image::Rgba([255, 255, 255, 255]);
//                     } else {
//                         match gt {
//                             0 => *pixel = image::Rgba([255, 0, 135, 255]),
//                             1 => {
//                                 *pixel = determine_line_color(
//                                     vcf,
//                                     &decoy_samples,
//                                     &longest_alleles,
//                                     y,
//                                     x,
//                                 )
//                             }
//                             _ => panic!(),
//                         }
//                     }
//                 }
//             }
//         }
//     }
//     Ok(image::imageops::resize(
//         &imgbuf,
//         s.width as u32,
//         s.height as u32,
//         image::imageops::FilterType::CatmullRom,
//     ))
// }

#[cfg(test)]
#[rustfmt::skip]
mod tests {
    use super::*;

    #[test]
    fn point_display() {
        let point = Point { x: 5.0, y: 6.0, gt: 1 };
        assert_eq!("Point { x: 5, y: 6 }".to_string(), format!("{point}"));
    }
}
