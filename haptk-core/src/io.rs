use std::ffi::OsStr;
use std::fs;
use std::fs::File;
use std::io;
use std::io::BufRead;
use std::path::Path;
use std::path::PathBuf;

use csv::QuoteStyle;
use csv::Reader;
use csv::ReaderBuilder;
use csv::Writer;
use csv::WriterBuilder;

use crate::error::Error;

#[allow(clippy::upper_case_acronyms)]
pub enum FileType {
    BED,
    VCF,
    CSV,
    TSV,
    HST,
    BAM,
    CRAM,
}

impl FileType {
    pub fn from_path(path: &Path) -> Result<Self, Error> {
        let ext = get_extension(path)?;

        Ok(match ext.as_str() {
            "vcf.gz" | "vcf" | "bcf" => Self::VCF,
            "cram" => Self::CRAM,
            "bam" => Self::BAM,
            "hst.gz" | "hst" => Self::HST,
            "bed.gz" | "bed" => Self::BED,
            "csv" | "ids" => Self::CSV,
            "tsv" => Self::TSV,
            _ => return Err(Error::FileNotSupported { ext }),
        })
    }
}

pub fn get_extension(path: &Path) -> Result<String, Error> {
    fn double_extension(path: &Path, e1: &str) -> Result<String, Error> {
        let stem = path
            .file_stem()
            .and_then(OsStr::to_str)
            .ok_or_else(|| Error::NoFileType {
                path: path.to_path_buf(),
            })?;

        let e2 = Path::new(&stem)
            .extension()
            .and_then(OsStr::to_str)
            .ok_or_else(|| Error::FileNotSupported {
                ext: format!("{path:?}"),
            })?;

        Ok(format!("{e2}.{e1}"))
    }

    let extension: &str = Path::new(&path)
        .extension()
        .and_then(OsStr::to_str)
        .ok_or_else(|| Error::NoFileType {
            path: path.to_path_buf(),
        })?;

    let ext = match extension {
        "gz" | "bgz" => double_extension(path, extension)?,
        _ => extension.to_string(),
    };

    Ok(ext)
}

pub fn get_csv_reader<R: io::Read>(input: R, has_headers: bool) -> Reader<R> {
    ReaderBuilder::new()
        .quote(b'"')
        .delimiter(b',')
        .has_headers(has_headers)
        .flexible(false)
        .from_reader(input)
}

pub fn get_csv_writer<W: io::Write>(output: W) -> Writer<W> {
    WriterBuilder::new()
        .quote(b'"')
        .delimiter(b',')
        .has_headers(false)
        .flexible(true)
        .from_writer(output)
}

pub fn get_vcf_writer<W: io::Write>(output: W) -> Writer<W> {
    WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .flexible(true)
        .double_quote(false)
        .quote_style(QuoteStyle::Never)
        .from_writer(output)
}

pub fn read_multiple_sample_ids(path: &Option<Vec<PathBuf>>) -> Result<Option<Vec<String>>, Error> {
    match path {
        Some(paths) => {
            let mut samples = vec![];

            for path in paths {
                for line in read_lines(path)?.map_while(Result::ok) {
                    let line = line.trim();
                    samples.push(line.to_string());
                }
            }
            Ok(Some(samples))
        }
        None => Ok(None),
    }
}

pub fn read_lines<P>(filename: P) -> Result<io::Lines<io::BufReader<File>>, Error>
where
    P: AsRef<Path> + Into<PathBuf>,
{
    let file = File::open(&filename).map_err(|e| Error::Io(filename.into(), e))?;

    Ok(io::BufReader::new(file).lines())
}

fn not_found(path: &Path, err: &str) -> Error {
    Error::Path(path.to_path_buf(), err.to_string())
}

pub fn get_input(filename: Option<PathBuf>) -> Result<Box<dyn io::Read>, Error> {
    let input: Box<dyn io::Read> = match filename {
        Some(name) => match name.to_str() {
            Some("-") => Box::new(std::io::stdin()),
            Some(path) => {
                let file_handle = compression::get_reader(path)?;

                Box::new(file_handle)
            }
            None => return Err(not_found(&name, "unreachable")),
        },
        None => Box::new(io::stdin()),
    };
    Ok(input)
}

pub fn get_output(filename: Option<PathBuf>) -> Result<Box<dyn io::Write>, Error> {
    let output: Box<dyn io::Write> = match filename {
        Some(name) => match name.to_str() {
            Some("-") => Box::new(io::stdout()),
            Some(path) => {
                let file_handle = fs::File::options()
                    .create(true)
                    .write(true)
                    .truncate(true)
                    .open(path)
                    .map_err(|e| not_found(&name, &e.to_string()))?;

                Box::new(file_handle)
            }
            None => unreachable!(),
        },
        None => Box::new(io::stdout()),
    };
    Ok(output)
}

mod compression {
    use crate::error;
    use std::io;
    use std::io::Read;

    #[derive(Debug, Eq, PartialEq, Clone, Copy)]
    pub enum Format {
        Gzip,
        Bzip,
        Lzma,
        Zstd,
        No,
    }

    pub fn get_reader<'a>(path: &str) -> Result<Box<dyn io::Read + 'a>, error::Error> {
        let mut in_stream = Box::new(io::BufReader::new(std::fs::File::open(path)?));

        // Check compression
        let mut first_bytes = [0u8; 5];
        in_stream
            .read_exact(&mut first_bytes)
            .map_err(|_| error::Error::FileTooShort)?;

        let compression = match first_bytes {
            [0x1f, 0x8b, ..] => Format::Gzip,
            [0x42, 0x5a, ..] => Format::Bzip,
            [0x28, 0xb5, 0x2f, 0xfd, ..] => Format::Zstd,
            [0xfd, 0x37, 0x7a, 0x58, 0x5a] => Format::Lzma,
            _ => Format::No,
        };

        let in_stream = Box::new(io::Cursor::new(first_bytes).chain(in_stream));

        // Return decoding stream
        match compression {
            Format::Gzip => Ok(Box::new(flate2::read::MultiGzDecoder::new(in_stream))),
            Format::Bzip => Ok(Box::new(bzip2::read::MultiBzDecoder::new(in_stream))),
            Format::Lzma => Ok(Box::new(liblzma::read::XzDecoder::new(in_stream))),
            Format::Zstd => Ok(Box::new(zstd::stream::read::Decoder::new(in_stream)?)),
            Format::No => Ok(in_stream),
        }
    }
}
