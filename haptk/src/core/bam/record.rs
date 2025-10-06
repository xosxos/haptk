use rust_htslib::bam::record::Cigar;

pub struct Record(rust_htslib::bam::Record);

impl Record {
    pub fn new(record: rust_htslib::bam::Record) -> Self {
        Self(record)
    }
    pub fn seq(&self) -> Vec<u8> {
        self.0.seq().as_bytes()
    }

    pub fn pos(&self) -> u64 {
        self.0.pos() as u64
    }

    pub fn aux(&self) -> Vec<String> {
        fn fold_func<S: ToString>(acc: String, cur: S) -> String {
            match acc.is_empty() {
                true => cur.to_string(),
                false => format!("{acc},{}", cur.to_string()),
            }
        }

        self.0
            .aux_iter()
            .map(|v| v.unwrap())
            .map(|(tag, v)| {
                let tag = String::from_utf8(tag.to_vec()).unwrap();
                match v {
                    rust_htslib::bam::record::Aux::Char(i) => i.to_string(),
                    rust_htslib::bam::record::Aux::I8(i) => format!("{tag}:i:{i}"),
                    rust_htslib::bam::record::Aux::U8(i) => format!("{tag}:i:{i}"),
                    rust_htslib::bam::record::Aux::I16(i) => format!("{tag}:i:{i}"),
                    rust_htslib::bam::record::Aux::U16(i) => format!("{tag}:i:{i}"),
                    rust_htslib::bam::record::Aux::I32(i) => format!("{tag}:i:{i}"),
                    rust_htslib::bam::record::Aux::U32(i) => format!("{tag}:i:{i}"),
                    rust_htslib::bam::record::Aux::Float(i) => format!("{tag}:f:{i}"),
                    rust_htslib::bam::record::Aux::Double(i) => format!("{tag}:f:{i}"),
                    rust_htslib::bam::record::Aux::String(i) => format!("{tag}:Z:{i}"),
                    rust_htslib::bam::record::Aux::HexByteArray(i) => format!("{tag}:H:{i}"),
                    rust_htslib::bam::record::Aux::ArrayI8(i) => {
                        format!("{tag}:B:i,{}", i.iter().fold(String::new(), fold_func))
                    }
                    rust_htslib::bam::record::Aux::ArrayU8(i) => {
                        format!("{tag}:B:i,{}", i.iter().fold(String::new(), fold_func))
                    }
                    rust_htslib::bam::record::Aux::ArrayI16(i) => {
                        format!("{tag}:B:i,{}", i.iter().fold(String::new(), fold_func))
                    }
                    rust_htslib::bam::record::Aux::ArrayU16(i) => {
                        format!("{tag}:B:i,{}", i.iter().fold(String::new(), fold_func))
                    }
                    rust_htslib::bam::record::Aux::ArrayI32(i) => {
                        format!("{tag}:B:i,{}", i.iter().fold(String::new(), fold_func))
                    }
                    rust_htslib::bam::record::Aux::ArrayU32(i) => {
                        format!("{tag}:B:i,{}", i.iter().fold(String::new(), fold_func))
                    }
                    rust_htslib::bam::record::Aux::ArrayFloat(i) => {
                        format!("{tag}:B:f,{}", i.iter().fold(String::new(), fold_func))
                    }
                }
            })
            .collect()
    }

    pub fn cigar(&self) -> Vec<Cigar> {
        self.0.cigar().take().to_vec()
    }

    pub fn cigar_string(&self) -> String {
        self.0.cigar().to_string()
    }

    pub fn _mtid(&self) -> i32 {
        self.0.mtid()
    }

    pub fn _mpos(&self) -> i64 {
        self.0.mpos()
    }

    pub fn qname(&self) -> String {
        unsafe { String::from_utf8_unchecked(self.0.qname().to_vec()) }
    }

    pub fn flags(&self) -> u16 {
        self.0.flags()
    }

    pub fn mapq(&self) -> u8 {
        self.0.mapq()
    }
}
