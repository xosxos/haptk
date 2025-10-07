use serde::{Deserialize, Serialize};

#[derive(Debug, Default, Clone, PartialEq, Serialize, Deserialize)]
pub enum Ploidy {
    Mixed,
    Haploid,
    #[default]
    Diploid,
}

impl From<Ploidy> for usize {
    fn from(ploidy: Ploidy) -> Self {
        match ploidy {
            Ploidy::Haploid => 1,
            Ploidy::Mixed => 1,
            Ploidy::Diploid => 2,
        }
    }
}

impl From<&Ploidy> for usize {
    fn from(ploidy: &Ploidy) -> Self {
        match ploidy {
            Ploidy::Haploid => 1,
            Ploidy::Mixed => 1,
            Ploidy::Diploid => 2,
        }
    }
}

impl std::ops::Deref for Ploidy {
    type Target = usize;
    fn deref(&self) -> &Self::Target {
        match self {
            Ploidy::Haploid => &1,
            Ploidy::Mixed => &1,
            Ploidy::Diploid => &2,
        }
    }
}
