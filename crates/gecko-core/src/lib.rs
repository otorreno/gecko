pub mod error;
pub mod extend;
pub mod formats;
pub mod hits;
pub mod kmer;
pub mod pipeline;
pub mod radix;
pub mod seq;

pub use error::{GeckoError, Result};
pub use extend::{CoverageMatrix, ExtendOptions, Scoring};
pub use formats::Fragment;
pub use hits::{Hit, HitOptions, HitStats};
pub use kmer::{DictOptions, DictStats, INDEX_WORD_SIZE, KmerDict};
pub use pipeline::{
    CompareConfig, CompareResult, CompareStats, compare_genomes, compare_with_dicts,
};
pub use seq::{Genome, SeqRecord};
