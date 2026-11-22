pub mod error;
pub mod extend;
pub mod formats;
pub mod hits;
pub mod kmer;
pub mod pipeline;
pub mod radix;
pub mod ref_index;
pub mod seq;
pub mod stream;

pub use error::{GeckoError, Result};
pub use extend::{CoverageMatrix, ExtendOptions, Scoring};
pub use formats::Fragment;
pub use hits::{Hit, HitOptions, HitStats};
pub use kmer::{DictOptions, DictStats, DiskKmerIndex, INDEX_WORD_SIZE, KmerDict, KmerIndex};
pub use pipeline::{
    CompareConfig, CompareResult, CompareStats, compare_genomes, compare_with_dicts,
    compare_with_ref_index,
};
pub use ref_index::{PosStrandVec, RefIndex};
pub use seq::{Genome, PackedSeq, SeqRecord};
