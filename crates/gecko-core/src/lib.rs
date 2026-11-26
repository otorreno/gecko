pub mod error;
pub mod extend;
pub mod formats;
pub mod hits;
pub mod kmer;
pub mod minimizer;
pub mod pipeline;
pub mod radix;
pub mod ref_index;
pub mod seq;
pub mod spaced;
pub mod stream;

pub use error::{GeckoError, Result};
pub use extend::{CoverageMatrix, ExtendOptions, Scoring};
pub use formats::Fragment;
pub use hits::{Hit, HitOptions, HitStats};
pub use kmer::{DictOptions, DictStats, DiskKmerIndex, INDEX_WORD_SIZE, KmerDict, KmerIndex};
pub use minimizer::{MINIMIZER_MIN_K, MinimizerIndex, compare_minimizer_streaming_with_coverage};
pub use pipeline::{
    CompareConfig, CompareResult, CompareStats, SeedingMode, compare_genomes,
    compare_genomes_with_coverage, compare_with_dicts, compare_with_ref_index,
    compare_with_ref_index_with_coverage, profile_enabled,
};
pub use ref_index::{PosStrandVec, RefIndex};
pub use seq::{Genome, PackedSeq, SeqRecord};
pub use spaced::SpacedSeed;
pub use stream::compare_streaming_with_coverage;
