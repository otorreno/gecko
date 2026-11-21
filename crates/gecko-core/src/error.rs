use std::path::PathBuf;
use thiserror::Error;

/// Errors returned by `gecko-core` operations.
#[derive(Debug, Error)]
pub enum GeckoError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Invalid parameter: {0}")]
    InvalidParameter(String),

    #[error("Invalid FASTA format in {path}: {message}")]
    InvalidFasta { path: PathBuf, message: String },

    #[error("Invalid binary format: {0}")]
    InvalidFormat(String),

    #[error("JSON serialization error: {0}")]
    Json(#[from] serde_json::Error),
}

pub type Result<T> = std::result::Result<T, GeckoError>;
