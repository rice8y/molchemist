//! # opensmiles
//!
//! A SMILES parser following the [OpenSMILES specification](http://opensmiles.org/opensmiles.html).
//!
//! ## Quick Start
//!
//! ```rust
//! use opensmiles::parse;
//!
//! let molecule = parse("CCO").unwrap(); // ethanol
//! println!("{}", molecule); // OCC
//! ```
//!
//! ## Parallel Parsing
//!
//! Enable the `parallel` feature for multi-threaded batch parsing:
//!
//! ```rust
//! # #[cfg(feature = "parallel")]
//! # {
//! use opensmiles::parse_batch;
//!
//! let smiles = vec!["CCO", "c1ccccc1", "CC(=O)O"];
//! let results = parse_batch(&smiles);
//! assert_eq!(results.len(), 3);
//! # }
//! ```
//!
//! ## Grammar
//!
//! This parser implements an LL(1) grammar based on the formal grammar
//! described at <https://depth-first.com/articles/2020/12/21/smiles-formal-grammar-revisited/>

/// Index type for atoms/nodes in a molecule graph.
///
/// Using `u32` allows molecules with up to ~4 billion atoms,
/// sufficient for any polymer simulation.
pub type NodeIndex = u32;

pub mod ast;
mod error;
pub mod parser;
#[cfg(feature = "parallel")]
pub mod parser_parallel;

// Re-export public API
pub use ast::*;
pub use error::*;
pub use parser::*;
#[cfg(feature = "parallel")]
pub use parser_parallel::*;
