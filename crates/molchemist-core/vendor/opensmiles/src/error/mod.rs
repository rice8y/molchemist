//! Error types for the SMILES parser.
//!
//! This module contains all error types used by the SMILES parser,
//! organized by component:
//!
//! - [`AtomError`] - Atom-related errors (charge, isotope, element)
//! - [`NodeError`] - Node-related errors (hydrogens, class, aromaticity)
//! - [`MoleculeError`] - Molecule construction errors
//! - [`ParserError`] - SMILES string parsing errors
//!
//! # Error hierarchy
//!
//! ```text
//! ParserError
//! ├── MoleculeError
//! │   ├── NodeError
//! │   │   └── AtomError
//! │   └── AtomError
//! └── NodeError
//!     └── AtomError
//! ```
//!
//! `From` conversions are implemented to allow using the `?` operator
//! throughout the hierarchy.

mod atom;
mod bond;
mod molecule;
mod node;
mod parser;

pub use atom::AtomError;
pub use bond::BondError;
pub use molecule::MoleculeError;
pub use node::NodeError;
pub use parser::ParserError;
