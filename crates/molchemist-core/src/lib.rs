//! Shared conversion engine used by molchemist's Typst plugins and CLI.

pub type NodeIndex = u32;

#[path = "../vendor/opensmiles/src/ast/mod.rs"]
pub mod ast;
#[path = "../vendor/opensmiles/src/error/mod.rs"]
mod error;
#[path = "../vendor/opensmiles/src/parser.rs"]
pub mod parser;

pub use ast::*;
pub use error::*;
pub use parser::*;

mod engine;
mod formatter;

#[cfg(feature = "native-layout")]
mod native_layout;

pub use engine::{
    sdf_to_ast, sdf_to_commands, smiles_layout_input, smiles_to_ast,
    smiles_to_commands_with_coords, smiles_to_full_layout_input, smiles_to_layout_input, Command,
    LinkData, RenderMode,
};
pub use formatter::{
    format_alchemist, format_standalone, format_standalone_code, StandaloneOptions,
    DEFAULT_ALCHEMIST_IMPORT,
};

#[cfg(feature = "native-layout")]
pub use native_layout::{layout_payload, smiles_to_commands};
