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
mod formats;
mod formatter;
mod semantic;
pub use formats::{mol2_to_sdf, parse_reaction, MoleculeInput, ReactionInput};
mod reaction;
mod rgroups;
mod stereo3d;
mod valence;
pub use reaction::{analyze_reaction, ReactionAnalysis};
pub use rgroups::{expand_superatoms, parse_rgroups};

#[cfg(feature = "native-layout")]
mod native_layout;

pub use engine::{
    inspect_sdf_record, sdf_depiction_positions, sdf_force_layout_input,
    sdf_projection_layout_input, sdf_record_to_ast, sdf_record_to_ast_with_coords,
    sdf_record_to_commands, sdf_record_to_commands_with_coords, sdf_record_to_inspection_cbor,
    sdf_record_to_layout_input, sdf_reoriented_ast, sdf_stereo3d_ast, sdf_to_ast, sdf_to_commands,
    smiles_layout_input, smiles_to_ast, smiles_to_commands_with_coords,
    smiles_to_full_layout_input, smiles_to_layout_input, AtomLabel, Command,
    CtfileAtomAnnotationDepiction, CtfileAtomQueryDepiction, CtfileBondQueryDepiction,
    CtfileHighlightBond, CtfileHighlightDepiction, CtfilePoint, CtfileSGroupDepiction,
    CtfileVariableAttachmentDepiction, LinkData, RenderMode,
};
pub use formatter::{composition_code, format_coordinate_code};
pub use semantic::{
    ChemicalAtom, ChemicalAtomQuery, ChemicalBond, ChemicalBondOrder, ChemicalBondStereo,
    ChemicalCollection, ChemicalFormat, ChemicalLinkNode, ChemicalLinkNodeConnection,
    ChemicalRecord, ChemicalSGroup, ChemicalSGroupAttachmentPoint, ChemicalStereoGroup,
    DiagnosticSeverity, FidelityDiagnostic, FidelityPolicy, SdfProperty,
    CHEMICAL_RECORD_SCHEMA_VERSION,
};

#[cfg(feature = "native-layout")]
pub use native_layout::{layout_payload, smiles_to_commands};
