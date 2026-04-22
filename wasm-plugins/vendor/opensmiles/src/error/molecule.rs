//! Molecule-related errors.

use thiserror::Error;

use super::{AtomError, NodeError};
use crate::NodeIndex;

/// Errors that can occur when building a molecule.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum MoleculeError {
    /// Error from a node.
    #[error(transparent)]
    NodeError(#[from] NodeError),

    /// Error from an atom.
    #[error(transparent)]
    AtomError(#[from] AtomError),

    /// Aromatic ring violates Hückel's rule (4n+2 pi electrons required).
    #[error("aromatic ring {ring:?} has {pi_electrons} pi electrons (Huckel: 4n+2 required)")]
    HuckelViolation {
        ring: Vec<NodeIndex>,
        pi_electrons: u8,
    },
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn node_error_conversion() {
        let node_err = NodeError::UndefinedHydrogen;
        let mol_err: MoleculeError = node_err.into();

        assert!(matches!(mol_err, MoleculeError::NodeError(_)));
        assert_eq!(mol_err.to_string(), "undefined hydrogen count");
    }

    #[test]
    fn atom_error_conversion() {
        let atom_err = AtomError::UnknownElement("Zz".to_string());
        let mol_err: MoleculeError = atom_err.into();

        assert!(matches!(mol_err, MoleculeError::AtomError(_)));
        assert_eq!(mol_err.to_string(), "unknown element: 'Zz'");
    }
}
