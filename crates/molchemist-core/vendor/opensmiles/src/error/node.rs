//! Node-related errors (atoms in the molecular graph).

use thiserror::Error;

use crate::AtomSymbol;

use super::AtomError;

/// Errors that can occur when creating or manipulating a node.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum NodeError {
    /// The specified hydrogen count is invalid.
    #[error("invalid hydrogen count: {0}")]
    InvalidHydrogen(u8),

    /// The hydrogen count was not defined when required.
    #[error("undefined hydrogen count")]
    UndefinedHydrogen,

    #[error("Invalid element {0} cannot be aromatic")]
    InvalidAromaticElement(AtomSymbol),

    /// The specified atom class is invalid.
    #[error("invalid atom class: {0}")]
    InvalidClass(u16),

    /// Aromaticity was not defined when required.
    #[error("undefined aromaticity")]
    UndefinedAromatic,

    /// Bond order is mandatory for organic atoms.
    #[error("bond order is mandatory for organic atoms")]
    BondOrderMandatoryForOrganicAtom,

    /// Error from the underlying atom.
    #[error(transparent)]
    AtomError(#[from] AtomError),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn error_messages_are_descriptive() {
        assert_eq!(
            NodeError::InvalidHydrogen(10).to_string(),
            "invalid hydrogen count: 10"
        );

        assert_eq!(
            NodeError::UndefinedHydrogen.to_string(),
            "undefined hydrogen count"
        );

        assert_eq!(
            NodeError::InvalidClass(65535).to_string(),
            "invalid atom class: 65535"
        );

        assert_eq!(
            NodeError::UndefinedAromatic.to_string(),
            "undefined aromaticity"
        );

        assert_eq!(
            NodeError::BondOrderMandatoryForOrganicAtom.to_string(),
            "bond order is mandatory for organic atoms"
        );
    }

    #[test]
    fn atom_error_conversion() {
        let atom_err = AtomError::InvalidCharge(20);
        let node_err: NodeError = atom_err.into();

        assert!(matches!(node_err, NodeError::AtomError(_)));
        assert_eq!(
            node_err.to_string(),
            "invalid charge: 20 (must be between -15 and +15)"
        );
    }
}
