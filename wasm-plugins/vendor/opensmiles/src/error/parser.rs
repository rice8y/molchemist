//! SMILES parsing errors.

use thiserror::Error;

use super::{BondError, MoleculeError, NodeError};
use crate::NodeIndex;

/// Errors that can occur when parsing a SMILES string.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum ParserError {
    /// Feature not yet implemented.
    #[error("feature not yet implemented")]
    NotYetImplemented,

    /// The molecule contains too many nodes (maximum 4294967295).
    #[error("too many nodes in molecule (maximum 4294967295)")]
    TooManyNodes,

    #[error("At least one node is necessary before creating a bond")]
    NoAtomToBond,

    /// Unexpected character in SMILES string.
    #[error("unexpected character '{0}' at position {1}")]
    UnexpectedCharacter(char, usize),

    /// Unexpected end of input.
    #[error("unexpected end of input, expected: {0}")]
    UnexpectedEndOfInput(String),

    /// Bracket Atom must have an element
    #[error("brackets atom must have an element")]
    MissingElementInBracketAtom,

    #[error("charge in bracket atom must have a sign")]
    ChargeWithoutSign,

    /// Missing closing parenthesis.
    #[error("missing closing parenthesis ')'")]
    UnclosedParenthesis,

    /// Missing opening parenthesis.
    #[error("missing opening parenthesis '('")]
    UnopenedParenthesis,

    /// Empty branch.
    #[error("empty branch detected")]
    EmptyBranch,

    /// Unclosed ring.
    #[error("unclosed ring(s): {0:?}")]
    UnclosedRing(Vec<u8>),

    /// Mismatched bond types for ring closure.
    #[error("mismatched bond types for ring {0}")]
    MismatchedRingBond(u8),

    /// Bond without preceding atom.
    #[error("bond without preceding atom")]
    BondWithoutPrecedingAtom,

    /// Bond without following atom.
    #[error("bond without following atom")]
    BondWithoutFollowingAtom,

    /// Hydrogens having hydrogens count is illegal
    #[error("hydrogens cannot have hydrogens count")]
    HydrogenWithHydrogenCount,

    #[error("charge should be between -15 and +15 {0}")]
    ChargeOutOfRange(String),

    #[error("hydrogen cannot be greater than 9 {0}")]
    HydrogenOutOfRange(String),

    #[error("invalid chirality class: {0} at position {1}")]
    InvalidChiralityClass(String, usize),

    #[error("invalid chirality specification: {0} at position {1}")]
    InvalidChiralitySpec(String, usize),

    /// Atom bonded to itself (e.g., C11).
    #[error("atom cannot be bonded to itself (ring {0})")]
    SelfBond(u8),

    /// Duplicate bond between the same pair of atoms (e.g., C12CCCCC12).
    #[error("duplicate bond between atoms {0} and {1}")]
    DuplicateBond(NodeIndex, NodeIndex),

    /// Error from molecule construction.
    #[error(transparent)]
    MoleculeError(#[from] MoleculeError),

    /// Error from a node.
    #[error(transparent)]
    NodeError(#[from] NodeError),

    /// Error from a bond.
    #[error(transparent)]
    BondError(#[from] BondError),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn error_messages_are_descriptive() {
        assert_eq!(
            ParserError::NotYetImplemented.to_string(),
            "feature not yet implemented"
        );

        assert_eq!(
            ParserError::TooManyNodes.to_string(),
            "too many nodes in molecule (maximum 4294967295)"
        );

        assert_eq!(
            ParserError::UnexpectedCharacter('$', 5).to_string(),
            "unexpected character '$' at position 5"
        );

        assert_eq!(
            ParserError::UnexpectedEndOfInput("element".to_string()).to_string(),
            "unexpected end of input, expected: element"
        );

        assert_eq!(
            ParserError::UnclosedParenthesis.to_string(),
            "missing closing parenthesis ')'"
        );

        assert_eq!(
            ParserError::UnclosedRing(vec!(1, 2, 5)).to_string(),
            "unclosed ring(s): [1, 2, 5]"
        );

        assert_eq!(
            ParserError::MismatchedRingBond(2).to_string(),
            "mismatched bond types for ring 2"
        );

        assert_eq!(
            ParserError::BondWithoutPrecedingAtom.to_string(),
            "bond without preceding atom"
        );

        assert_eq!(
            ParserError::BondWithoutFollowingAtom.to_string(),
            "bond without following atom"
        );
    }

    #[test]
    fn molecule_error_conversion() {
        let mol_err = MoleculeError::NodeError(super::super::NodeError::UndefinedHydrogen);
        let parser_err: ParserError = mol_err.into();

        assert!(matches!(parser_err, ParserError::MoleculeError(_)));
    }

    #[test]
    fn node_error_conversion() {
        let node_err = super::super::NodeError::InvalidHydrogen(99);
        let parser_err: ParserError = node_err.into();

        assert!(matches!(parser_err, ParserError::NodeError(_)));
        assert_eq!(parser_err.to_string(), "invalid hydrogen count: 99");
    }

    #[test]
    fn bond_error_conversion() {
        let bond_err = super::super::BondError::UnknownBond('c');
        let parser_err: ParserError = bond_err.into();

        assert!(matches!(parser_err, ParserError::BondError(_)));
        assert_eq!(parser_err.to_string(), "unknown bond: 'c'");
    }
}
