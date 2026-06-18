//! Atom-related errors.

use thiserror::Error;

/// Errors that can occur when creating or manipulating an atom.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum AtomError {
    /// The specified charge is out of range (must be between -15 and +15).
    #[error("invalid charge: {0} (must be between -15 and +15)")]
    InvalidCharge(i8),

    /// The specified isotope is too large (maximum 999).
    #[error("invalid isotope: {0} (maximum 999)")]
    InvalidIsotope(u16),

    /// Bond order is required but was not provided.
    #[error("missing bond order to compute implicit hydrogens")]
    MissingBondOrder,

    /// The specified element is not recognized.
    #[error("unknown element: '{0}'")]
    UnknownElement(String),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn error_messages_are_descriptive() {
        assert_eq!(
            AtomError::InvalidCharge(-20).to_string(),
            "invalid charge: -20 (must be between -15 and +15)"
        );

        assert_eq!(
            AtomError::InvalidIsotope(1000).to_string(),
            "invalid isotope: 1000 (maximum 999)"
        );

        assert_eq!(
            AtomError::MissingBondOrder.to_string(),
            "missing bond order to compute implicit hydrogens"
        );

        assert_eq!(
            AtomError::UnknownElement("Xy".to_string()).to_string(),
            "unknown element: 'Xy'"
        );
    }
}
