use thiserror::Error;

/// Errors that can occur when creating or manipulating an atom.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum BondError {
    /// The specified bond is not recognized.
    #[error("unknown bond: '{0}'")]
    UnknownBond(char),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn error_messages_are_descriptive() {
        assert_eq!(BondError::UnknownBond('X').to_string(), "unknown bond: 'X'");
    }
}
