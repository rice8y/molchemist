//! Parallel parsing utilities for SMILES strings.
//!
//! This module provides parallel parsing capabilities using Rayon.
//! Enable with the `parallel` feature flag.
//!
//! ## Usage
//!
//! ```rust
//! use opensmiles::parse_batch;
//!
//! let smiles = vec!["CCO", "c1ccccc1", "CC(=O)O"];
//! let results = parse_batch(&smiles);
//! assert_eq!(results.len(), 3);
//! assert!(results.iter().all(|r| r.is_ok()));
//! ```

use rayon::prelude::*;

use crate::error::ParserError;
use crate::parser::parse;
use crate::Molecule;

/// Parse multiple SMILES strings in parallel.
///
/// This is the most effective parallelization strategy for SMILES parsing,
/// as each molecule is completely independent and can be parsed on a separate thread.
///
/// # Arguments
///
/// * `inputs` - A slice of SMILES strings to parse
///
/// # Returns
///
/// A vector of Results, one for each input string, in the same order as the input.
///
/// # Example
///
/// ```rust
/// use opensmiles::parse_batch;
///
/// let smiles = vec!["CCO", "c1ccccc1", "CC(=O)O", "invalid["];
/// let results = parse_batch(&smiles);
///
/// assert!(results[0].is_ok()); // ethanol
/// assert!(results[1].is_ok()); // benzene
/// assert!(results[2].is_ok()); // acetic acid
/// assert!(results[3].is_err()); // invalid
/// ```
pub fn parse_batch(inputs: &[&str]) -> Vec<Result<Molecule, ParserError>> {
    inputs.par_iter().map(|s| parse(s)).collect()
}

/// Parse multiple SMILES strings in parallel, returning only successful results.
///
/// This is useful when you want to process a large dataset and filter out invalid entries.
///
/// # Arguments
///
/// * `inputs` - A slice of SMILES strings to parse
///
/// # Returns
///
/// A vector of successfully parsed Molecules. Invalid SMILES are silently skipped.
///
/// # Example
///
/// ```rust
/// use opensmiles::parse_batch_ok;
///
/// let smiles = vec!["CCO", "invalid[", "c1ccccc1"];
/// let molecules = parse_batch_ok(&smiles);
///
/// assert_eq!(molecules.len(), 2); // only valid ones
/// ```
pub fn parse_batch_ok(inputs: &[&str]) -> Vec<Molecule> {
    inputs.par_iter().filter_map(|s| parse(s).ok()).collect()
}

/// Parse multiple SMILES strings in parallel with their indices.
///
/// This preserves the association between input index and result,
/// useful for tracking which inputs failed.
///
/// # Arguments
///
/// * `inputs` - A slice of SMILES strings to parse
///
/// # Returns
///
/// A vector of tuples (index, Result) for each input.
pub fn parse_batch_indexed(inputs: &[&str]) -> Vec<(usize, Result<Molecule, ParserError>)> {
    inputs
        .par_iter()
        .enumerate()
        .map(|(i, s)| (i, parse(s)))
        .collect()
}

/// Statistics from a batch parse operation.
#[derive(Debug, Clone, Default)]
pub struct BatchParseStats {
    /// Number of successfully parsed molecules
    pub success_count: usize,
    /// Number of failed parses
    pub error_count: usize,
    /// Total number of inputs
    pub total_count: usize,
}

impl BatchParseStats {
    /// Success rate as a percentage (0.0 to 100.0)
    pub fn success_rate(&self) -> f64 {
        if self.total_count == 0 {
            0.0
        } else {
            (self.success_count as f64 / self.total_count as f64) * 100.0
        }
    }
}

/// Parse multiple SMILES strings in parallel and return statistics.
///
/// # Arguments
///
/// * `inputs` - A slice of SMILES strings to parse
///
/// # Returns
///
/// A tuple of (successful molecules, errors with their indices, statistics)
pub fn parse_batch_with_stats(
    inputs: &[&str],
) -> (Vec<Molecule>, Vec<(usize, ParserError)>, BatchParseStats) {
    let results: Vec<(usize, Result<Molecule, ParserError>)> = inputs
        .par_iter()
        .enumerate()
        .map(|(i, s)| (i, parse(s)))
        .collect();

    let mut molecules = Vec::new();
    let mut errors = Vec::new();

    for (i, result) in results {
        match result {
            Ok(mol) => molecules.push(mol),
            Err(e) => errors.push((i, e)),
        }
    }

    let stats = BatchParseStats {
        success_count: molecules.len(),
        error_count: errors.len(),
        total_count: inputs.len(),
    };

    (molecules, errors, stats)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_batch() {
        let inputs = vec!["C", "CC", "CCC", "CCCC"];
        let results = parse_batch(&inputs);

        assert_eq!(results.len(), 4);
        for result in &results {
            assert!(result.is_ok());
        }
    }

    #[test]
    fn test_parse_batch_with_errors() {
        let inputs = vec!["C", "invalid[", "CC"];
        let results = parse_batch(&inputs);

        assert_eq!(results.len(), 3);
        assert!(results[0].is_ok());
        assert!(results[1].is_err());
        assert!(results[2].is_ok());
    }

    #[test]
    fn test_parse_batch_ok() {
        let inputs = vec!["C", "invalid[", "CC", "also-bad"];
        let molecules = parse_batch_ok(&inputs);

        assert_eq!(molecules.len(), 2);
    }

    #[test]
    fn test_parse_batch_with_stats() {
        let inputs = vec!["C", "CC", "bad[", "CCC"];
        let (molecules, errors, stats) = parse_batch_with_stats(&inputs);

        assert_eq!(molecules.len(), 3);
        assert_eq!(errors.len(), 1);
        assert_eq!(errors[0].0, 2); // index of "bad["
        assert_eq!(stats.success_count, 3);
        assert_eq!(stats.error_count, 1);
        assert_eq!(stats.total_count, 4);
        assert!((stats.success_rate() - 75.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_parse_complex_molecules() {
        // Test with more complex real-world SMILES
        let inputs = vec![
            "c1ccccc1",       // benzene
            "CC(=O)O",        // acetic acid
            "CCO",            // ethanol
            "C1CC1",          // cyclopropane
            "CC(C)(C)C",      // neopentane
            "c1ccc2ccccc2c1", // naphthalene
            "CC(C)CC(C)(C)C", // 2,2,4-trimethylpentane
        ];

        let results = parse_batch(&inputs);
        for (i, result) in results.iter().enumerate() {
            assert!(result.is_ok(), "Failed to parse index {}: {:?}", i, result);
        }
    }

    #[test]
    fn test_parse_batch_indexed() {
        let inputs = vec!["C", "bad[", "CC", "also-bad", "CCC"];
        let results = parse_batch_indexed(&inputs);

        assert_eq!(results.len(), 5);

        // Check that indices are preserved
        let (idx0, res0) = &results[0];
        assert_eq!(*idx0, 0);
        assert!(res0.is_ok());

        let (idx1, res1) = &results[1];
        assert_eq!(*idx1, 1);
        assert!(res1.is_err());

        let (idx2, res2) = &results[2];
        assert_eq!(*idx2, 2);
        assert!(res2.is_ok());

        let (idx3, res3) = &results[3];
        assert_eq!(*idx3, 3);
        assert!(res3.is_err());

        let (idx4, res4) = &results[4];
        assert_eq!(*idx4, 4);
        assert!(res4.is_ok());
    }

    #[test]
    fn test_parse_batch_empty_input() {
        let inputs: Vec<&str> = vec![];
        let results = parse_batch(&inputs);
        assert!(results.is_empty());

        let molecules = parse_batch_ok(&inputs);
        assert!(molecules.is_empty());

        let (molecules, errors, stats) = parse_batch_with_stats(&inputs);
        assert!(molecules.is_empty());
        assert!(errors.is_empty());
        assert_eq!(stats.total_count, 0);
        assert_eq!(stats.success_rate(), 0.0);
    }

    #[test]
    fn test_parse_batch_all_invalid() {
        let inputs = vec!["bad[", "also-bad", "[invalid"];
        let (molecules, errors, stats) = parse_batch_with_stats(&inputs);

        assert!(molecules.is_empty());
        assert_eq!(errors.len(), 3);
        assert_eq!(stats.success_count, 0);
        assert_eq!(stats.error_count, 3);
        assert_eq!(stats.success_rate(), 0.0);
    }

    #[test]
    fn test_batch_parse_stats_success_rate() {
        let stats = BatchParseStats {
            success_count: 3,
            error_count: 1,
            total_count: 4,
        };
        assert!((stats.success_rate() - 75.0).abs() < f64::EPSILON);

        let empty_stats = BatchParseStats::default();
        assert_eq!(empty_stats.success_rate(), 0.0);
    }
}
