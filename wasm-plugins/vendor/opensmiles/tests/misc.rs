//! Miscellaneous tests
//!
//! These tests cover:
//! - Wildcards (`*`)
//! - Complex real molecules
//! - Edge cases

use opensmiles::{parse, AtomSymbol};

// ============================================================================
// Wildcards
// ============================================================================

#[test]
fn parse_wildcard_atom() {
    // * = wildcard atom (any atom)
    let molecule = parse("*").expect("Failed to parse wildcard");

    assert_eq!(molecule.nodes().len(), 1);
    assert_eq!(*molecule.nodes()[0].atom().element(), AtomSymbol::Wildcard);
}

#[test]
fn parse_wildcard_in_chain() {
    // C*C = chain with wildcard in the middle
    let molecule = parse("C*C").expect("Failed to parse C*C");

    assert_eq!(molecule.nodes().len(), 3);
    assert_eq!(molecule.bonds().len(), 2);

    assert_eq!(*molecule.nodes()[1].atom().element(), AtomSymbol::Wildcard);
}

#[test]
fn parse_wildcard_in_ring() {
    // C1*CC1 = ring with wildcard
    let molecule = parse("C1*CC1").expect("Failed to parse ring with wildcard");

    assert_eq!(molecule.nodes().len(), 4);
    assert_eq!(molecule.bonds().len(), 4);
}

// ============================================================================
// Complex real molecules
// ============================================================================

#[test]
fn parse_aspirin() {
    // CC(=O)Oc1ccccc1C(=O)O = aspirin
    let molecule = parse("CC(=O)Oc1ccccc1C(=O)O").expect("Failed to parse aspirin");

    assert_eq!(molecule.nodes().len(), 13);
}

#[test]
fn parse_caffeine() {
    // Cn1cnc2c1c(=O)n(c(=O)n2C)C = caffeine
    let molecule = parse("Cn1cnc2c1c(=O)n(c(=O)n2C)C").expect("Failed to parse caffeine");

    // Caffeine has 14 heavy atoms (without hydrogens)
    assert_eq!(molecule.nodes().len(), 14);
}

#[test]
fn parse_glucose() {
    // Simplified version without stereochemistry: OCC1OC(O)C(O)C(O)C1O
    let molecule = parse("OCC1OC(O)C(O)C(O)C1O").expect("Failed to parse glucose");

    assert_eq!(molecule.nodes().len(), 12); // 6 C + 6 O
}

#[test]
fn parse_ethanol_explicit() {
    // [CH3][CH2][OH] = ethanol with everything explicit
    let molecule = parse("[CH3][CH2][OH]").expect("Failed to parse explicit ethanol");

    assert_eq!(molecule.nodes().len(), 3);
    assert_eq!(molecule.bonds().len(), 2);

    assert_eq!(molecule.nodes()[0].hydrogens(), 3);
    assert_eq!(molecule.nodes()[1].hydrogens(), 2);
    assert_eq!(molecule.nodes()[2].hydrogens(), 1);
}

#[test]
fn parse_acetone() {
    // CC(=O)C = acetone
    let molecule = parse("CC(=O)C").expect("Failed to parse acetone");

    assert_eq!(molecule.nodes().len(), 4);
    assert_eq!(molecule.bonds().len(), 3);
}

#[test]
fn parse_acetic_acid() {
    // CC(=O)O = acetic acid
    let molecule = parse("CC(=O)O").expect("Failed to parse acetic acid");

    assert_eq!(molecule.nodes().len(), 4);
}

#[test]
fn parse_benzaldehyde() {
    // c1ccccc1C=O = benzaldehyde
    let molecule = parse("c1ccccc1C=O").expect("Failed to parse benzaldehyde");

    assert_eq!(molecule.nodes().len(), 8);
}

#[test]
fn parse_aniline() {
    // Nc1ccccc1 = aniline
    let molecule = parse("Nc1ccccc1").expect("Failed to parse aniline");

    assert_eq!(molecule.nodes().len(), 7);
}

#[test]
fn parse_phenol() {
    // Oc1ccccc1 = phenol
    let molecule = parse("Oc1ccccc1").expect("Failed to parse phenol");

    assert_eq!(molecule.nodes().len(), 7);
}

// ============================================================================
// Edge cases and validation
// ============================================================================

#[test]
fn parse_empty_string() {
    // Empty string should produce an error or an empty molecule
    let result = parse("");
    // Depending on implementation, this can be Ok with 0 atoms or Err
    if let Ok(molecule) = result {
        assert_eq!(molecule.nodes().len(), 0);
    }
}

#[test]
fn parse_single_bond_at_start() {
    // -C should be invalid (bond without preceding atom)
    let result = parse("-C");
    assert!(result.is_err(), "Bond at start should be invalid");
}

#[test]
fn parse_double_bond_at_end() {
    // C= should be invalid (bond without following atom)
    let result = parse("C=");
    assert!(result.is_err(), "Bond at end should be invalid");
}

#[test]
fn parse_unclosed_bracket() {
    // [C should be invalid
    let result = parse("[C");
    assert!(result.is_err(), "Unclosed bracket should be invalid");
}

#[test]
fn parse_unclosed_ring() {
    // C1CC should be invalid (unclosed ring)
    let result = parse("C1CC");
    assert!(result.is_err(), "Unclosed ring should be invalid");
}

#[test]
fn parse_mismatched_ring_bonds() {
    // C=1CC1 vs C1CC=1 - bond types at closures must match
    // C1CC=1 is valid (simple then double = double)
    // C=1CC1 is valid (double then simple = double)
    // C=1CC-1 should be invalid (double vs explicit simple)
    let result = parse("C=1CC-1");
    assert!(
        result.is_err(),
        "Mismatched ring bond types should be invalid"
    );
}
