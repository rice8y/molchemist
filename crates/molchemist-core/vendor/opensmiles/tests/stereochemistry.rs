//! Stereochemistry tests
//!
//! Tests for:
//! - Tetrahedral chirality (`@`, `@@`, `@TH1`, `@TH2`)
//! - Extended chirality classes (`@AL`, `@SP`, `@TB`, `@OH`)
//! - Double bond geometry (`/`, `\`)
//! - Interactions with branches and ring closures
//! - Error cases

use opensmiles::{parse, BondType, Chirality, ParserError};

// ============================================================================
// Tetrahedral chirality — basic
// ============================================================================

#[test]
fn parse_at_anticlockwise() {
    let mol = parse("[C@](F)(Cl)(Br)I").expect("Failed to parse");
    assert_eq!(mol.nodes().len(), 5);
    assert_eq!(mol.nodes()[0].chirality(), Some(Chirality::TH1));
}

#[test]
fn parse_at_at_clockwise() {
    let mol = parse("[C@@](F)(Cl)(Br)I").expect("Failed to parse");
    assert_eq!(mol.nodes().len(), 5);
    assert_eq!(mol.nodes()[0].chirality(), Some(Chirality::TH2));
}

#[test]
fn parse_chirality_with_implicit_hydrogen() {
    // L-alanine: N[C@@H](C)C(=O)O
    let mol = parse("N[C@@H](C)C(=O)O").expect("Failed to parse L-alanine");
    assert_eq!(mol.nodes()[1].chirality(), Some(Chirality::TH2));
    assert_eq!(mol.nodes()[1].hydrogens(), 1);
}

#[test]
fn parse_chirality_with_single_h() {
    let mol = parse("[C@H](F)(Cl)Br").expect("Failed to parse");
    assert_eq!(mol.nodes()[0].chirality(), Some(Chirality::TH1));
    assert_eq!(mol.nodes()[0].hydrogens(), 1);
}

#[test]
fn no_chirality_on_organic_atoms() {
    let mol = parse("CCO").expect("Failed to parse");
    for node in mol.nodes() {
        assert_eq!(node.chirality(), None);
    }
}

#[test]
fn no_chirality_on_bracket_without_at() {
    let mol = parse("[CH4]").expect("Failed to parse");
    assert_eq!(mol.nodes()[0].chirality(), None);
}

// ============================================================================
// Tetrahedral chirality — explicit TH
// ============================================================================

#[test]
fn parse_explicit_th1() {
    let mol = parse("[C@TH1](F)(Cl)(Br)I").expect("Failed to parse");
    assert_eq!(mol.nodes()[0].chirality(), Some(Chirality::TH1));
}

#[test]
fn parse_explicit_th2() {
    let mol = parse("[C@TH2](F)(Cl)(Br)I").expect("Failed to parse");
    assert_eq!(mol.nodes()[0].chirality(), Some(Chirality::TH2));
}

// ============================================================================
// Allenal chirality
// ============================================================================

#[test]
fn parse_al1() {
    let mol = parse("[C@AL1](F)(Cl)C").expect("Failed to parse");
    assert_eq!(mol.nodes()[0].chirality(), Some(Chirality::AL1));
}

#[test]
fn parse_al2() {
    let mol = parse("[C@AL2](F)(Cl)C").expect("Failed to parse");
    assert_eq!(mol.nodes()[0].chirality(), Some(Chirality::AL2));
}

// ============================================================================
// Square planar chirality
// ============================================================================

#[test]
fn parse_sp1() {
    let mol = parse("[Pt@SP1](F)(Cl)(Br)I").expect("Failed to parse");
    assert_eq!(mol.nodes()[0].chirality(), Some(Chirality::SP1));
}

#[test]
fn parse_sp2() {
    let mol = parse("[Pt@SP2](F)(Cl)(Br)I").expect("Failed to parse");
    assert_eq!(mol.nodes()[0].chirality(), Some(Chirality::SP2));
}

#[test]
fn parse_sp3() {
    let mol = parse("[Pt@SP3](F)(Cl)(Br)I").expect("Failed to parse");
    assert_eq!(mol.nodes()[0].chirality(), Some(Chirality::SP3));
}

// ============================================================================
// Trigonal bipyramidal chirality
// ============================================================================

#[test]
fn parse_tb1() {
    let mol = parse("[As@TB1](F)(Cl)(Br)(N)S").expect("Failed to parse");
    assert_eq!(mol.nodes()[0].chirality(), Some(Chirality::TB1));
}

#[test]
fn parse_tb20() {
    let mol = parse("[As@TB20](F)(Cl)(Br)(N)S").expect("Failed to parse");
    assert_eq!(mol.nodes()[0].chirality(), Some(Chirality::TB20));
}

#[test]
fn parse_tb_single_digit() {
    let mol = parse("[As@TB5](F)(Cl)(Br)(N)S").expect("Failed to parse");
    assert_eq!(mol.nodes()[0].chirality(), Some(Chirality::TB5));
}

// ============================================================================
// Octahedral chirality
// ============================================================================

#[test]
fn parse_oh1() {
    let mol = parse("[Co@OH1](F)(Cl)(Br)(I)(N)S").expect("Failed to parse");
    assert_eq!(mol.nodes()[0].chirality(), Some(Chirality::OH1));
}

#[test]
fn parse_oh30() {
    let mol = parse("[Co@OH30](F)(Cl)(Br)(I)(N)S").expect("Failed to parse");
    assert_eq!(mol.nodes()[0].chirality(), Some(Chirality::OH30));
}

#[test]
fn parse_oh_single_digit() {
    let mol = parse("[Co@OH5](F)(Cl)(Br)(I)(N)S").expect("Failed to parse");
    assert_eq!(mol.nodes()[0].chirality(), Some(Chirality::OH5));
}

// ============================================================================
// Double bond geometry — directional bonds
// ============================================================================

#[test]
fn parse_trans_difluoroethylene() {
    // F/C=C/F — trans
    let mol = parse("F/C=C/F").expect("Failed to parse");
    assert_eq!(mol.nodes().len(), 4);
    // Bond 0: F-C (Up)
    assert_eq!(mol.bonds()[0].kind(), BondType::Up);
    // Bond 1: C=C (Double)
    assert_eq!(mol.bonds()[1].kind(), BondType::Double);
    // Bond 2: C-F (Up)
    assert_eq!(mol.bonds()[2].kind(), BondType::Up);
}

#[test]
fn parse_cis_difluoroethylene() {
    // F/C=C\F — cis
    let mol = parse(r"F/C=C\F").expect("Failed to parse");
    assert_eq!(mol.nodes().len(), 4);
    assert_eq!(mol.bonds()[0].kind(), BondType::Up);
    assert_eq!(mol.bonds()[1].kind(), BondType::Double);
    assert_eq!(mol.bonds()[2].kind(), BondType::Down);
}

#[test]
fn parse_backslash_backslash_trans() {
    let mol = parse(r"F\C=C\F").expect("Failed to parse");
    assert_eq!(mol.bonds()[0].kind(), BondType::Down);
    assert_eq!(mol.bonds()[1].kind(), BondType::Double);
    assert_eq!(mol.bonds()[2].kind(), BondType::Down);
}

#[test]
fn parse_backslash_slash_cis() {
    let mol = parse(r"F\C=C/F").expect("Failed to parse");
    assert_eq!(mol.bonds()[0].kind(), BondType::Down);
    assert_eq!(mol.bonds()[1].kind(), BondType::Double);
    assert_eq!(mol.bonds()[2].kind(), BondType::Up);
}

#[test]
fn directional_bonds_contribute_single_bond_order() {
    // Fluorine has valence 1, so /F should give 0 implicit H
    let mol = parse("F/C=C/F").expect("Failed to parse");
    // F atoms should have 0 implicit H
    assert_eq!(mol.nodes()[0].hydrogens(), 0); // F
    assert_eq!(mol.nodes()[3].hydrogens(), 0); // F
                                               // C atoms each have 1 double + 1 single (directional) = bond order 3 → 1 implicit H
    assert_eq!(mol.nodes()[1].hydrogens(), 1); // C
    assert_eq!(mol.nodes()[2].hydrogens(), 1); // C
}

// ============================================================================
// Chirality + ring closures
// ============================================================================

#[test]
fn chirality_with_ring_closure() {
    let mol = parse("[C@H]1(F)CCCC1").expect("Failed to parse");
    assert_eq!(mol.nodes()[0].chirality(), Some(Chirality::TH1));
    assert_eq!(mol.nodes()[0].hydrogens(), 1);
}

// ============================================================================
// Chirality + branches
// ============================================================================

#[test]
fn chirality_in_branch() {
    // Molecule with a chiral center inside a branch
    let mol = parse("C([C@H](F)Cl)O").expect("Failed to parse");
    assert_eq!(mol.nodes()[1].chirality(), Some(Chirality::TH1));
}

// ============================================================================
// Real molecules with stereochemistry
// ============================================================================

#[test]
fn parse_l_alanine() {
    let mol = parse("N[C@@H](C)C(=O)O").expect("Failed to parse L-alanine");
    assert_eq!(mol.nodes().len(), 6);
    assert_eq!(mol.nodes()[1].chirality(), Some(Chirality::TH2));
}

#[test]
fn parse_d_alanine() {
    let mol = parse("N[C@H](C)C(=O)O").expect("Failed to parse D-alanine");
    assert_eq!(mol.nodes().len(), 6);
    assert_eq!(mol.nodes()[1].chirality(), Some(Chirality::TH1));
}

// ============================================================================
// Error cases
// ============================================================================

#[test]
fn error_tb0_invalid() {
    let result = parse("[As@TB0](F)(Cl)(Br)(N)S");
    assert!(result.is_err());
    assert!(matches!(
        result.unwrap_err(),
        ParserError::InvalidChiralityClass(..)
    ));
}

#[test]
fn error_tb21_out_of_range() {
    let result = parse("[As@TB21](F)(Cl)(Br)(N)S");
    assert!(result.is_err());
    assert!(matches!(
        result.unwrap_err(),
        ParserError::InvalidChiralityClass(..)
    ));
}

#[test]
fn error_oh0_invalid() {
    let result = parse("[Co@OH0](F)(Cl)(Br)(I)(N)S");
    assert!(result.is_err());
}

#[test]
fn error_oh31_out_of_range() {
    let result = parse("[Co@OH31](F)(Cl)(Br)(I)(N)S");
    assert!(result.is_err());
}

#[test]
fn error_sp0_out_of_range() {
    let result = parse("[Pt@SP0](F)(Cl)(Br)I");
    assert!(result.is_err());
}

#[test]
fn error_sp4_out_of_range() {
    let result = parse("[Pt@SP4](F)(Cl)(Br)I");
    assert!(result.is_err());
}

#[test]
fn error_unknown_chirality_class() {
    let result = parse("[C@XY1](F)(Cl)(Br)I");
    assert!(result.is_err());
}

#[test]
fn error_al0_out_of_range() {
    let result = parse("[C@AL0](F)(Cl)C");
    assert!(result.is_err());
}

#[test]
fn error_al3_out_of_range() {
    let result = parse("[C@AL3](F)(Cl)C");
    assert!(result.is_err());
}

#[test]
fn error_th0_out_of_range() {
    let result = parse("[C@TH0](F)(Cl)(Br)I");
    assert!(result.is_err());
}

#[test]
fn error_th3_out_of_range() {
    let result = parse("[C@TH3](F)(Cl)(Br)I");
    assert!(result.is_err());
}
