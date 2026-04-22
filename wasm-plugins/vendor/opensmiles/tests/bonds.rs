//! Explicit bond tests
//!
//! These tests verify the parsing of different bond types:
//! - Simple (`-`)
//! - Double (`=`)
//! - Triple (`#`)
//! - Quadruple (`$`)
//! - Aromatic (`:`)

use opensmiles::{parse, BondType};

#[test]
fn parse_explicit_single_bond() {
    // C-C = ethane with explicit single bond
    let molecule = parse("C-C").expect("Failed to parse C-C");

    assert_eq!(molecule.nodes().len(), 2);
    assert_eq!(molecule.bonds().len(), 1);

    let bond = &molecule.bonds()[0];
    assert_eq!(bond.kind(), BondType::Simple);
    assert_eq!(bond.source(), 0);
    assert_eq!(bond.target(), 1);

    // Check implicit hydrogens (CH3-CH3)
    assert_eq!(molecule.nodes()[0].hydrogens(), 3);
    assert_eq!(molecule.nodes()[1].hydrogens(), 3);
}

#[test]
fn parse_double_bond() {
    // C=C = ethene (ethylene): double bond
    let molecule = parse("C=C").expect("Failed to parse ethene");

    assert_eq!(molecule.nodes().len(), 2);
    assert_eq!(molecule.bonds().len(), 1);

    let bond = &molecule.bonds()[0];
    assert_eq!(bond.kind(), BondType::Double);
    assert_eq!(bond.source(), 0);
    assert_eq!(bond.target(), 1);

    // CH2=CH2: each carbon has 2 hydrogens
    assert_eq!(molecule.nodes()[0].hydrogens(), 2);
    assert_eq!(molecule.nodes()[1].hydrogens(), 2);
}

#[test]
fn parse_triple_bond() {
    // C#C = ethyne (acetylene): triple bond
    let molecule = parse("C#C").expect("Failed to parse ethyne");

    assert_eq!(molecule.nodes().len(), 2);
    assert_eq!(molecule.bonds().len(), 1);

    let bond = &molecule.bonds()[0];
    assert_eq!(bond.kind(), BondType::Triple);
    assert_eq!(bond.source(), 0);
    assert_eq!(bond.target(), 1);

    // HC≡CH: each carbon has 1 hydrogen
    assert_eq!(molecule.nodes()[0].hydrogens(), 1);
    assert_eq!(molecule.nodes()[1].hydrogens(), 1);
}

#[test]
fn parse_quadruple_bond() {
    // C$C = quadruple bond (rare, used for some metal complexes)
    let molecule = parse("C$C").expect("Failed to parse quadruple bond");

    assert_eq!(molecule.nodes().len(), 2);
    assert_eq!(molecule.bonds().len(), 1);

    let bond = &molecule.bonds()[0];
    assert_eq!(bond.kind(), BondType::Quadruple);
    assert_eq!(bond.source(), 0);
    assert_eq!(bond.target(), 1);

    // No hydrogens with a quadruple bond
    assert_eq!(molecule.nodes()[0].hydrogens(), 0);
    assert_eq!(molecule.nodes()[1].hydrogens(), 0);
}

#[test]
fn parse_mixed_bonds() {
    // C=C-C#N = acrylonitrile: mix of bonds
    let molecule = parse("C=C-C#N").expect("Failed to parse acrylonitrile");

    assert_eq!(molecule.nodes().len(), 4);
    assert_eq!(molecule.bonds().len(), 3);

    // Double bond C=C
    assert_eq!(molecule.bonds()[0].kind(), BondType::Double);
    assert_eq!(molecule.bonds()[0].source(), 0);
    assert_eq!(molecule.bonds()[0].target(), 1);

    // Single bond C-C
    assert_eq!(molecule.bonds()[1].kind(), BondType::Simple);
    assert_eq!(molecule.bonds()[1].source(), 1);
    assert_eq!(molecule.bonds()[1].target(), 2);

    // Triple bond C#N
    assert_eq!(molecule.bonds()[2].kind(), BondType::Triple);
    assert_eq!(molecule.bonds()[2].source(), 2);
    assert_eq!(molecule.bonds()[2].target(), 3);
}

#[test]
fn parse_explicit_then_implicit_bond() {
    // C=CC = propene
    let molecule = parse("C=CC").expect("Failed to parse propene");

    assert_eq!(molecule.nodes().len(), 3);
    assert_eq!(molecule.bonds().len(), 2);

    // Double bond C=C
    assert_eq!(molecule.bonds()[0].kind(), BondType::Double);
    assert_eq!(molecule.bonds()[0].source(), 0);
    assert_eq!(molecule.bonds()[0].target(), 1);

    // Single bond C-C
    assert_eq!(molecule.bonds()[1].kind(), BondType::Simple);
    assert_eq!(molecule.bonds()[1].source(), 1);
    assert_eq!(molecule.bonds()[1].target(), 2);

    assert_eq!(molecule.nodes()[0].hydrogens(), 2);
    assert_eq!(molecule.nodes()[1].hydrogens(), 1);
    assert_eq!(molecule.nodes()[2].hydrogens(), 3);
}

#[test]
fn parse_aromatic_bond() {
    // c:c = explicit aromatic bond
    let molecule = parse("c:c").expect("Failed to parse aromatic bond");

    assert_eq!(molecule.nodes().len(), 2);
    assert_eq!(molecule.bonds().len(), 1);

    let bond = &molecule.bonds()[0];
    assert_eq!(bond.kind(), BondType::Aromatic);

    // Both atoms must be aromatic
    assert!(molecule.nodes()[0].aromatic());
    assert!(molecule.nodes()[1].aromatic());
}
