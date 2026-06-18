//! Basic SMILES parsing tests
//!
//! These tests verify the basic parser functionality:
//! - Simple atoms (methane)
//! - Linear chains (ethane, ethanol)
//! - Two-letter atoms (chloromethane)

use opensmiles::{parse, AtomSymbol, BondType, OrganicAtom};

#[test]
fn parse_methane() {
    // C = methane: 1 carbon, 0 bonds
    let molecule = parse("C").expect("Failed to parse methane");

    // Check node and bond count
    assert_eq!(molecule.nodes().len(), 1);
    assert_eq!(molecule.bonds().len(), 0);

    // Check the first atom
    let node = &molecule.nodes()[0];
    assert_eq!(*node.atom().element(), AtomSymbol::Organic(OrganicAtom::C));
    assert_eq!(node.atom().charge(), 0);
    assert_eq!(node.atom().isotope(), None);
    assert!(!node.aromatic());
    assert_eq!(node.class(), None);
    assert_eq!(node.hydrogens(), 4);
}

#[test]
fn parse_ethane() {
    // CC = ethane: 2 carbons, 1 simple bond
    let molecule = parse("CC").expect("Failed to parse ethane");

    // Check node and bond count
    assert_eq!(molecule.nodes().len(), 2);
    assert_eq!(molecule.bonds().len(), 1);

    // Check the first carbon (CH3-)
    let node0 = &molecule.nodes()[0];
    assert_eq!(*node0.atom().element(), AtomSymbol::Organic(OrganicAtom::C));
    assert_eq!(node0.atom().charge(), 0);
    assert_eq!(node0.atom().isotope(), None);
    assert!(!node0.aromatic());
    assert_eq!(node0.class(), None);
    assert_eq!(node0.hydrogens(), 3);

    // Check the second carbon (-CH3)
    let node1 = &molecule.nodes()[1];
    assert_eq!(*node1.atom().element(), AtomSymbol::Organic(OrganicAtom::C));
    assert_eq!(node1.atom().charge(), 0);
    assert_eq!(node1.atom().isotope(), None);
    assert!(!node1.aromatic());
    assert_eq!(node1.class(), None);
    assert_eq!(node1.hydrogens(), 3);

    // Check the bond: C(0) - C(1)
    let bond = &molecule.bonds()[0];
    assert_eq!(bond.kind(), BondType::Simple);
    assert_eq!(bond.source(), 0);
    assert_eq!(bond.target(), 1);
}

#[test]
fn parse_ethanol() {
    // CCO = ethanol: 2 carbons + 1 oxygen, 2 simple bonds
    let molecule = parse("CCO").expect("Failed to parse ethanol");

    // Check node and bond count
    assert_eq!(molecule.nodes().len(), 3);
    assert_eq!(molecule.bonds().len(), 2);

    // Check the first carbon (CH3-)
    let node0 = &molecule.nodes()[0];
    assert_eq!(*node0.atom().element(), AtomSymbol::Organic(OrganicAtom::C));
    assert_eq!(node0.atom().charge(), 0);
    assert_eq!(node0.atom().isotope(), None);
    assert!(!node0.aromatic());
    assert_eq!(node0.class(), None);
    assert_eq!(node0.hydrogens(), 3);

    // Check the second carbon (-CH2-)
    let node1 = &molecule.nodes()[1];
    assert_eq!(*node1.atom().element(), AtomSymbol::Organic(OrganicAtom::C));
    assert_eq!(node1.atom().charge(), 0);
    assert_eq!(node1.atom().isotope(), None);
    assert!(!node1.aromatic());
    assert_eq!(node1.class(), None);
    assert_eq!(node1.hydrogens(), 2);

    // Check the oxygen (-OH)
    let node2 = &molecule.nodes()[2];
    assert_eq!(*node2.atom().element(), AtomSymbol::Organic(OrganicAtom::O));
    assert_eq!(node2.atom().charge(), 0);
    assert_eq!(node2.atom().isotope(), None);
    assert!(!node2.aromatic());
    assert_eq!(node2.class(), None);
    assert_eq!(node2.hydrogens(), 1);

    // Check the bonds: C(0) - C(1) - O(2)
    let bond0 = &molecule.bonds()[0];
    assert_eq!(bond0.kind(), BondType::Simple);
    assert_eq!(bond0.source(), 0);
    assert_eq!(bond0.target(), 1);

    let bond1 = &molecule.bonds()[1];
    assert_eq!(bond1.kind(), BondType::Simple);
    assert_eq!(bond1.source(), 1);
    assert_eq!(bond1.target(), 2);
}

#[test]
fn parse_chloromethane() {
    // CCl = chloromethane: tests two-letter atoms
    let molecule = parse("CCl").expect("Failed to parse chloromethane");

    // Check node and bond count
    assert_eq!(molecule.nodes().len(), 2);
    assert_eq!(molecule.bonds().len(), 1);

    // Check the carbon (CH3-)
    let node0 = &molecule.nodes()[0];
    assert_eq!(*node0.atom().element(), AtomSymbol::Organic(OrganicAtom::C));
    assert_eq!(node0.atom().charge(), 0);
    assert_eq!(node0.atom().isotope(), None);
    assert!(!node0.aromatic());
    assert_eq!(node0.class(), None);
    assert_eq!(node0.hydrogens(), 3);

    // Check the chlorine (-Cl)
    let node1 = &molecule.nodes()[1];
    assert_eq!(
        *node1.atom().element(),
        AtomSymbol::Organic(OrganicAtom::Cl)
    );
    assert_eq!(node1.atom().charge(), 0);
    assert_eq!(node1.atom().isotope(), None);
    assert!(!node1.aromatic());
    assert_eq!(node1.class(), None);
    assert_eq!(node1.hydrogens(), 0);

    // Check the bond: C(0) - Cl(1)
    let bond = &molecule.bonds()[0];
    assert_eq!(bond.kind(), BondType::Simple);
    assert_eq!(bond.source(), 0);
    assert_eq!(bond.target(), 1);
}
