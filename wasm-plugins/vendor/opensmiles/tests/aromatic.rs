//! Aromatic atom tests
//!
//! These tests verify the parsing of aromatic atoms (lowercase):
//! - Aromatic carbons (`c`)
//! - Aromatic heteroatoms (`n`, `o`, `s`)
//! - Common aromatic rings (pyridine, furan, thiophene)

use opensmiles::{parse, AtomSymbol, BondType, OrganicAtom};

#[test]
fn parse_aromatic_carbon() {
    // c = single aromatic carbon
    let molecule = parse("c").expect("Failed to parse aromatic carbon");

    assert_eq!(molecule.nodes().len(), 1);
    assert!(molecule.nodes()[0].aromatic());
    assert_eq!(
        *molecule.nodes()[0].atom().element(),
        AtomSymbol::Organic(OrganicAtom::C)
    );
}

#[test]
fn parse_aromatic_nitrogen() {
    // n = aromatic nitrogen
    let molecule = parse("n").expect("Failed to parse aromatic nitrogen");

    assert!(molecule.nodes()[0].aromatic());
    assert_eq!(
        *molecule.nodes()[0].atom().element(),
        AtomSymbol::Organic(OrganicAtom::N)
    );
}

#[test]
fn parse_aromatic_oxygen() {
    // o = aromatic oxygen
    let molecule = parse("o").expect("Failed to parse aromatic oxygen");

    assert!(molecule.nodes()[0].aromatic());
    assert_eq!(
        *molecule.nodes()[0].atom().element(),
        AtomSymbol::Organic(OrganicAtom::O)
    );
}

#[test]
fn parse_aromatic_sulfur() {
    // s = aromatic sulfur
    let molecule = parse("s").expect("Failed to parse aromatic sulfur");

    assert!(molecule.nodes()[0].aromatic());
    assert_eq!(
        *molecule.nodes()[0].atom().element(),
        AtomSymbol::Organic(OrganicAtom::S)
    );
}

#[test]
fn parse_benzene() {
    // c1ccccc1 = benzene: 6 aromatic carbons in a ring
    let molecule = parse("c1ccccc1").expect("Failed to parse benzene");

    // Check node and bond count
    assert_eq!(molecule.nodes().len(), 6);
    assert_eq!(molecule.bonds().len(), 6); // closed ring

    // All atoms must be aromatic with the same properties
    for node in molecule.nodes() {
        assert_eq!(*node.atom().element(), AtomSymbol::Organic(OrganicAtom::C));
        assert_eq!(node.atom().charge(), 0);
        assert_eq!(node.atom().isotope(), None);
        assert!(node.aromatic());
        assert_eq!(node.class(), None);
        assert_eq!(node.hydrogens(), 1); // each aromatic carbon has 1 hydrogen
    }

    // Check that all bonds are aromatic
    for bond in molecule.bonds() {
        assert_eq!(bond.kind(), BondType::Aromatic);
    }
}

#[test]
fn parse_pyridine() {
    // c1ccncc1 = pyridine (aromatic ring with nitrogen)
    let molecule = parse("c1ccncc1").expect("Failed to parse pyridine");

    assert_eq!(molecule.nodes().len(), 6);
    assert_eq!(molecule.bonds().len(), 6);

    // 5 aromatic carbons + 1 aromatic nitrogen
    let carbons: Vec<_> = molecule
        .nodes()
        .iter()
        .filter(|n| *n.atom().element() == AtomSymbol::Organic(OrganicAtom::C))
        .collect();
    let nitrogens: Vec<_> = molecule
        .nodes()
        .iter()
        .filter(|n| *n.atom().element() == AtomSymbol::Organic(OrganicAtom::N))
        .collect();

    assert_eq!(carbons.len(), 5);
    assert_eq!(nitrogens.len(), 1);

    // All must be aromatic
    for node in molecule.nodes() {
        assert!(node.aromatic());
    }
}

#[test]
fn parse_furan() {
    // c1ccoc1 = furan (aromatic ring with oxygen)
    let molecule = parse("c1ccoc1").expect("Failed to parse furan");

    assert_eq!(molecule.nodes().len(), 5);
    assert_eq!(molecule.bonds().len(), 5);

    // Check that there is an aromatic oxygen
    let oxygen = molecule
        .nodes()
        .iter()
        .find(|n| *n.atom().element() == AtomSymbol::Organic(OrganicAtom::O));
    assert!(oxygen.is_some());
    assert!(oxygen.unwrap().aromatic());
}

#[test]
fn parse_thiophene() {
    // c1ccsc1 = thiophene (aromatic ring with sulfur)
    let molecule = parse("c1ccsc1").expect("Failed to parse thiophene");

    assert_eq!(molecule.nodes().len(), 5);

    // Check that there is an aromatic sulfur
    let sulfur = molecule
        .nodes()
        .iter()
        .find(|n| *n.atom().element() == AtomSymbol::Organic(OrganicAtom::S));
    assert!(sulfur.is_some());
    assert!(sulfur.unwrap().aromatic());
}

#[test]
fn parse_pyrrole() {
    // c1cc[nH]c1 = pyrrole (aromatic ring with NH)
    let molecule = parse("c1cc[nH]c1").expect("Failed to parse pyrrole");

    assert_eq!(molecule.nodes().len(), 5);
    assert_eq!(molecule.bonds().len(), 5);

    // Check the nitrogen with its hydrogen
    let nitrogen = molecule
        .nodes()
        .iter()
        .find(|n| *n.atom().element() == AtomSymbol::Organic(OrganicAtom::N));
    assert!(nitrogen.is_some());
    assert_eq!(nitrogen.unwrap().hydrogens(), 1);
}

#[test]
fn parse_imidazole() {
    // c1cnc[nH]1 = imidazole (two nitrogens in the ring)
    let molecule = parse("c1cnc[nH]1").expect("Failed to parse imidazole");

    assert_eq!(molecule.nodes().len(), 5);

    // Count the nitrogens
    let nitrogens: Vec<_> = molecule
        .nodes()
        .iter()
        .filter(|n| *n.atom().element() == AtomSymbol::Organic(OrganicAtom::N))
        .collect();
    assert_eq!(nitrogens.len(), 2);
}

#[test]
fn parse_indole() {
    // c1ccc2[nH]ccc2c1 = indole (benzene fused with pyrrole)
    let molecule = parse("c1ccc2[nH]ccc2c1").expect("Failed to parse indole");

    assert_eq!(molecule.nodes().len(), 9);
    // 9 atoms with 10 bonds (two fused rings)
    assert_eq!(molecule.bonds().len(), 10);
}

#[test]
fn parse_mixed_aromatic_aliphatic() {
    // Cc1ccccc1 = toluene (methyl + benzene)
    let molecule = parse("Cc1ccccc1").expect("Failed to parse toluene");

    assert_eq!(molecule.nodes().len(), 7);

    // The first carbon (methyl) is not aromatic
    assert!(!molecule.nodes()[0].aromatic());

    // The other 6 carbons are aromatic
    for i in 1..7 {
        assert!(molecule.nodes()[i].aromatic());
    }
}

#[test]
fn parse_biphenyl() {
    // c1ccccc1-c2ccccc2 = biphenyl: two benzenes connected by a SIMPLE bond
    let molecule = parse("c1ccccc1-c2ccccc2").expect("Failed to parse biphenyl");

    // 12 aromatic carbons
    assert_eq!(molecule.nodes().len(), 12);

    // 13 bonds: 6 + 6 (rings) + 1 (between rings)
    assert_eq!(molecule.bonds().len(), 13);

    // All atoms are aromatic
    for node in molecule.nodes() {
        assert!(node.aromatic());
    }

    // Count bond types
    let aromatic_bonds = molecule
        .bonds()
        .iter()
        .filter(|b| b.kind() == BondType::Aromatic)
        .count();
    let simple_bonds = molecule
        .bonds()
        .iter()
        .filter(|b| b.kind() == BondType::Simple)
        .count();

    // 12 aromatic bonds (6 per ring) + 1 simple bond (between rings)
    assert_eq!(aromatic_bonds, 12);
    assert_eq!(
        simple_bonds, 1,
        "The bond between the two benzenes must be simple, not aromatic"
    );
}
