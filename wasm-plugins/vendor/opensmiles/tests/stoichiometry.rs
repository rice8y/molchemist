//! Stoichiometry tests (disconnected structures)
//!
//! These tests verify the parsing of molecules with the `.` separator:
//! - Simple disconnected molecules `C.C`
//! - Ionic compounds `[Na+].[Cl-]`
//! - Complex mixtures

use opensmiles::{parse, AtomSymbol, OrganicAtom};

#[test]
fn parse_disconnected_simple() {
    // C.C = two separate methane molecules
    let molecule = parse("C.C").expect("Failed to parse C.C");

    assert_eq!(molecule.nodes().len(), 2);
    assert_eq!(molecule.bonds().len(), 0); // No bond between the two
}

#[test]

fn parse_disconnected_ionic() {
    // [Na+].[Cl-] = sodium chloride (salt)
    let molecule = parse("[Na+].[Cl-]").expect("Failed to parse NaCl");

    assert_eq!(molecule.nodes().len(), 2);
    assert_eq!(molecule.bonds().len(), 0);

    assert_eq!(*molecule.nodes()[0].atom().element(), AtomSymbol::Na);
    assert_eq!(molecule.nodes()[0].atom().charge(), 1);

    assert_eq!(
        *molecule.nodes()[1].atom().element(),
        AtomSymbol::Organic(OrganicAtom::Cl)
    );
    assert_eq!(molecule.nodes()[1].atom().charge(), -1);
}

#[test]

fn parse_disconnected_multiple() {
    // C.C.C = three methane molecules
    let molecule = parse("C.C.C").expect("Failed to parse C.C.C");

    assert_eq!(molecule.nodes().len(), 3);
    assert_eq!(molecule.bonds().len(), 0);
}

#[test]

fn parse_disconnected_complex() {
    // CC.CC = two ethane molecules
    let molecule = parse("CC.CC").expect("Failed to parse CC.CC");

    assert_eq!(molecule.nodes().len(), 4);
    assert_eq!(molecule.bonds().len(), 2); // One bond per ethane molecule

    // Check that bonds are in the correct molecules
    let bond0 = &molecule.bonds()[0];
    let bond1 = &molecule.bonds()[1];

    // First bond: between atoms 0 and 1
    assert_eq!(bond0.source(), 0);
    assert_eq!(bond0.target(), 1);

    // Second bond: between atoms 2 and 3
    assert_eq!(bond1.source(), 2);
    assert_eq!(bond1.target(), 3);
}

#[test]

fn parse_disconnected_with_branches() {
    // CC(C)C.CC = isobutane + ethane
    let molecule = parse("CC(C)C.CC").expect("Failed to parse CC(C)C.CC");

    assert_eq!(molecule.nodes().len(), 6);
    assert_eq!(molecule.bonds().len(), 4); // 3 for isobutane + 1 for ethane
}

#[test]

fn parse_water_hydronium() {
    // [OH2].[H+] = water + proton (acid)
    let molecule = parse("[OH2].[H+]").expect("Failed to parse water + hydronium");

    assert_eq!(molecule.nodes().len(), 2);
    assert_eq!(molecule.bonds().len(), 0);

    // Water
    assert_eq!(
        *molecule.nodes()[0].atom().element(),
        AtomSymbol::Organic(OrganicAtom::O)
    );
    assert_eq!(molecule.nodes()[0].hydrogens(), 2);

    // Proton
    assert_eq!(*molecule.nodes()[1].atom().element(), AtomSymbol::H);
    assert_eq!(molecule.nodes()[1].atom().charge(), 1);
}

#[test]

fn parse_hydrate() {
    // CCO.[OH2] = ethanol + water (representation of a hydrate)
    let molecule = parse("CCO.[OH2]").expect("Failed to parse ethanol hydrate");

    assert_eq!(molecule.nodes().len(), 4); // 3 for ethanol + 1 for water
    assert_eq!(molecule.bonds().len(), 2); // Only ethanol bonds
}

#[test]

fn parse_metal_complex() {
    // [Cu+2].[O-].[O-] = copper(II) oxide
    let molecule = parse("[Cu+2].[O-].[O-]").expect("Failed to parse copper oxide");

    assert_eq!(molecule.nodes().len(), 3);
    assert_eq!(molecule.bonds().len(), 0);

    assert_eq!(*molecule.nodes()[0].atom().element(), AtomSymbol::Cu);
    assert_eq!(molecule.nodes()[0].atom().charge(), 2);

    for i in 1..=2 {
        assert_eq!(
            *molecule.nodes()[i].atom().element(),
            AtomSymbol::Organic(OrganicAtom::O)
        );
        assert_eq!(molecule.nodes()[i].atom().charge(), -1);
    }
}
