//! Bracket atom tests
//!
//! These tests verify the parsing of bracket atoms with:
//! - Explicit hydrogens (`[CH4]`, `[CH3]`, etc.)
//! - Charges (`[NH4+]`, `[O-]`, `[Fe+2]`, etc.)
//! - Isotopes (`[13C]`, `[2H]`, etc.)
//! - Classes/atom mapping (`[C:1]`, `[N:2]`, etc.)
//! - All possible combinations of these attributes

use opensmiles::{parse, AtomSymbol, MoleculeError, NodeError, OrganicAtom, ParserError};

// ============================================================================
// Explicit hydrogens
// ============================================================================

#[test]
fn parse_bracket_atom_no_hydrogen() {
    // [C] = carbon without explicit hydrogen (radical)
    let molecule = parse("[C]").expect("Failed to parse [C]");

    assert_eq!(molecule.nodes().len(), 1);

    let node = &molecule.nodes()[0];
    assert_eq!(*node.atom().element(), AtomSymbol::Organic(OrganicAtom::C));
    assert_eq!(node.hydrogens(), 0);
    assert_eq!(node.atom().charge(), 0);
    assert_eq!(node.atom().isotope(), None);
    assert_eq!(node.class(), None);
}

#[test]
fn parse_bracket_atom_with_hydrogen_count() {
    // [CH4] = methane with explicit hydrogens
    let molecule = parse("[CH4]").expect("Failed to parse [CH4]");

    assert_eq!(molecule.nodes().len(), 1);

    let node = &molecule.nodes()[0];
    assert_eq!(*node.atom().element(), AtomSymbol::Organic(OrganicAtom::C));
    assert_eq!(node.hydrogens(), 4);
}

#[test]
fn parse_bracket_atom_hydrogen_implicit_one() {
    // [CH] = carbon with 1 hydrogen (H without number = 1)
    let molecule = parse("[CH]").expect("Failed to parse [CH]");

    assert_eq!(molecule.nodes().len(), 1);
    assert_eq!(molecule.nodes()[0].hydrogens(), 1);
}

#[test]
fn parse_bracket_hydrogen_variations() {
    // Test different hydrogen values
    let test_cases = [
        ("[CH]", 1),
        ("[CH2]", 2),
        ("[CH3]", 3),
        ("[NH2]", 2),
        ("[OH]", 1),
        ("[SH]", 1),
    ];

    for (smiles, expected_h) in test_cases {
        let molecule = parse(smiles).unwrap_or_else(|_| panic!("Failed to parse {}", smiles));
        assert_eq!(
            molecule.nodes()[0].hydrogens(),
            expected_h,
            "Wrong hydrogen count for {}",
            smiles
        );
    }
}

// ============================================================================
// Charges
// ============================================================================

#[test]
fn parse_bracket_atom_positive_charge() {
    // [NH4+] = ammonium ion
    let molecule = parse("[NH4+]").expect("Failed to parse [NH4+]");

    assert_eq!(molecule.nodes().len(), 1);

    let node = &molecule.nodes()[0];
    assert_eq!(*node.atom().element(), AtomSymbol::Organic(OrganicAtom::N));
    assert_eq!(node.atom().charge(), 1);
    assert_eq!(node.hydrogens(), 4);
}

#[test]
fn parse_bracket_atom_negative_charge() {
    // [O-] = oxide ion
    let molecule = parse("[O-]").expect("Failed to parse [O-]");

    assert_eq!(molecule.nodes().len(), 1);

    let node = &molecule.nodes()[0];
    assert_eq!(*node.atom().element(), AtomSymbol::Organic(OrganicAtom::O));
    assert_eq!(node.atom().charge(), -1);
}

#[test]
fn parse_bracket_atom_multiple_charge_plus() {
    // [Fe++] or [Fe+2] = iron(II) ion
    let molecule1 = parse("[Fe++]").expect("Failed to parse [Fe++]");
    let molecule2 = parse("[Fe+2]").expect("Failed to parse [Fe+2]");

    assert_eq!(molecule1.nodes()[0].atom().charge(), 2);
    assert_eq!(molecule2.nodes()[0].atom().charge(), 2);
}

#[test]
fn parse_bracket_atom_multiple_charge_minus() {
    // [O--] or [O-2] = doubly charged oxide ion
    let molecule1 = parse("[O--]").expect("Failed to parse [O--]");
    let molecule2 = parse("[O-2]").expect("Failed to parse [O-2]");

    assert_eq!(molecule1.nodes()[0].atom().charge(), -2);
    assert_eq!(molecule2.nodes()[0].atom().charge(), -2);
}

#[test]
fn parse_bracket_atom_high_charge() {
    // High charges (up to +/-15 according to specification)
    let molecule_pos = parse("[Fe+3]").expect("Failed to parse [Fe+3]");
    let molecule_neg = parse("[P-3]").expect("Failed to parse [P-3]");

    assert_eq!(molecule_pos.nodes()[0].atom().charge(), 3);
    assert_eq!(molecule_neg.nodes()[0].atom().charge(), -3);
}

// ============================================================================
// Isotopes
// ============================================================================

#[test]
fn parse_bracket_atom_isotope() {
    // [13C] = carbon-13
    let molecule = parse("[13C]").expect("Failed to parse [13C]");

    assert_eq!(molecule.nodes().len(), 1);

    let node = &molecule.nodes()[0];
    assert_eq!(*node.atom().element(), AtomSymbol::Organic(OrganicAtom::C));
    assert_eq!(node.atom().isotope(), Some(13));
}

#[test]
fn parse_bracket_atom_deuterium() {
    // [2H] = deuterium
    let molecule = parse("[2H]").expect("Failed to parse [2H]");

    assert_eq!(molecule.nodes().len(), 1);
    assert_eq!(*molecule.nodes()[0].atom().element(), AtomSymbol::H);
    assert_eq!(molecule.nodes()[0].atom().isotope(), Some(2));
}

#[test]
fn parse_bracket_atom_tritium() {
    // [3H] = tritium
    let molecule = parse("[3H]").expect("Failed to parse [3H]");

    assert_eq!(molecule.nodes()[0].atom().isotope(), Some(3));
}

#[test]
fn parse_bracket_atom_oxygen_18() {
    // [18O] = oxygen-18
    let molecule = parse("[18O]").expect("Failed to parse [18O]");

    assert_eq!(
        *molecule.nodes()[0].atom().element(),
        AtomSymbol::Organic(OrganicAtom::O)
    );
    assert_eq!(molecule.nodes()[0].atom().isotope(), Some(18));
}

// ============================================================================
// Classes (atom mapping)
// ============================================================================

#[test]
fn parse_bracket_atom_class() {
    // [C:1] = carbon with class 1
    let molecule = parse("[C:1]").expect("Failed to parse [C:1]");

    assert_eq!(molecule.nodes().len(), 1);

    let node = &molecule.nodes()[0];
    assert_eq!(*node.atom().element(), AtomSymbol::Organic(OrganicAtom::C));
    assert_eq!(node.class(), Some(1));
}

#[test]
fn parse_bracket_atom_class_variations() {
    // Test different class values
    let test_cases = [("[C:0]", 0), ("[C:1]", 1), ("[C:42]", 42), ("[C:999]", 999)];

    for (smiles, expected_class) in test_cases {
        let molecule = parse(smiles).unwrap_or_else(|_| panic!("Failed to parse {}", smiles));
        assert_eq!(
            molecule.nodes()[0].class(),
            Some(expected_class),
            "Wrong class for {}",
            smiles
        );
    }
}

// ============================================================================
// Combinations of 2 attributes
// ============================================================================

#[test]
fn parse_bracket_atom_hydrogen_and_charge() {
    // [NH4+] = ammonium: hydrogens + charge
    let molecule = parse("[NH4+]").expect("Failed to parse [NH4+]");

    let node = &molecule.nodes()[0];
    assert_eq!(node.hydrogens(), 4);
    assert_eq!(node.atom().charge(), 1);
}

#[test]
fn parse_bracket_atom_isotope_and_hydrogen() {
    // [13CH4] = methane labeled with carbon-13
    let molecule = parse("[13CH4]").expect("Failed to parse [13CH4]");

    let node = &molecule.nodes()[0];
    assert_eq!(node.atom().isotope(), Some(13));
    assert_eq!(node.hydrogens(), 4);
}

#[test]
fn parse_bracket_atom_isotope_and_charge() {
    // [13C-] = anionic carbon-13
    let molecule = parse("[13C-]").expect("Failed to parse [13C-]");

    let node = &molecule.nodes()[0];
    assert_eq!(node.atom().isotope(), Some(13));
    assert_eq!(node.atom().charge(), -1);
}

#[test]
fn parse_bracket_atom_hydrogen_and_class() {
    // [CH3:1] = methyl with class
    let molecule = parse("[CH3:1]").expect("Failed to parse [CH3:1]");

    let node = &molecule.nodes()[0];
    assert_eq!(node.hydrogens(), 3);
    assert_eq!(node.class(), Some(1));
}

#[test]
fn parse_bracket_atom_out_of_range_hydrogen() {
    let result = parse("[CH256]");
    dbg!(&result);
    assert!(matches!(&result, Err(ParserError::HydrogenOutOfRange(s)) if s == "256"));
}

#[test]
fn parse_bracket_atom_too_much_hydrogen() {
    let result = parse("[CH10]");
    dbg!(&result);
    assert!(matches!(
        &result,
        Err(ParserError::MoleculeError(MoleculeError::NodeError(
            NodeError::InvalidHydrogen(10)
        )))
    ));
}

#[test]
fn parse_bracket_atom_charge_and_class() {
    // [O-:2] = oxide with class
    let molecule = parse("[O-:2]").expect("Failed to parse [O-:2]");

    let node = &molecule.nodes()[0];
    assert_eq!(node.atom().charge(), -1);
    assert_eq!(node.class(), Some(2));
}

#[test]
fn parse_bracket_atom_isotope_and_class() {
    // [13C:1] = carbon-13 with class
    let molecule = parse("[13C:1]").expect("Failed to parse [13C:1]");

    let node = &molecule.nodes()[0];
    assert_eq!(node.atom().isotope(), Some(13));
    assert_eq!(node.class(), Some(1));
}

// ============================================================================
// Combinations of 3 attributes
// ============================================================================

#[test]
fn parse_bracket_atom_isotope_hydrogen_charge() {
    // [13CH3+] = labeled methyl cation
    let molecule = parse("[13CH3+]").expect("Failed to parse [13CH3+]");

    let node = &molecule.nodes()[0];
    assert_eq!(node.atom().isotope(), Some(13));
    assert_eq!(node.hydrogens(), 3);
    assert_eq!(node.atom().charge(), 1);
}

#[test]
fn parse_bracket_atom_isotope_hydrogen_class() {
    // [13CH4:1] = labeled methane with class
    let molecule = parse("[13CH4:1]").expect("Failed to parse [13CH4:1]");

    let node = &molecule.nodes()[0];
    assert_eq!(node.atom().isotope(), Some(13));
    assert_eq!(node.hydrogens(), 4);
    assert_eq!(node.class(), Some(1));
}

#[test]
fn parse_bracket_atom_isotope_charge_class() {
    // [13C+:1] = carbon-13 cation with class
    let molecule = parse("[13C+:1]").expect("Failed to parse [13C+:1]");

    let node = &molecule.nodes()[0];
    assert_eq!(node.atom().isotope(), Some(13));
    assert_eq!(node.atom().charge(), 1);
    assert_eq!(node.class(), Some(1));
}

#[test]
fn parse_bracket_atom_hydrogen_charge_class() {
    // [NH4+:1] = ammonium with class
    let molecule = parse("[NH4+:1]").expect("Failed to parse [NH4+:1]");

    let node = &molecule.nodes()[0];
    assert_eq!(node.hydrogens(), 4);
    assert_eq!(node.atom().charge(), 1);
    assert_eq!(node.class(), Some(1));
}

// ============================================================================
// Combination of all attributes (4)
// ============================================================================

#[test]
fn parse_bracket_atom_all_attributes() {
    // [13CH4+:1] = all properties combined
    let molecule = parse("[13CH4+:1]").expect("Failed to parse [13CH4+:1]");

    let node = &molecule.nodes()[0];
    assert_eq!(*node.atom().element(), AtomSymbol::Organic(OrganicAtom::C));
    assert_eq!(node.atom().isotope(), Some(13));
    assert_eq!(node.hydrogens(), 4);
    assert_eq!(node.atom().charge(), 1);
    assert_eq!(node.class(), Some(1));
}

#[test]
fn parse_bracket_atom_all_attributes_negative() {
    // [18OH-:5] = labeled hydroxide with class
    let molecule = parse("[18OH-:5]").expect("Failed to parse [18OH-:5]");

    let node = &molecule.nodes()[0];
    assert_eq!(*node.atom().element(), AtomSymbol::Organic(OrganicAtom::O));
    assert_eq!(node.atom().isotope(), Some(18));
    assert_eq!(node.hydrogens(), 1);
    assert_eq!(node.atom().charge(), -1);
    assert_eq!(node.class(), Some(5));
}

// ============================================================================
// Non-organic elements
// ============================================================================

#[test]
fn parse_bracket_metal() {
    // [Fe] = iron
    let molecule = parse("[Fe]").expect("Failed to parse [Fe]");

    assert_eq!(molecule.nodes().len(), 1);
    assert_eq!(*molecule.nodes()[0].atom().element(), AtomSymbol::Fe);
}

#[test]
fn parse_bracket_noble_gas() {
    // [He] = helium
    let molecule = parse("[He]").expect("Failed to parse [He]");

    assert_eq!(*molecule.nodes()[0].atom().element(), AtomSymbol::He);
}

#[test]
fn parse_bracket_lanthanide() {
    // [La] = lanthanum
    let molecule = parse("[La]").expect("Failed to parse [La]");

    assert_eq!(*molecule.nodes()[0].atom().element(), AtomSymbol::La);
}

// ============================================================================
// Errors
// ============================================================================

#[test]
fn parse_bracket_unexpected_char() {
    // [Fe] = iron
    match parse("[C+X]") {
        Err(ParserError::UnexpectedCharacter(c, pos)) => {
            assert_eq!(c, 'X');
            assert_eq!(pos, 4);
        }
        other => panic!("Expected UnexpectedCharacter, got {:?}", other),
    }
}
