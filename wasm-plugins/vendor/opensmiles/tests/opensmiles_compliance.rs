//! OpenSMILES specification compliance tests.
//!
//! Tests for features required by the OpenSMILES spec that were previously
//! missing or incorrect.

use opensmiles::{parse, AtomSymbol, BondType, ParserError};

// ============================================================================
// Aromatic two-letter bracket symbols: [se], [as], [te]
// OpenSMILES spec: aromatic_symbols ::= 'b' | 'c' | 'n' | 'o' | 'p' | 's' | 'se' | 'as'
// ============================================================================

#[test]
fn parse_aromatic_selenium_bracket() {
    // [se] is aromatic selenium per OpenSMILES
    // c1cc[se]cc1 is a 6-membered ring with Se → 5×1 + 2 = 7 pi electrons (not 4n+2)
    // Without huckel-validation: parses successfully (syntax is valid)
    // With huckel-validation: rejected as chemically invalid
    let result = parse("c1cc[se]cc1");
    #[cfg(not(feature = "huckel-validation"))]
    {
        let mol = result.expect("Failed to parse c1cc[se]cc1");
        assert_eq!(mol.nodes().len(), 6);
        assert_eq!(*mol.nodes()[3].atom().element(), AtomSymbol::Se);
        assert!(mol.nodes()[3].aromatic());
    }
    #[cfg(feature = "huckel-validation")]
    {
        assert!(
            result.is_err(),
            "c1cc[se]cc1 should fail Hückel validation (7 pi electrons)"
        );
    }
}

#[test]
fn parse_aromatic_arsenic_bracket() {
    // [as] is aromatic arsenic per OpenSMILES
    let mol = parse("c1cc[as]cc1").expect("Failed to parse c1cc[as]cc1");
    assert_eq!(mol.nodes().len(), 6);
    assert_eq!(*mol.nodes()[3].atom().element(), AtomSymbol::As);
    assert!(mol.nodes()[3].aromatic());
}

#[test]
fn parse_aromatic_tellurium_bracket() {
    // [te] is aromatic tellurium (Te has can_be_aromatic = true)
    let mol = parse("[te]").expect("Failed to parse [te]");
    assert_eq!(mol.nodes().len(), 1);
    assert_eq!(*mol.nodes()[0].atom().element(), AtomSymbol::Te);
    assert!(mol.nodes()[0].aromatic());
}

#[test]
fn parse_selenophene() {
    // Selenophene: aromatic ring with selenium
    let mol = parse("c1cc[se]c1").expect("Failed to parse selenophene");
    assert_eq!(mol.nodes().len(), 5);
}

// ============================================================================
// [HH1] rejection: hydrogen atom cannot have hydrogen count
// OpenSMILES spec: "A hydrogen atom cannot have a hydrogen count"
// ============================================================================

#[test]
fn reject_hydrogen_with_hydrogen_count() {
    let result = parse("[HH1]");
    assert!(result.is_err());
    assert!(matches!(
        result.unwrap_err(),
        ParserError::HydrogenWithHydrogenCount
    ));
}

#[test]
fn reject_hydrogen_with_implicit_h1() {
    // [HH] means H with H1 (implicit 1), should be illegal
    let result = parse("[HH]");
    assert!(result.is_err());
    assert!(matches!(
        result.unwrap_err(),
        ParserError::HydrogenWithHydrogenCount
    ));
}

#[test]
fn accept_hydrogen_with_zero_h() {
    // [HH0] or just [H] should be fine (explicit H0)
    let mol = parse("[H]").expect("Failed to parse [H]");
    assert_eq!(mol.nodes().len(), 1);
    assert_eq!(*mol.nodes()[0].atom().element(), AtomSymbol::H);
    assert_eq!(mol.nodes()[0].hydrogens(), 0);
}

#[test]
fn accept_molecular_hydrogen() {
    // [H][H] — molecular hydrogen (two explicit H atoms bonded)
    let mol = parse("[H][H]").expect("Failed to parse [H][H]");
    assert_eq!(mol.nodes().len(), 2);
    assert_eq!(mol.bonds().len(), 1);
}

// ============================================================================
// Self-bond rejection: C11 is illegal
// OpenSMILES spec: "an atom cannot be bonded to itself"
// ============================================================================

#[test]
fn reject_self_bond() {
    let result = parse("C11");
    assert!(result.is_err());
    assert!(matches!(result.unwrap_err(), ParserError::SelfBond(1)));
}

#[test]
fn reject_self_bond_zero() {
    let result = parse("C00");
    assert!(result.is_err());
    assert!(matches!(result.unwrap_err(), ParserError::SelfBond(0)));
}

// ============================================================================
// Duplicate bond rejection: C12CCCCC12 is illegal
// OpenSMILES spec: "Two atoms cannot be joined by more than one bond"
// ============================================================================

#[test]
fn reject_duplicate_bond() {
    let result = parse("C12CCCCC12");
    assert!(result.is_err());
    assert!(matches!(
        result.unwrap_err(),
        ParserError::DuplicateBond(_, _)
    ));
}

#[test]
fn reject_duplicate_bond_two_char() {
    // Same with two-digit ring numbers
    let result = parse("C%12%13CCCCC%12%13");
    assert!(result.is_err());
    assert!(matches!(
        result.unwrap_err(),
        ParserError::DuplicateBond(_, _)
    ));
}

// ============================================================================
// Whitespace terminator
// OpenSMILES spec: terminator ::= SPACE | TAB | LINEFEED | CARRIAGE_RETURN | END_OF_STRING
// ============================================================================

#[test]
fn parse_with_trailing_space() {
    let mol = parse("CCO ").expect("Failed to parse 'CCO '");
    assert_eq!(mol.nodes().len(), 3);
}

#[test]
fn parse_with_trailing_tab() {
    let mol = parse("CCO\t").expect("Failed to parse 'CCO\\t'");
    assert_eq!(mol.nodes().len(), 3);
}

#[test]
fn parse_with_trailing_newline() {
    let mol = parse("CCO\n").expect("Failed to parse 'CCO\\n'");
    assert_eq!(mol.nodes().len(), 3);
}

#[test]
fn parse_with_trailing_cr() {
    let mol = parse("CCO\r").expect("Failed to parse 'CCO\\r'");
    assert_eq!(mol.nodes().len(), 3);
}

#[test]
fn parse_with_trailing_text_after_space() {
    // Content after whitespace should be ignored (e.g., SMILES with name)
    let mol = parse("CCO ethanol").expect("Failed to parse 'CCO ethanol'");
    assert_eq!(mol.nodes().len(), 3);
}

#[test]
fn parse_with_tab_then_name() {
    let mol = parse("c1ccccc1\tbenzene").expect("Failed to parse 'c1ccccc1\\tbenzene'");
    assert_eq!(mol.nodes().len(), 6);
}

// ============================================================================
// Disconnected structures: dot creates no bond
// OpenSMILES spec: dot is NOT a bond type in the grammar
// ============================================================================

#[test]
fn parse_dot_creates_no_bond() {
    let mol = parse("C.C").expect("Failed to parse C.C");
    assert_eq!(mol.nodes().len(), 2);
    assert_eq!(mol.bonds().len(), 0);
}

#[test]
fn parse_dot_in_branch_creates_no_bond() {
    // C(.C) = carbon with a disconnected carbon in branch
    let mol = parse("C(.C)").expect("Failed to parse C(.C)");
    assert_eq!(mol.nodes().len(), 2);
    assert_eq!(mol.bonds().len(), 0);
}

#[test]
fn parse_dot_ionic_compound() {
    let mol = parse("[Na+].[Cl-]").expect("Failed to parse [Na+].[Cl-]");
    assert_eq!(mol.nodes().len(), 2);
    assert_eq!(mol.bonds().len(), 0);
    assert_eq!(mol.nodes()[0].atom().charge(), 1);
    assert_eq!(mol.nodes()[1].atom().charge(), -1);
}

#[test]
fn parse_dot_multiple_fragments() {
    // Three separate methane molecules
    let mol = parse("C.C.C").expect("Failed to parse C.C.C");
    assert_eq!(mol.nodes().len(), 3);
    assert_eq!(mol.bonds().len(), 0);
}

#[test]
fn parse_dot_ethane_pair() {
    // CC.CC = two ethane molecules
    let mol = parse("CC.CC").expect("Failed to parse CC.CC");
    assert_eq!(mol.nodes().len(), 4);
    assert_eq!(mol.bonds().len(), 2);
    // First ethane: bond 0-1
    assert_eq!(mol.bonds()[0].source(), 0);
    assert_eq!(mol.bonds()[0].target(), 1);
    // Second ethane: bond 2-3
    assert_eq!(mol.bonds()[1].source(), 2);
    assert_eq!(mol.bonds()[1].target(), 3);
    // No bond crosses the dot
    for bond in mol.bonds() {
        assert_ne!(
            (bond.source(), bond.target()),
            (1, 2),
            "Bond should not cross the dot separator"
        );
    }
}

// ============================================================================
// Ring number matching: %01 == 1, ring number 0
// OpenSMILES spec: digits are interpreted as numbers, not symbols
// ============================================================================

#[test]
fn ring_number_zero() {
    let mol = parse("C0CCCCC0").expect("Failed to parse C0CCCCC0");
    assert_eq!(mol.nodes().len(), 6);
    assert_eq!(mol.bonds().len(), 6); // 5 chain + 1 ring closure
}

#[test]
fn ring_percent_01_matches_1() {
    let mol = parse("C1CCCCC%01").expect("Failed to parse C1CCCCC%01");
    assert_eq!(mol.nodes().len(), 6);
    assert_eq!(mol.bonds().len(), 6);
}

// ============================================================================
// Deprecated charge notation
// OpenSMILES spec: ++ and -- should be accepted for backwards compatibility
// ============================================================================

#[test]
fn parse_deprecated_plus_plus() {
    let mol = parse("[Cu++]").expect("Failed to parse [Cu++]");
    assert_eq!(mol.nodes()[0].atom().charge(), 2);
}

#[test]
fn parse_deprecated_minus_minus() {
    let mol = parse("[O--]").expect("Failed to parse [O--]");
    assert_eq!(mol.nodes()[0].atom().charge(), -2);
}

// ============================================================================
// Bond at ring open/close
// OpenSMILES spec: bond symbol can be on either or both sides
// ============================================================================

#[test]
fn bond_at_ring_open() {
    let mol = parse("C=1CCCCC1").expect("Failed to parse C=1CCCCC1");
    assert_eq!(mol.nodes().len(), 6);
    // Should have a double bond for the ring closure
    let ring_bond = mol
        .bonds()
        .iter()
        .find(|b| (b.source() == 0 && b.target() == 5) || (b.source() == 5 && b.target() == 0));
    assert!(ring_bond.is_some());
    assert_eq!(ring_bond.unwrap().kind(), BondType::Double);
}

#[test]
fn bond_at_ring_close() {
    let mol = parse("C1CCCCC=1").expect("Failed to parse C1CCCCC=1");
    assert_eq!(mol.nodes().len(), 6);
    let ring_bond = mol
        .bonds()
        .iter()
        .find(|b| (b.source() == 0 && b.target() == 5) || (b.source() == 5 && b.target() == 0));
    assert!(ring_bond.is_some());
    assert_eq!(ring_bond.unwrap().kind(), BondType::Double);
}

#[test]
fn mismatched_ring_bond_rejected() {
    let result = parse("C-1CCCCC=1");
    assert!(result.is_err());
    assert!(matches!(
        result.unwrap_err(),
        ParserError::MismatchedRingBond(1)
    ));
}

// ============================================================================
// Wildcard atom
// ============================================================================

#[test]
fn wildcard_bare() {
    let mol = parse("*CC*").expect("Failed to parse *CC*");
    assert_eq!(mol.nodes().len(), 4);
    assert_eq!(*mol.nodes()[0].atom().element(), AtomSymbol::Wildcard);
    assert_eq!(*mol.nodes()[3].atom().element(), AtomSymbol::Wildcard);
}

#[test]
fn wildcard_bracket() {
    let mol = parse("[*]").expect("Failed to parse [*]");
    assert_eq!(mol.nodes().len(), 1);
    assert_eq!(*mol.nodes()[0].atom().element(), AtomSymbol::Wildcard);
}
