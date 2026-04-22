//! Branch tests
//!
//! These tests verify the parsing of branches in molecules:
//! - Simple branches `CC(C)C`
//! - Multiple branches `CC(C)(C)C`
//! - Nested branches `CC(C(C)C)C`
//! - Branches with different bond types

use opensmiles::{parse, BondType, ParserError};

#[test]
fn parse_simple_branch() {
    // CC(C)C = isobutane (2-methylpropane)
    let molecule = parse("CC(C)C").expect("Failed to parse isobutane");

    assert_eq!(molecule.nodes().len(), 4);
    assert_eq!(molecule.bonds().len(), 3);

    // Check connectivity: atom 1 is connected to 0, 2, and 3
    let bonds: Vec<_> = molecule.bonds().iter().collect();

    // Bond 0: C(0) - C(1)
    assert_eq!(bonds[0].source(), 0);
    assert_eq!(bonds[0].target(), 1);

    // Bond 1: C(1) - C(2) (branch)
    assert_eq!(bonds[1].source(), 1);
    assert_eq!(bonds[1].target(), 2);

    // Bond 2: C(1) - C(3)
    assert_eq!(bonds[2].source(), 1);
    assert_eq!(bonds[2].target(), 3);
}

#[test]
fn parse_multiple_branches() {
    // CC(C)(C)C = neopentane (2,2-dimethylpropane)
    let molecule = parse("CC(C)(C)C").expect("Failed to parse neopentane");

    assert_eq!(molecule.nodes().len(), 5);
    assert_eq!(molecule.bonds().len(), 4);

    // The central atom (index 1) must have 4 bonds
    let central_bonds: Vec<_> = molecule
        .bonds()
        .iter()
        .filter(|b| b.source() == 1 || b.target() == 1)
        .collect();
    assert_eq!(central_bonds.len(), 4);
}

#[test]
fn parse_nested_branches() {
    // CC(C(C)C)C = 2,3-dimethylbutane
    let molecule = parse("CC(C(C)C)C").expect("Failed to parse 2,3-dimethylbutane");

    assert_eq!(molecule.nodes().len(), 6);
    assert_eq!(molecule.bonds().len(), 5);
}

#[test]
fn parse_branch_with_double_bond() {
    // CC(=O)O = acetic acid
    let molecule = parse("CC(=O)O").expect("Failed to parse acetic acid");

    assert_eq!(molecule.nodes().len(), 4);
    assert_eq!(molecule.bonds().len(), 3);

    // Find the double bond C=O
    let double_bond = molecule
        .bonds()
        .iter()
        .find(|b| b.kind() == BondType::Double);
    assert!(double_bond.is_some(), "Should have a double bond");
}

#[test]
fn parse_branch_with_triple_bond() {
    // CC(C#N)C = isobutyronitrile
    let molecule = parse("CC(C#N)C").expect("Failed to parse molecule with triple bond in branch");

    // Check that there is a triple bond
    let triple_bond = molecule
        .bonds()
        .iter()
        .find(|b| b.kind() == BondType::Triple);
    assert!(triple_bond.is_some(), "Should have a triple bond");
}

#[test]
fn parse_long_branch() {
    // C(CCCC)C = hexane with long branch
    let molecule = parse("C(CCCC)C").expect("Failed to parse molecule with long branch");

    assert_eq!(molecule.nodes().len(), 6);
    assert_eq!(molecule.bonds().len(), 5);
}

#[test]
fn parse_branch_at_start() {
    // (C)CC = propane (branch at start, equivalent to CCC)
    let molecule = parse("(C)CC").expect("Failed to parse branch at start");

    assert_eq!(molecule.nodes().len(), 3);
    assert_eq!(molecule.bonds().len(), 2);
}

#[test]
fn parse_empty_branch() {
    // C()C = should be equivalent to CC or error depending on implementation
    // According to OpenSMILES, empty branches are not valid
    let result = parse("C()C");
    assert!(result.is_err(), "Empty branches should not be allowed");
}

#[test]
fn error_position_in_branch() {
    // C([C+X])C - X is an invalid character in the bracket atom
    // Position: C=1, (=2, [=3, C=4, +=5, X=6
    // The error should indicate position 6 (absolute position, not relative to branch)
    match parse("C([C+X])C") {
        Err(ParserError::UnexpectedCharacter(c, pos)) => {
            assert_eq!(c, 'X', "Expected unexpected character 'X'");
            assert_eq!(
                pos, 6,
                "Position should be absolute (6), not relative to branch"
            );
        }
        Ok(_) => panic!("Expected error, but parsing succeeded"),
        Err(other) => panic!("Expected UnexpectedCharacter, got {:?}", other),
    }
}

// ============================================================================
// Regression tests for deeply nested and complex branched structures
// These tests verify the fix for the add_branch index calculation bug
// ============================================================================

#[test]
fn parse_dendrimer_generation_1() {
    // C(C)(C) - simple 3-atom star
    let molecule = parse("C(C)(C)").expect("Failed to parse dendrimer gen 1");

    assert_eq!(molecule.nodes().len(), 3, "Should have 3 atoms");
    assert_eq!(molecule.bonds().len(), 2, "Should have 2 bonds");

    // Central atom (0) connected to both branches (1 and 2)
    let bonds: Vec<_> = molecule.bonds().iter().collect();
    assert_eq!(bonds[0].source(), 0);
    assert_eq!(bonds[0].target(), 1);
    assert_eq!(bonds[1].source(), 0);
    assert_eq!(bonds[1].target(), 2);
}

#[test]
fn parse_dendrimer_generation_2() {
    // C(C(C)(C))(C(C)(C)) - 7 atoms, tree structure
    let molecule = parse("C(C(C)(C))(C(C)(C))").expect("Failed to parse dendrimer gen 2");

    assert_eq!(molecule.nodes().len(), 7, "Should have 7 atoms (2^3 - 1)");
    assert_eq!(molecule.bonds().len(), 6, "Should have 6 bonds");
}

#[test]
fn parse_dendrimer_generation_3() {
    // C(C(C(C)(C))(C(C)(C)))(C(C(C)(C))(C(C)(C))) - 15 atoms
    let smiles = "C(C(C(C)(C))(C(C)(C)))(C(C(C)(C))(C(C)(C)))";
    let molecule = parse(smiles).expect("Failed to parse dendrimer gen 3");

    assert_eq!(molecule.nodes().len(), 15, "Should have 15 atoms (2^4 - 1)");
    assert_eq!(molecule.bonds().len(), 14, "Should have 14 bonds");
}

#[test]
fn parse_dendrimer_generation_4() {
    // 31 atoms dendrimer
    let smiles = "C(C(C(C(C)(C))(C(C)(C)))(C(C(C)(C))(C(C)(C))))(C(C(C(C)(C))(C(C)(C)))(C(C(C)(C))(C(C)(C))))";
    let molecule = parse(smiles).expect("Failed to parse dendrimer gen 4");

    assert_eq!(molecule.nodes().len(), 31, "Should have 31 atoms (2^5 - 1)");
    assert_eq!(molecule.bonds().len(), 30, "Should have 30 bonds");
}

#[test]
fn parse_deeply_nested_linear_branch() {
    // C(C(C(C(C(C)))))  - deeply nested but linear
    let molecule = parse("C(C(C(C(C(C)))))").expect("Failed to parse deeply nested branch");

    assert_eq!(molecule.nodes().len(), 6, "Should have 6 atoms");
    assert_eq!(molecule.bonds().len(), 5, "Should have 5 bonds");

    // Verify linear connectivity
    for (i, bond) in molecule.bonds().iter().enumerate() {
        assert_eq!(bond.source() as usize, i, "Source should be sequential");
        assert_eq!(bond.target() as usize, i + 1, "Target should be source + 1");
    }
}

#[test]
fn parse_star_molecule_10_branches() {
    // C(C)(C)(C)(C)(C)(C)(C)(C)(C)(C) - 11 atoms, 10 branches from central
    let mut smiles = String::from("C");
    for _ in 0..10 {
        smiles.push_str("(C)");
    }

    let molecule = parse(&smiles).expect("Failed to parse star with 10 branches");

    assert_eq!(molecule.nodes().len(), 11, "Should have 11 atoms");
    assert_eq!(molecule.bonds().len(), 10, "Should have 10 bonds");

    // All bonds should originate from atom 0
    for bond in molecule.bonds() {
        assert_eq!(bond.source(), 0, "All bonds should start from central atom");
    }
}

#[test]
fn parse_comb_polymer_structure() {
    // CC(C)C(C)C(C)C(C)C - main chain with pendant groups
    let molecule = parse("CC(C)C(C)C(C)C(C)C").expect("Failed to parse comb structure");

    assert_eq!(molecule.nodes().len(), 10, "Should have 10 atoms");
    assert_eq!(molecule.bonds().len(), 9, "Should have 9 bonds");
}

#[test]
fn parse_branch_with_internal_bonds() {
    // C(CCC)C - branch with internal chain
    let molecule = parse("C(CCC)C").expect("Failed to parse branch with internal bonds");

    assert_eq!(molecule.nodes().len(), 5, "Should have 5 atoms");
    assert_eq!(molecule.bonds().len(), 4, "Should have 4 bonds");

    // Verify internal branch connectivity
    let bonds: Vec<_> = molecule.bonds().iter().collect();
    // Bond 0: main(0) -> branch_start(1)
    assert_eq!(bonds[0].source(), 0);
    assert_eq!(bonds[0].target(), 1);
    // Bond 1: branch internal 1->2
    assert_eq!(bonds[1].source(), 1);
    assert_eq!(bonds[1].target(), 2);
    // Bond 2: branch internal 2->3
    assert_eq!(bonds[2].source(), 2);
    assert_eq!(bonds[2].target(), 3);
    // Bond 3: main(0) -> after_branch(4)
    assert_eq!(bonds[3].source(), 0);
    assert_eq!(bonds[3].target(), 4);
}
