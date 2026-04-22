use opensmiles::parse;

#[test]
fn display_simple_methane() {
    let mol = parse("C").unwrap();
    assert_eq!(format!("{}", mol), "C");
}

#[test]
fn display_bracket_methane() {
    let mol = parse("[CH4]").unwrap();
    assert_eq!(format!("{}", mol), "C");
}

#[test]
fn display_ethane() {
    let mol = parse("CC").unwrap();
    assert_eq!(format!("{}", mol), "CC");
}

#[test]
fn display_ethanol() {
    let mol = parse("CCO").unwrap();
    assert_eq!(format!("{}", mol), "OCC");
}

#[test]
fn display_methanone() {
    let mol = parse("C=O").unwrap();
    assert_eq!(format!("{}", mol), "O=C");
}

#[test]
fn display_methylamine() {
    let mol = parse("C#N").unwrap();
    assert_eq!(format!("{}", mol), "N#C");
}

#[test]
fn display_branch() {
    let mol = parse("CC(C)C").unwrap();
    assert_eq!(format!("{}", mol), "CC(C)C");
}

#[test]
fn display_branch_with_bond() {
    let mol = parse("CC(=O)O").unwrap();
    assert_eq!(format!("{}", mol), "O=C(C)O");
}

#[test]
fn display_benzene() {
    let mol = parse("C1CCCCC1").unwrap();
    assert_eq!(format!("{}", mol), "C1CCCCC1");
}

#[test]
fn display_benzene_with_unstandard_cycle_number() {
    let mol = parse("C4CCCCC4").unwrap();
    assert_eq!(format!("{}", mol), "C1CCCCC1");
}

#[test]
fn display_aromatic_cycle() {
    let mol = parse("c1ccccc1").unwrap();
    assert_eq!(format!("{}", mol), "c1ccccc1");
}

#[test]
fn display_toluene() {
    let mol = parse("Cc1ccccc1").unwrap();
    assert_eq!(format!("{}", mol), "Cc1ccccc1");
}

#[test]
fn display_neopentane() {
    let mol = parse("CC(C)(C)C").unwrap();
    assert_eq!(format!("{}", mol), "CC(C)(C)C");
}

#[test]
fn display_acetate_ion() {
    let mol = parse("CC(=O)[O-]").unwrap();
    assert_eq!(format!("{}", mol), "O=C(C)[O-]");
}

#[test]
fn display_sodium_ion() {
    let mol = parse("[Na+]").unwrap();
    assert_eq!(format!("{}", mol), "[Na+]");
}

#[test]
fn display_cyclohexanone() {
    let mol = parse("O=C1CCCCC1").unwrap();
    assert_eq!(format!("{}", mol), "O=C1CCCCC1");
}

#[test]
fn display_naphthalene() {
    let mol = parse("c1ccc2ccccc2c1").unwrap();
    assert_eq!(format!("{}", mol), "c1ccc2ccccc2c1");
}

#[test]
fn display_aspirin() {
    let mol = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
    assert_eq!(format!("{}", mol), "O=C(C)Oc1ccccc1C(=O)O");
}

#[test]
fn display_trans_2_butene() {
    let mol = parse("C/C=C/C").unwrap();
    assert_eq!(format!("{}", mol), "C/C=C/C");
}

#[test]
fn display_cis_2_butene() {
    let mol = parse("C/C=C\\C").unwrap();
    assert_eq!(format!("{}", mol), "C/C=C\\C");
}

#[test]
fn display_chiral_carbon() {
    let mol = parse("[C@@H](F)(Cl)Br").unwrap();
    assert_eq!(format!("{}", mol), "F[C@@H](Cl)Br");
}

#[test]
fn display_isotope() {
    let mol = parse("[13C]").unwrap();
    assert_eq!(format!("{}", mol), "[13C]");
}

#[test]
fn display_two_digits_ring_number_must_begin_with_percent_symbol() {
    let mol = parse("c1ccccc1Cc2ccccc2Cc3ccccc3Cc4ccccc4Cc5ccccc5Cc6ccccc6Cc7ccccc7Cc8ccccc8Cc9ccccc9Cc%10ccccc%10").unwrap();
    assert_eq!(format!("{}", mol), "c1ccccc1Cc2ccccc2Cc3ccccc3Cc4ccccc4Cc5ccccc5Cc6ccccc6Cc7ccccc7Cc8ccccc8Cc9ccccc9Cc%10ccccc%10");
}

#[test]
fn display_standard_form_atoms() {
    // Write atoms in the "organic subset" as bare atomic symbols whenever possible.
    let mol1 = parse("[CH3][CH3]").unwrap();
    assert_eq!(format!("{}", mol1), "CC");

    // If the charge is +1 or -1, leave off the digit.
    let mol2 = parse("[CH3-1]").unwrap();
    assert_eq!(format!("{}", mol2), "[CH3-]");

    // If the hydrogen count is 1, leave off the digit.
    let mol3 = parse("C[13CH1](C)C").unwrap();
    assert_eq!(format!("{}", mol3), "C[13CH](C)C");

    // Always write the atom properties in the order: Chirality, hydrogen-count, charge.
    let mol4 = parse("[C-H3]").unwrap();
    assert_eq!(format!("{}", mol4), "[CH3-]");

    let mol5 = parse("F[CH@](Br)Cl").unwrap();
    assert_eq!(format!("{}", mol5), "F[C@H](Br)Cl");

    // Represent hydrogens as a property of the heavy atom rather than as explicit atoms, unless other rules (e.g. [2H]) require that the hydrogen be explicit.
    let mol6 = parse("[H][C-]([H])[H]").unwrap();
    assert_eq!(format!("{}", mol6), "[CH3-]");
}

#[test]
fn display_standard_form_bonds() {
    // Only write '-' (single bond) when it is between two aromatic atoms.
    let mol1 = parse("C-C").unwrap();
    assert_eq!(format!("{}", mol1), "CC");

    // Never write the ':' (aromatic bond) symbol. Bonds are single or aromatic by default (as appropriate).
    let mol2 = parse("c:1:c:c:c:c:c:1").unwrap();
    assert_eq!(format!("{}", mol2), "c1ccccc1");

    let mol3 = parse("c1ccccc1c2ccccc2").unwrap();
    assert_eq!(format!("{}", mol3), "c1ccccc1-c2ccccc2");
}

#[test]
fn display_standard_form_cycles() {
    // Don’t reuse ring-closure digits.
    let mol1 = parse("c1ccccc1C1CCCC1").unwrap();
    assert_eq!(format!("{}", mol1), "c1ccccc1C2CCCC2");

    // Begin ring numbering with 1, not zero (or any other number)
    let mol2 = parse("c0ccccc0C1CCCC1").unwrap();
    assert_eq!(format!("{}", mol2), "c1ccccc1C2CCCC2");

    // Avoid making a ring-closure on a double or triple bond. For the ring-closure digits, choose a single bond whenever possible.
    let mol3 = parse("CC=1CCCCC=1").unwrap();
    assert_eq!(format!("{}", mol3), "CC1=CCCCC1");

    // Avoid starting a ring system on an atom that is in two or more rings, such that two ring-closure bonds will be on the same atom.
    let mol4 = parse("C12(CCCCC1)CCCCC2").unwrap();
    assert_eq!(format!("{}", mol4), "C1C2(CCCC1)CCCCC2");

    // Use the simpler single-digit form for rnums less than 10.
    let mol5 = parse("C%01CCCCC%01").unwrap();
    assert_eq!(format!("{}", mol5), "C1CCCCC1");
}

#[test]
fn display_standard_form_branches() {
    // Start on a terminal atom if possible.
    let mol1 = parse("c1cc(CO)ccc1").unwrap();
    assert_eq!(format!("{}", mol1), "OCc1ccccc1");

    // Try to make "side chains" short; pick the longest chains as the "main branch" of the SMILES.
    let mol2 = parse("CC(CCCCCC)C").unwrap();
    assert_eq!(format!("{}", mol2), "CC(C)CCCCCC");

    // Start on a heteroatom if possible.
    let mol3 = parse("CCCO").unwrap();
    assert_eq!(format!("{}", mol3), "OCCC");

    // Only use dots for disconnected components.
    let mol4 = parse("C1.C1").unwrap();
    assert_eq!(format!("{}", mol4), "CC");
}

#[test]
fn display_standard_form_aromaticity() {
    // Write the aromatic form in preference to the Kekule form.
    let mol1 = parse("C1=CC=CC=C1").unwrap();
    assert_eq!(format!("{}", mol1), "c1ccccc1");
}

#[test]
fn display_standard_form_chirality() {
    // Remove chiral markings for atoms that are not chiral.
    let mol1 = parse("Br[C@H](Br)C").unwrap();
    assert_eq!(format!("{}", mol1), "BrC(Br)C");

    // Remove cis/trans markings for double bonds that are not cis or trans.
    let mol2 = parse("F/C(/F)=C/F").unwrap();
    assert_eq!(format!("{}", mol2), "FC(F)=CF");
}
