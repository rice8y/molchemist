use super::graph::Ring;
use super::molecule::Molecule;
use crate::{BondType, MoleculeError, NodeIndex};
use std::collections::{HashSet, VecDeque};

/// Result of aromaticity validation for a single ring.
#[derive(Debug, Clone, PartialEq)]
pub struct AromaticityCheck {
    pub ring: Ring,
    pub pi_electrons: Option<u8>,
    pub is_valid: bool,
}

/// Counts the number of sigma bonds for a given atom.
///
/// For aromatic pi-electron calculation, we count:
/// - Each bond in the molecule as 1 sigma bond (regardless of bond type)
/// - All hydrogens (now correctly calculated by the fixed implicit H logic)
fn count_sigma_bonds(molecule: &Molecule, node_idx: NodeIndex) -> u8 {
    let node = &molecule.nodes()[node_idx as usize];

    // Count bonds to other atoms (each bond = 1 sigma)
    let bond_count = molecule
        .bonds()
        .iter()
        .filter(|bond| bond.source() == node_idx || bond.target() == node_idx)
        .count() as u8;

    // Add hydrogens (now correctly calculated with aromatic subvalence rule)
    bond_count + node.hydrogens()
}

/// Determines the pi electron contribution of an atom in an aromatic ring.
///
/// Calculates the contribution based on valence electrons, sigma bonds, and charge.
///
/// Returns `None` if the contribution cannot be determined (e.g., Wildcard).
pub(crate) fn pi_electron_contribution(molecule: &Molecule, node_idx: NodeIndex) -> Option<u8> {
    let node = &molecule.nodes()[node_idx as usize];
    let element = node.atom().element();
    let charge = node.atom().charge();
    let valence_electrons = element.valence_electrons();

    // Wildcard atoms have valence_electrons = 0, cannot be determined
    if valence_electrons == 0 {
        return None;
    }

    // An aromatic carbon contributes one electron in its neutral state even
    // when a mixed aromatic/Kekulé spelling makes one adjacent double bond
    // explicit and therefore suppresses its implicit hydrogen.
    if element.atomic_number() == 6 {
        return Some(match charge.cmp(&0) {
            std::cmp::Ordering::Less => 2,
            std::cmp::Ordering::Equal => 1,
            std::cmp::Ordering::Greater => 0,
        });
    }

    let sigma_bonds = count_sigma_bonds(molecule, node_idx);

    // Calculate electrons remaining after forming sigma bonds, adjusted for charge
    // Negative charge adds electrons, positive charge removes them
    let electrons_after_sigma = (valence_electrons as i16) - (sigma_bonds as i16) - (charge as i16);

    // Can't have negative electrons
    if electrons_after_sigma < 0 {
        return Some(0);
    }

    let electrons_after_sigma = electrons_after_sigma as u8;

    // For aromatic atoms (sp2 hybridized), determine pi contribution:
    // - 0 electrons: empty p orbital → 0 pi electrons (e.g., B, C+)
    // - 1 electron: half-filled p orbital → 1 pi electron (e.g., C neutral)
    // - 2 electrons: filled p orbital or lone pair → 2 pi electrons (e.g., [nH], O, C-)
    // - 3 electrons: depends on element group
    //   * Group 14 (C): 2 pi electrons (anion with lone pair)
    //   * Group 15 (N, P, As): 1 pi electron (one in-plane lone pair + one p electron)
    // - 4+ electrons: 2 pi electrons (one lone pair participates, rest stay in-plane)

    let group_number = element.valence_electrons();

    match electrons_after_sigma {
        0 => Some(0),
        1 => Some(1),
        2 => Some(2),
        3 => {
            // For group 14 (C with 4 valence electrons), 3 electrons means anion → 2 pi
            // For group 15 (N, P, As with 5 valence electrons), 3 electrons → 1 pi (pyridine-type)
            if group_number == 4 {
                Some(2)
            } else {
                Some(1)
            }
        }
        _ => Some(2), // 4+ electrons: max 2 can participate in pi system
    }
}

/// Checks whether a pi electron count satisfies Hückel's rule (4n+2).
pub(crate) fn satisfies_huckel(pi_electrons: u8) -> bool {
    pi_electrons >= 2 && (pi_electrons - 2) % 4 == 0
}

/// Validates aromaticity of all aromatic rings in a molecule.
///
/// Returns a vector of `AromaticityCheck` results, one per aromatic ring.
pub fn validate_aromaticity(molecule: &Molecule) -> Vec<AromaticityCheck> {
    let rings = molecule.aromatic_rings();
    let mut results = Vec::with_capacity(rings.len());

    for ring in rings {
        let mut total_pi: u8 = 0;
        let mut determinable = true;

        for &node_idx in &ring.nodes {
            match pi_electron_contribution(molecule, node_idx) {
                Some(contrib) => {
                    total_pi = total_pi.saturating_add(contrib);
                }
                None => {
                    determinable = false;
                    break;
                }
            }
        }

        if determinable {
            results.push(AromaticityCheck {
                ring,
                pi_electrons: Some(total_pi),
                is_valid: satisfies_huckel(total_pi),
            });
        } else {
            results.push(AromaticityCheck {
                ring,
                pi_electrons: None,
                is_valid: true, // undeterminable -> skip validation
            });
        }
    }

    results
}

/// Validates that all aromatic rings in the molecule satisfy Hückel's rule.
///
/// Returns `Ok(())` if all rings are valid or undeterminable.
/// Returns `Err(MoleculeError::HuckelViolation)` on the first invalid ring.
pub fn require_valid_aromaticity(molecule: &Molecule) -> Result<(), MoleculeError> {
    aromatic_kekule_bonds(molecule)?;
    Ok(())
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum AromaticDoubleRole {
    Required,
    Forbidden,
    Flexible,
}

/// Returns the aromatic bond indices that must be rendered as double bonds in
/// a valence-compatible Kekule assignment.
pub(crate) fn aromatic_kekule_bonds(molecule: &Molecule) -> Result<HashSet<usize>, MoleculeError> {
    let atom_count = molecule.nodes().len();
    let mut adjacency = vec![Vec::<(NodeIndex, usize)>::new(); atom_count];

    for (bond_idx, bond) in molecule.bonds().iter().enumerate() {
        if bond.kind() != BondType::Aromatic {
            continue;
        }
        if !molecule.nodes()[bond.source() as usize].aromatic()
            || !molecule.nodes()[bond.target() as usize].aromatic()
        {
            return Err(MoleculeError::InvalidAromaticBond {
                atom1: bond.source(),
                atom2: bond.target(),
            });
        }
        adjacency[bond.source() as usize].push((bond.target(), bond_idx));
        adjacency[bond.target() as usize].push((bond.source(), bond_idx));
    }

    let rings = molecule.aromatic_rings();
    let ring_atoms = rings
        .iter()
        .flat_map(|ring| ring.nodes.iter().copied())
        .collect::<HashSet<_>>();
    for (atom_idx, node) in molecule.nodes().iter().enumerate() {
        if node.aromatic() && !ring_atoms.contains(&(atom_idx as NodeIndex)) {
            return Err(MoleculeError::AromaticAtomOutsideRing {
                atom: atom_idx as NodeIndex,
            });
        }
    }

    for check in validate_aromaticity(molecule) {
        if !check.is_valid {
            return Err(MoleculeError::HuckelViolation {
                ring: check.ring.nodes,
                pi_electrons: check.pi_electrons.unwrap_or(0),
            });
        }
    }

    let roles = (0..atom_count)
        .map(|atom_idx| {
            if !molecule.nodes()[atom_idx].aromatic() {
                return AromaticDoubleRole::Forbidden;
            }
            match pi_electron_contribution(molecule, atom_idx as NodeIndex) {
                Some(1) => AromaticDoubleRole::Required,
                Some(_) => AromaticDoubleRole::Forbidden,
                None => AromaticDoubleRole::Flexible,
            }
        })
        .collect::<Vec<_>>();

    let mut used_atoms = vec![false; atom_count];
    for bond in molecule.bonds() {
        if bond.kind() != BondType::Double
            || !molecule.nodes()[bond.source() as usize].aromatic()
            || !molecule.nodes()[bond.target() as usize].aromatic()
        {
            continue;
        }

        let source = bond.source() as usize;
        let target = bond.target() as usize;
        if roles[source] == AromaticDoubleRole::Forbidden
            || roles[target] == AromaticDoubleRole::Forbidden
            || used_atoms[source]
            || used_atoms[target]
        {
            return Err(MoleculeError::AromaticKekulizationFailed {
                atoms: vec![bond.source(), bond.target()],
            });
        }
        used_atoms[source] = true;
        used_atoms[target] = true;
    }

    let mut visited = vec![false; atom_count];
    let mut selected_bonds = HashSet::new();
    for start in 0..atom_count {
        if visited[start] || adjacency[start].is_empty() {
            continue;
        }
        let mut queue = VecDeque::from([start as NodeIndex]);
        let mut component = Vec::new();
        visited[start] = true;
        while let Some(atom) = queue.pop_front() {
            component.push(atom);
            for &(neighbor, _) in &adjacency[atom as usize] {
                if !visited[neighbor as usize] {
                    visited[neighbor as usize] = true;
                    queue.push_back(neighbor);
                }
            }
        }

        let mut component_bonds = Vec::new();
        if !match_required_aromatic_atoms(
            &component,
            &adjacency,
            &roles,
            &mut used_atoms,
            &mut component_bonds,
        ) {
            return Err(MoleculeError::AromaticKekulizationFailed { atoms: component });
        }
        selected_bonds.extend(component_bonds);
    }

    let unmatched = (0..atom_count)
        .filter(|&atom| roles[atom] == AromaticDoubleRole::Required && !used_atoms[atom])
        .map(|atom| atom as NodeIndex)
        .collect::<Vec<_>>();
    if !unmatched.is_empty() {
        return Err(MoleculeError::AromaticKekulizationFailed { atoms: unmatched });
    }

    Ok(selected_bonds)
}

fn match_required_aromatic_atoms(
    component: &[NodeIndex],
    adjacency: &[Vec<(NodeIndex, usize)>],
    roles: &[AromaticDoubleRole],
    used_atoms: &mut [bool],
    selected_bonds: &mut Vec<usize>,
) -> bool {
    let Some(atom) = component.iter().copied().find(|atom| {
        roles[*atom as usize] == AromaticDoubleRole::Required && !used_atoms[*atom as usize]
    }) else {
        return true;
    };

    let mut choices = adjacency[atom as usize].clone();
    choices.sort_unstable_by_key(|(neighbor, bond_idx)| {
        (
            roles[*neighbor as usize] != AromaticDoubleRole::Required,
            *bond_idx,
        )
    });
    for (neighbor, bond_idx) in choices {
        if used_atoms[neighbor as usize]
            || roles[neighbor as usize] == AromaticDoubleRole::Forbidden
        {
            continue;
        }
        used_atoms[atom as usize] = true;
        used_atoms[neighbor as usize] = true;
        selected_bonds.push(bond_idx);
        if match_required_aromatic_atoms(component, adjacency, roles, used_atoms, selected_bonds) {
            return true;
        }
        selected_bonds.pop();
        used_atoms[atom as usize] = false;
        used_atoms[neighbor as usize] = false;
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parse;

    #[test]
    fn benzene_is_valid() {
        let mol = parse("c1ccccc1").unwrap();
        let checks = validate_aromaticity(&mol);
        assert_eq!(checks.len(), 1);
        assert_eq!(checks[0].pi_electrons, Some(6));
        assert!(checks[0].is_valid);
    }

    #[test]
    fn pyridine_is_valid() {
        let mol = parse("c1ccncc1").unwrap();
        let checks = validate_aromaticity(&mol);
        assert_eq!(checks.len(), 1);
        assert_eq!(checks[0].pi_electrons, Some(6));
        assert!(checks[0].is_valid);
    }

    #[test]
    fn pyrrole_is_valid() {
        let mol = parse("c1cc[nH]c1").unwrap();
        let checks = validate_aromaticity(&mol);
        assert_eq!(checks.len(), 1);
        assert_eq!(checks[0].pi_electrons, Some(6));
        assert!(checks[0].is_valid);
    }

    #[test]
    fn furan_is_valid() {
        let mol = parse("c1ccoc1").unwrap();
        let checks = validate_aromaticity(&mol);
        assert_eq!(checks.len(), 1);
        assert_eq!(checks[0].pi_electrons, Some(6));
        assert!(checks[0].is_valid);
    }

    #[test]
    fn thiophene_is_valid() {
        let mol = parse("c1ccsc1").unwrap();
        let checks = validate_aromaticity(&mol);
        assert_eq!(checks.len(), 1);
        assert_eq!(checks[0].pi_electrons, Some(6));
        assert!(checks[0].is_valid);
    }

    #[test]
    fn imidazole_is_valid() {
        let mol = parse("c1cnc[nH]1").unwrap();
        let checks = validate_aromaticity(&mol);
        assert_eq!(checks.len(), 1);
        assert_eq!(checks[0].pi_electrons, Some(6));
        assert!(checks[0].is_valid);
    }

    #[test]
    fn cyclopentadienyl_anion_is_valid() {
        let mol = parse("[c-]1cccc1").unwrap();
        let checks = validate_aromaticity(&mol);
        assert_eq!(checks.len(), 1);
        assert_eq!(checks[0].pi_electrons, Some(6));
        assert!(checks[0].is_valid);
    }

    #[test]
    fn huckel_validation_passes_for_benzene() {
        let mol = parse("c1ccccc1").unwrap();
        assert!(require_valid_aromaticity(&mol).is_ok());
    }

    #[test]
    fn satisfies_huckel_values() {
        assert!(!satisfies_huckel(0));
        assert!(!satisfies_huckel(1));
        assert!(satisfies_huckel(2));
        assert!(!satisfies_huckel(3));
        assert!(!satisfies_huckel(4));
        assert!(!satisfies_huckel(5));
        assert!(satisfies_huckel(6));
        assert!(!satisfies_huckel(7));
        assert!(!satisfies_huckel(8));
        assert!(!satisfies_huckel(9));
        assert!(satisfies_huckel(10));
        assert!(satisfies_huckel(14));
        assert!(satisfies_huckel(18));
        assert!(satisfies_huckel(22));
    }

    #[test]
    fn non_aromatic_molecule_returns_empty() {
        let mol = parse("CCCCCC").unwrap();
        let checks = validate_aromaticity(&mol);
        assert!(checks.is_empty());
    }

    #[test]
    fn naphthalene_both_rings_valid() {
        let mol = parse("c1ccc2ccccc2c1").unwrap();
        let checks = validate_aromaticity(&mol);
        assert_eq!(checks.len(), 2);
        for check in &checks {
            assert!(check.is_valid);
        }
    }
}
