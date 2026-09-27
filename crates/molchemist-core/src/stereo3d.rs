//! Preserve oriented tetrahedral volume when projecting a 3D CTAB into 2D.
//! Wedges encode local orientation; explicit stereochemical annotations take
//! precedence over orientation inferred from coordinates.
use crate::{ChemicalBondOrder, ChemicalRecord};
use std::collections::{BTreeMap, BTreeSet};

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn volume(p: &[[f64; 3]; 4]) -> f64 {
    let a = sub(p[0], p[3]);
    let b = sub(p[1], p[3]);
    let c = sub(p[2], p[3]);
    a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
}

fn ligand_classes(record: &ChemicalRecord, center: usize) -> Vec<usize> {
    let initial = record
        .atoms
        .iter()
        .map(|a| format!("{}:{:?}:{}", a.element, a.isotope, a.formal_charge))
        .collect::<Vec<_>>();
    let rank = |keys: &[String]| {
        let dictionary = keys
            .iter()
            .cloned()
            .collect::<BTreeSet<_>>()
            .into_iter()
            .enumerate()
            .map(|(i, k)| (k, i))
            .collect::<BTreeMap<_, _>>();
        keys.iter().map(|k| dictionary[k]).collect::<Vec<_>>()
    };
    let mut colors = rank(&initial);
    for _ in 0..record.atoms.len() {
        let keys = record
            .atoms
            .iter()
            .enumerate()
            .map(|(i, _)| {
                let mut neighbors = Vec::new();
                for b in &record.bonds {
                    let other = if b.atom1_index == i {
                        b.atom2_index
                    } else if b.atom2_index == i {
                        b.atom1_index
                    } else {
                        continue;
                    };
                    if i != center && other != center {
                        neighbors.push(format!("{:?}:{}", b.order, colors[other]));
                    }
                }
                neighbors.sort();
                format!("{}:{}:{}", initial[i], colors[i], neighbors.join(","))
            })
            .collect::<Vec<_>>();
        let next = rank(&keys);
        let stable = (0..colors.len())
            .all(|i| (0..colors.len()).all(|j| (colors[i] == colors[j]) == (next[i] == next[j])));
        colors = next;
        if stable {
            break;
        }
    }
    colors
}

pub fn tetrahedral_wedges(
    record: &ChemicalRecord,
    xy: &[[f64; 2]],
) -> Result<Vec<(usize, usize, u8)>, String> {
    oriented_wedges(record, xy, None)
}

pub fn reorient_explicit_center(
    record: &ChemicalRecord,
    xy: &[[f64; 2]],
    center: usize,
) -> Result<(usize, usize, u8), String> {
    oriented_wedges(record, xy, Some(center))?.into_iter().next().ok_or_else(|| format!("Cannot preserve explicit wedge at atom {} during relayout: local geometry is degenerate or unsupported", record.atoms[center].source_id))
}

fn oriented_wedges(
    record: &ChemicalRecord,
    xy: &[[f64; 2]],
    explicit_center: Option<usize>,
) -> Result<Vec<(usize, usize, u8)>, String> {
    if xy.len() != record.atoms.len() {
        return Err("3D stereo projection coordinate count mismatch".into());
    }
    let mut output = Vec::new();
    let mut used = BTreeSet::new();
    for (center, atom) in record.atoms.iter().enumerate() {
        if explicit_center.is_some_and(|c| c != center) {
            continue;
        }
        let incident = record
            .bonds
            .iter()
            .filter_map(|b| {
                if b.atom1_index == center {
                    Some((b, b.atom2_index))
                } else if b.atom2_index == center {
                    Some((b, b.atom1_index))
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();
        if !(3..=4).contains(&incident.len())
            || incident
                .iter()
                .any(|(b, _)| b.order != ChemicalBondOrder::Single)
        {
            continue;
        }
        // Three ligands imply a fourth hydrogen only at a neutral tetravalent
        // carbon/silicon. Other inferred centers require four explicit ligands.
        if incident.len() == 3
            && explicit_center.is_none()
            && (!matches!(atom.element.as_str(), "C" | "Si") || atom.formal_charge != 0)
        {
            continue;
        }
        let classes = ligand_classes(record, center);
        if explicit_center.is_none()
            && incident
                .iter()
                .map(|(_, i)| classes[*i])
                .collect::<BTreeSet<_>>()
                .len()
                != incident.len()
        {
            continue;
        }
        if incident.len() == 3
            && explicit_center.is_none()
            && incident
                .iter()
                .any(|(_, i)| record.atoms[*i].element == "H")
        {
            continue;
        }
        let mut original = [[0.0; 3]; 4];
        let mut planar = [[0.0; 3]; 4];
        for (i, (_, neighbor)) in incident.iter().enumerate() {
            original[i] = sub(record.atoms[*neighbor].coordinates, atom.coordinates);
            planar[i] = [
                xy[*neighbor][0] - xy[center][0],
                xy[*neighbor][1] - xy[center][1],
                0.0,
            ];
        }
        if incident.len() == 3 {
            for axis in 0..3 {
                original[3][axis] = -(original[0][axis] + original[1][axis] + original[2][axis]);
                planar[3][axis] = -(planar[0][axis] + planar[1][axis] + planar[2][axis]);
            }
        }
        let desired = volume(&original);
        let size = original
            .iter()
            .flatten()
            .fold(0.0_f64, |s, v| s.max(v.abs()));
        if desired.abs() <= size.powi(3) * 1e-8 {
            continue;
        }
        let mut candidates = Vec::new();
        for (i, (bond, neighbor)) in incident.iter().enumerate() {
            if used.contains(&bond.index) {
                continue;
            }
            let mut candidate = planar;
            candidate[i][2] = 1.0;
            if incident.len() == 3 {
                candidate[3][2] = -1.0;
            }
            let projected = volume(&candidate);
            if projected.abs() > 1e-10 {
                candidates.push((
                    projected.abs(),
                    *neighbor,
                    bond.index,
                    if desired * projected > 0.0 { 1 } else { 6 },
                ));
            }
        }
        candidates.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.2.cmp(&b.2)));
        if let Some(&(_, neighbor, bond, style)) = candidates.first() {
            used.insert(bond);
            output.push((center, neighbor, style));
        } else {
            return Err(format!(
                "Cannot project 3D stereocenter {} onto collinear 2D ligands; request reflow",
                atom.source_id
            ));
        }
    }
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn mirror_reverses_wedge_and_planar_input_is_not_invented() {
        let mol = "tetra\n\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 5 4 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 2 F 1 1 1 0\nM  V30 3 Cl -1 -1 1 0\nM  V30 4 Br -1 1 -1 0\nM  V30 5 I 1 -1 -1 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 2 1 1 3\nM  V30 3 1 1 4\nM  V30 4 1 1 5\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n";
        let mut record = crate::inspect_sdf_record(mol, 1).unwrap();
        let xy = [[0., 0.], [1., 0.], [0., 1.], [-1., 0.], [0., -1.]];
        let first = tetrahedral_wedges(&record, &xy).unwrap();
        assert_eq!(first.len(), 1);
        for a in &mut record.atoms {
            a.coordinates[2] *= -1.;
        }
        let mirror = tetrahedral_wedges(&record, &xy).unwrap();
        assert_ne!(first[0].2, mirror[0].2);
        for a in &mut record.atoms {
            a.coordinates[2] = 0.;
        }
        assert!(tetrahedral_wedges(&record, &xy).unwrap().is_empty());
    }
}
