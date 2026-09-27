//! Implicit hydrogens for unambiguous CTfile atom environments.
use crate::{ChemicalBondOrder, ChemicalFormat, ChemicalRecord};

pub(crate) fn implicit_hydrogens(record: &ChemicalRecord, index: usize) -> u8 {
    let atom = &record.atoms[index];
    // An attachment site has an external bond whose order is not in this CTAB.
    if !atom.attachment_points.is_empty() || !atom.rgroup_labels.is_empty() {
        return 0;
    }
    if atom.query.as_ref().is_some_and(|q| {
        q.elements.is_some()
            || q.hydrogen_count.is_some()
            || q.substitution_count.is_some()
            || q.unsaturated.is_some()
            || q.ring_bond_count.is_some()
    }) {
        return 0;
    }
    let mut bonded = 0u16;
    for bond in record
        .bonds
        .iter()
        .filter(|b| b.atom1_index == index || b.atom2_index == index)
    {
        if !bond.endpoint_source_ids.is_empty() || bond.attachment_mode.is_some() {
            return 0;
        }
        bonded += match bond.order {
            ChemicalBondOrder::Single => 1,
            ChemicalBondOrder::Double => 2,
            ChemicalBondOrder::Triple => 3,
            // Aromatic and query orders need a valence assignment before an H
            // count can be chosen. Do not turn an unknown count into a label.
            _ => return 0,
        };
    }
    let radical = match atom.radical.unwrap_or(0) {
        0 => 0,
        2 => 1,
        1 | 3 => 2,
        _ => return 0,
    };
    let occupied = bonded + radical;
    let explicit = atom
        .query
        .as_ref()
        .and_then(|q| q.valence)
        .or(atom.valence.map(|v| v as i8))
        .unwrap_or(0);
    let target =
        if explicit == -1 || (explicit == 15 && record.format == ChemicalFormat::MolfileV2000) {
            0
        } else if explicit > 0 {
            explicit as u16
        } else {
            let defaults: &[u16] = match (atom.element.as_str(), atom.formal_charge) {
                ("B", 0) => &[3],
                ("B", -1) | ("C" | "Si", 0) | ("N", 1) => &[4],
                ("C", -1 | 1) | ("O", 1) => &[3],
                ("N", 0) | ("P", 0) => &[3, 5],
                ("N", -1) | ("O", 0) => &[2],
                ("O", -1) | ("F" | "Cl" | "Br" | "I", 0) => &[1],
                ("S", 0) => &[2, 4, 6],
                _ => return 0,
            };
            match defaults.iter().find(|&&v| v >= occupied) {
                Some(&v) => v,
                None => return 0,
            }
        };
    target.saturating_sub(occupied).min(u8::MAX as u16) as u8
}

#[cfg(test)]
mod tests {
    use super::*;
    fn oxygen(attributes: &str, order: u8, extra_h: bool) -> ChemicalRecord {
        let (count, bonds, hydrogen) = if extra_h {
            (3, 2, "M  V30 3 H 2 0 0 0\n")
        } else {
            (2, 1, "")
        };
        crate::inspect_sdf_record(&format!("valence\n\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS {count} {bonds} 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 2 O 1 0 0 0 {attributes}\n{hydrogen}M  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 {order} 1 2\n{}M  V30 END BOND\nM  V30 END CTAB\nM  END\n", if extra_h { "M  V30 2 1 2 3\n" } else { "" }), 1).unwrap()
    }
    #[test]
    fn oxygen_hydrogens_respect_charge_valence_queries_and_explicit_atoms() {
        for (attributes, order, extra_h, expected) in [
            ("", 1, false, 1),
            ("", 2, false, 0),
            ("CHG=-1", 1, false, 0),
            ("", 1, true, 0),
            ("VAL=-1", 1, false, 0),
            ("VAL=3", 1, false, 2),
            ("HCOUNT=1", 1, false, 0),
            ("HCOUNT=-1", 1, false, 0),
            ("RAD=2", 1, false, 0),
            ("", 8, false, 0),
            ("ATTCHPT=1", 1, false, 0),
        ] {
            assert_eq!(
                implicit_hydrogens(&oxygen(attributes, order, extra_h), 1),
                expected,
                "{attributes}, order {order}"
            );
        }
    }
}
