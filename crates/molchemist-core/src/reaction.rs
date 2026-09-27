//! Deterministic atom correspondence and graph-edit reaction centers.
//! Results report search completion and tied optima so callers can distinguish
//! a unique correspondence from a provisional assignment.
use crate::{MoleculeInput, ReactionInput};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ReactionAtom {
    pub molecule: usize,
    pub atom: usize,
    pub element: String,
    pub isotope: Option<u16>,
    pub map: Option<u32>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct AtomCorrespondence {
    pub reactant: ReactionAtom,
    pub product: ReactionAtom,
    pub map: u32,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct BondChange {
    pub maps: [u32; 2],
    pub before: Option<String>,
    pub after: Option<String>,
    pub reactant: Option<[usize; 3]>,
    pub product: Option<[usize; 3]>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ReactionAnalysis {
    pub reaction: ReactionInput,
    pub mapping: Vec<AtomCorrespondence>,
    pub changes: Vec<BondChange>,
    pub reactant_atoms: Vec<ReactionAtom>,
    pub product_atoms: Vec<ReactionAtom>,
    pub ambiguous: bool,
    pub search_complete: bool,
    pub unmapped_reactants: Vec<ReactionAtom>,
    pub unmapped_products: Vec<ReactionAtom>,
}

type Edges = BTreeMap<(usize, usize), String>;
fn side_graph(molecules: &[MoleculeInput]) -> Result<(Vec<ReactionAtom>, Edges), String> {
    let mut atoms = Vec::new();
    let mut edges = Edges::new();
    let mut maps = BTreeSet::new();
    for (molecule, input) in molecules.iter().enumerate() {
        let offset = atoms.len();
        if input.format == "smiles" {
            let graph = crate::parse(&input.data).map_err(|e| e.to_string())?;
            for (atom, node) in graph.nodes().iter().enumerate() {
                atoms.push(ReactionAtom {
                    molecule,
                    atom,
                    element: node.atom().element().to_string(),
                    isotope: node.atom().isotope(),
                    map: node.class().filter(|&m| m > 0).map(u32::from),
                });
            }
            for bond in graph.bonds() {
                use crate::BondType::*;
                let kind = match bond.kind() {
                    Disconnected => continue,
                    Simple | Up | Down => "single",
                    Double => "double",
                    Triple => "triple",
                    Quadruple => "quadruple",
                    Aromatic => "aromatic",
                };
                let a = bond.source() as usize + offset;
                let b = bond.target() as usize + offset;
                edges.insert((a.min(b), a.max(b)), kind.into());
            }
        } else {
            let record = crate::inspect_sdf_record(&input.data, 1)?;
            for atom in &record.atoms {
                atoms.push(ReactionAtom {
                    molecule,
                    atom: atom.index,
                    element: atom.element.clone(),
                    isotope: atom.isotope,
                    map: atom.atom_map.filter(|&m| m > 0),
                });
            }
            for bond in &record.bonds {
                let a = bond.atom1_index + offset;
                let b = bond.atom2_index + offset;
                edges.insert(
                    (a.min(b), a.max(b)),
                    format!("{:?}", bond.order).to_lowercase(),
                );
            }
        }
    }
    for atom in &atoms {
        if let Some(map) = atom.map {
            if !maps.insert(map) {
                return Err(format!(
                    "Duplicate atom-map number {map} on one reaction side"
                ));
            }
        }
    }
    Ok((atoms, edges))
}

struct Search<'a> {
    left: &'a [ReactionAtom],
    right: &'a [ReactionAtom],
    le: &'a Edges,
    re: &'a Edges,
    candidates: Vec<Vec<usize>>,
    order: Vec<usize>,
    best: Vec<Option<usize>>,
    best_score: (usize, usize),
    visited: usize,
    limit: usize,
    ambiguous: bool,
}
impl Search<'_> {
    fn visit(
        &mut self,
        depth: usize,
        mapping: &mut [Option<usize>],
        used: &mut [bool],
        matched: usize,
    ) {
        if self.visited >= self.limit {
            return;
        }
        self.visited += 1;
        if matched + self.order.len() - depth < self.best_score.0 {
            return;
        }
        if depth == self.order.len() {
            let conserved = self
                .le
                .iter()
                .filter(|((a, b), kind)| {
                    if let (Some(x), Some(y)) = (mapping[*a], mapping[*b]) {
                        self.re.get(&(x.min(y), x.max(y))) == Some(*kind)
                    } else {
                        false
                    }
                })
                .count();
            let score = (matched, conserved);
            if score > self.best_score {
                self.best_score = score;
                self.best = mapping.to_vec();
                self.ambiguous = false;
            } else if score == self.best_score && self.best != mapping {
                self.ambiguous = true;
            }
            return;
        }
        let a = self.order[depth];
        for b in self.candidates[a].clone() {
            if used[b] {
                continue;
            }
            mapping[a] = Some(b);
            used[b] = true;
            self.visit(depth + 1, mapping, used, matched + 1);
            mapping[a] = None;
            used[b] = false;
        }
        // Explicit common map numbers are mandatory constraints.
        let mandatory = self.left[a]
            .map
            .is_some_and(|m| self.right.iter().any(|a| a.map == Some(m)));
        if !mandatory {
            self.visit(depth + 1, mapping, used, matched);
        }
    }
}

pub fn analyze_reaction(
    reaction: ReactionInput,
    infer: bool,
    search_limit: usize,
) -> Result<ReactionAnalysis, String> {
    if search_limit == 0 {
        return Err("Atom-mapping search limit must be positive".into());
    }
    let (left, le) = side_graph(&reaction.reactants)?;
    let (right, re) = side_graph(&reaction.products)?;
    for a in &left {
        if let Some(b) = a.map.and_then(|m| right.iter().find(|b| b.map == Some(m))) {
            if a.element != b.element || a.isotope != b.isotope {
                return Err("Mapped atoms have different elements or isotopes".into());
            }
        }
    }
    let candidates = left
        .iter()
        .map(|a| {
            right
                .iter()
                .enumerate()
                .filter_map(|(i, b)| {
                    if a.element != b.element || a.isotope != b.isotope {
                        return None;
                    }
                    let same = a.map.is_some() && a.map == b.map;
                    if !infer && !same {
                        return None;
                    }
                    if a.map
                        .is_some_and(|m| right.iter().any(|r| r.map == Some(m)))
                        && !same
                    {
                        return None;
                    }
                    if b.map.is_some_and(|m| left.iter().any(|l| l.map == Some(m))) && !same {
                        return None;
                    }
                    if matches!((a.map, b.map), (Some(x), Some(y)) if x != y) {
                        return None;
                    }
                    Some(i)
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let mut order = (0..left.len()).collect::<Vec<_>>();
    order.sort_by_key(|&a| (candidates[a].len(), a));
    let mut search = Search {
        left: &left,
        right: &right,
        le: &le,
        re: &re,
        candidates,
        order,
        best: vec![None; left.len()],
        best_score: (0, 0),
        visited: 0,
        limit: search_limit,
        ambiguous: false,
    };
    search.visit(
        0,
        &mut vec![None; left.len()],
        &mut vec![false; right.len()],
        0,
    );
    let mut reserved = left
        .iter()
        .chain(&right)
        .filter_map(|a| a.map)
        .collect::<BTreeSet<_>>();
    let mut next = 1u32;
    let mut allocate = || {
        while reserved.contains(&next) {
            next = next
                .checked_add(1)
                .ok_or("Atom-map identifier space exhausted")?;
        }
        reserved.insert(next);
        Ok::<_, String>(next)
    };
    let mut mapping = Vec::new();
    let mut lm = BTreeMap::new();
    let mut rm = BTreeMap::new();
    for (a, b) in search.best.iter().enumerate() {
        if let Some(b) = b {
            let map = match left[a].map.or(right[*b].map) {
                Some(map) => map,
                None => allocate()?,
            };
            lm.insert(a, map);
            rm.insert(*b, map);
            mapping.push(AtomCorrespondence {
                reactant: left[a].clone(),
                product: right[*b].clone(),
                map,
            });
        }
    }
    let unmapped_reactants = left
        .iter()
        .enumerate()
        .filter(|(i, _)| !lm.contains_key(i))
        .map(|(_, a)| a.clone())
        .collect();
    let unmapped_products = right
        .iter()
        .enumerate()
        .filter(|(i, _)| !rm.contains_key(i))
        .map(|(_, a)| a.clone())
        .collect();
    // Distinct identifiers retain bonds to entering and leaving groups in the
    // reaction graph difference.
    for (atoms, maps) in [(&left, &mut lm), (&right, &mut rm)] {
        for (i, atom) in atoms.iter().enumerate() {
            if let std::collections::btree_map::Entry::Vacant(entry) = maps.entry(i) {
                entry.insert(match atom.map {
                    Some(map) => map,
                    None => allocate()?,
                });
            }
        }
    }
    let annotated_atoms = |atoms: &[ReactionAtom], maps: &BTreeMap<usize, u32>| {
        atoms
            .iter()
            .enumerate()
            .map(|(i, a)| {
                let mut a = a.clone();
                a.map = Some(maps[&i]);
                a
            })
            .collect::<Vec<_>>()
    };
    let reactant_atoms = annotated_atoms(&left, &lm);
    let product_atoms = annotated_atoms(&right, &rm);
    let mapped_edges = |edges: &Edges, maps: &BTreeMap<usize, u32>| {
        edges
            .iter()
            .filter_map(|(&(a, b), k)| {
                let (x, y) = (*maps.get(&a)?, *maps.get(&b)?);
                Some(((x.min(y), x.max(y)), k.clone()))
            })
            .collect::<BTreeMap<_, _>>()
    };
    let before = mapped_edges(&le, &lm);
    let after = mapped_edges(&re, &rm);
    let keys = before
        .keys()
        .chain(after.keys())
        .copied()
        .collect::<BTreeSet<_>>();
    let location = |a: u32, b: u32, atoms: &[ReactionAtom], maps: &BTreeMap<usize, u32>| {
        let x = maps
            .iter()
            .find(|(_, m)| **m == a)
            .map(|(i, _)| &atoms[*i])?;
        let y = maps
            .iter()
            .find(|(_, m)| **m == b)
            .map(|(i, _)| &atoms[*i])?;
        (x.molecule == y.molecule).then_some([x.molecule, x.atom, y.atom])
    };
    let changes = keys
        .into_iter()
        .filter_map(|(a, b)| {
            let x = before.get(&(a, b)).cloned();
            let y = after.get(&(a, b)).cloned();
            (x != y).then(|| BondChange {
                maps: [a, b],
                reactant: if x.is_some() {
                    location(a, b, &left, &lm)
                } else {
                    None
                },
                product: if y.is_some() {
                    location(a, b, &right, &rm)
                } else {
                    None
                },
                before: x,
                after: y,
            })
        })
        .collect();
    Ok(ReactionAnalysis {
        reaction,
        mapping,
        changes,
        reactant_atoms,
        product_atoms,
        ambiguous: search.ambiguous,
        search_complete: search.visited < search.limit,
        unmapped_reactants,
        unmapped_products,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn mapped_bond_order_change() {
        let r = crate::parse_reaction("[CH3:1][OH:2]>>[CH2:1]=[O:2]", "reaction-smiles").unwrap();
        let a = analyze_reaction(r, false, 10000).unwrap();
        assert_eq!(a.mapping.len(), 2);
        assert_eq!(a.changes.len(), 1);
        assert_eq!(a.changes[0].before.as_deref(), Some("single"));
        assert_eq!(a.changes[0].after.as_deref(), Some("double"));
        assert!(a.search_complete && !a.ambiguous);
    }
    #[test]
    fn symmetric_mapping_is_reported() {
        let r = crate::parse_reaction("CC>>CC", "reaction-smiles").unwrap();
        let a = analyze_reaction(r.clone(), true, 10000).unwrap();
        assert!(a.ambiguous);
        assert_eq!(a.mapping.len(), 2);
        assert!(!analyze_reaction(r, true, 1).unwrap().search_complete);
    }
    #[test]
    fn leaving_and_entering_groups_are_not_dropped() {
        let r = crate::parse_reaction("[CH3:1]Br>>[CH3:1]O", "reaction-smiles").unwrap();
        let a = analyze_reaction(r, false, 1000).unwrap();
        assert_eq!(a.changes.len(), 2);
        assert!(a.changes.iter().any(|c| c.before.is_none()));
        assert!(a.changes.iter().any(|c| c.after.is_none()));
    }

    #[test]
    fn substitution_preserves_the_mapped_carbon_backbone() {
        let r = crate::parse_reaction(
            "[CH3:1][CH2:2][CH2:3]Br>>[CH3:1][CH2:2][CH2:3]O",
            "reaction-smiles",
        )
        .unwrap();
        let a = analyze_reaction(r, false, 1000).unwrap();
        assert_eq!(a.mapping.len(), 3);
        assert_eq!(a.changes.len(), 2);
        assert!(a
            .changes
            .iter()
            .all(|c| c.maps.contains(&3) && !c.maps.contains(&1) && !c.maps.contains(&2)));
        let broken = a.changes.iter().find(|c| c.after.is_none()).unwrap();
        let formed = a.changes.iter().find(|c| c.before.is_none()).unwrap();
        assert_eq!(broken.reactant, Some([0, 2, 3]));
        assert_eq!(formed.product, Some([0, 2, 3]));
    }
}
