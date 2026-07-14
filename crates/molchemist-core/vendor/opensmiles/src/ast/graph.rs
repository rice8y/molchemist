use super::bond::BondType;
use super::molecule::Molecule;
use crate::NodeIndex;
use std::collections::{HashSet, VecDeque};

/// A ring represented by the ordered sequence of node indices.
#[derive(Debug, Clone, PartialEq)]
pub struct Ring {
    pub nodes: Vec<NodeIndex>,
}

impl Ring {
    pub fn size(&self) -> usize {
        self.nodes.len()
    }
}

/// Finds all aromatic rings in a molecule.
///
/// Uses the shortest-cycle-per-edge approach:
/// 1. Build a subgraph of only aromatic nodes connected by aromatic bonds
/// 2. For each edge, find the shortest cycle containing it (by removing the edge
///    and finding the shortest path between its endpoints)
/// 3. Deduplicate rings by their node sets
pub fn find_aromatic_rings(molecule: &Molecule) -> Vec<Ring> {
    let n = molecule.nodes().len();
    if n == 0 {
        return Vec::new();
    }

    // Build aromatic subgraph adjacency list
    let mut adj: Vec<Vec<NodeIndex>> = vec![Vec::new(); n];
    let mut is_aromatic = vec![false; n];
    let mut aromatic_edges: HashSet<(NodeIndex, NodeIndex)> = HashSet::new();

    for (i, node) in molecule.nodes().iter().enumerate() {
        if node.aromatic() {
            is_aromatic[i] = true;
        }
    }

    for bond in molecule.bonds() {
        let s = bond.source() as usize;
        let t = bond.target() as usize;
        if is_aromatic[s] && is_aromatic[t] && bond.kind() == BondType::Aromatic {
            adj[s].push(bond.target());
            adj[t].push(bond.source());
            let edge = if bond.source() < bond.target() {
                (bond.source(), bond.target())
            } else {
                (bond.target(), bond.source())
            };
            aromatic_edges.insert(edge);
        }
    }

    find_rings_in_subgraph(&adj, &aromatic_edges, n)
}

/// Finds all minimal rings in a subgraph defined by an adjacency list and edge set.
/// For each edge, finds the shortest cycle containing it (BFS), then deduplicates.
/// Shared by `find_aromatic_rings` and Kekulé ring detection in display.
pub(crate) fn find_rings_in_subgraph(
    adj: &[Vec<NodeIndex>],
    edges: &HashSet<(NodeIndex, NodeIndex)>,
    n: usize,
) -> Vec<Ring> {
    let mut rings: Vec<Vec<NodeIndex>> = Vec::new();

    for &(u, v) in edges {
        if let Some(path) = shortest_path_excluding_edge(u, v, adj, n) {
            let mut sorted = path.clone();
            sorted.sort();
            let is_duplicate = rings.iter().any(|existing| {
                let mut existing_sorted = existing.clone();
                existing_sorted.sort();
                existing_sorted == sorted
            });
            if !is_duplicate {
                rings.push(path);
            }
        }
    }

    rings.sort_by_key(|r| r.len());
    rings.into_iter().map(|nodes| Ring { nodes }).collect()
}

/// BFS from `u` to `v` without using the direct edge (u, v).
/// Returns the ring as the path from u to v (which, combined with the
/// excluded edge, forms a cycle).
pub(crate) fn shortest_path_excluding_edge(
    u: NodeIndex,
    v: NodeIndex,
    adj: &[Vec<NodeIndex>],
    n: usize,
) -> Option<Vec<NodeIndex>> {
    let mut visited = vec![false; n];
    let mut parent: Vec<Option<NodeIndex>> = vec![None; n];
    let mut queue = VecDeque::new();

    visited[u as usize] = true;
    queue.push_back(u);

    let mut found = false;

    'outer: while let Some(curr) = queue.pop_front() {
        for &next in &adj[curr as usize] {
            // Skip the direct edge u→v (only from u's side)
            if curr == u && next == v {
                continue;
            }

            if !visited[next as usize] {
                visited[next as usize] = true;
                parent[next as usize] = Some(curr);
                if next == v {
                    found = true;
                    break 'outer;
                }
                queue.push_back(next);
            }
        }
    }

    if !found {
        return None;
    }

    // Reconstruct path from u to v
    let mut path = Vec::new();
    let mut node = v;
    loop {
        path.push(node);
        if node == u {
            break;
        }
        node = parent[node as usize]?;
    }
    path.reverse();
    Some(path)
}

impl Molecule {
    /// Returns all aromatic rings found in this molecule.
    pub fn aromatic_rings(&self) -> Vec<Ring> {
        find_aromatic_rings(self)
    }
}

#[cfg(test)]
mod tests {
    use crate::parse;

    #[test]
    fn benzene_has_one_ring_of_size_6() {
        let mol = parse("c1ccccc1").unwrap();
        let rings = mol.aromatic_rings();
        assert_eq!(rings.len(), 1);
        assert_eq!(rings[0].size(), 6);
    }

    #[test]
    fn naphthalene_has_two_rings() {
        let mol = parse("c1ccc2ccccc2c1").unwrap();
        let rings = mol.aromatic_rings();
        assert_eq!(rings.len(), 2);
        for ring in &rings {
            assert!(ring.size() == 6, "expected size 6, got {}", ring.size());
        }
    }

    #[test]
    fn biphenyl_has_two_separate_rings() {
        let mol = parse("c1ccccc1c1ccccc1").unwrap();
        let rings = mol.aromatic_rings();
        assert_eq!(rings.len(), 2);
        for ring in &rings {
            assert_eq!(ring.size(), 6);
        }
    }

    #[test]
    fn cyclohexane_has_no_aromatic_rings() {
        let mol = parse("C1CCCCC1").unwrap();
        let rings = mol.aromatic_rings();
        assert_eq!(rings.len(), 0);
    }

    #[test]
    fn toluene_has_one_ring() {
        let mol = parse("Cc1ccccc1").unwrap();
        let rings = mol.aromatic_rings();
        assert_eq!(rings.len(), 1);
        assert_eq!(rings[0].size(), 6);
    }

    #[test]
    fn indole_has_two_rings() {
        let mol = parse("c1ccc2[nH]ccc2c1").unwrap();
        let rings = mol.aromatic_rings();
        assert_eq!(rings.len(), 2);
        let mut sizes: Vec<usize> = rings.iter().map(|r| r.size()).collect();
        sizes.sort();
        assert_eq!(sizes, vec![5, 6]);
    }

    #[test]
    fn pyrrole_has_one_ring_of_size_5() {
        let mol = parse("c1cc[nH]c1").unwrap();
        let rings = mol.aromatic_rings();
        assert_eq!(rings.len(), 1);
        assert_eq!(rings[0].size(), 5);
    }

    #[test]
    fn empty_molecule_has_no_rings() {
        let mol = parse("C").unwrap();
        let rings = mol.aromatic_rings();
        assert_eq!(rings.len(), 0);
    }
}
