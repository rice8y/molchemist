use opensmiles::{parse as parse_smiles, BondType as SmilesBondType};
use sdfrust::{parse_sdf_string, BondOrder};
use serde::Serialize;
use std::collections::{HashMap, HashSet, VecDeque};
use wasm_minimal_protocol::*;

initiate_protocol!();

const LAYOUT_MAGIC: &[u8; 4] = b"MCG2";
const COORD_MAGIC: &[u8; 4] = b"MCC2";

#[derive(Serialize)]
#[serde(tag = "type")]
enum Command {
    #[serde(rename = "fragment")]
    Fragment {
        element: String,
        name: String,
        links: Vec<LinkData>,
        #[serde(skip_serializing_if = "Option::is_none")]
        annotation: Option<String>,
    },
    #[serde(rename = "bond")]
    Bond {
        #[serde(rename = "bondType")]
        bond_type: String,
        angle: f64,
        offset: Option<String>,
        #[serde(rename = "lengthScale")]
        length_scale: f64,
    },
    #[serde(rename = "branch")]
    Branch { body: Vec<Command> },
}

#[derive(Serialize)]
struct LinkData {
    target: String,
    #[serde(rename = "bondType")]
    bond_type: String,
    angle: f64,
    offset: Option<String>,
    #[serde(rename = "lengthScale")]
    length_scale: f64,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum BondKind {
    Single,
    Double,
    Triple,
}

#[derive(Clone)]
struct RenderAtom {
    element: String,
    x: f64,
    y: f64,
    hydrogens: u8,
    charge: i8,
    stereo_annotation: Option<String>,
}

#[derive(Clone)]
struct RenderBond {
    atom1: usize,
    atom2: usize,
    kind: BondKind,
}

struct RenderMol {
    atoms: Vec<RenderAtom>,
    bonds: Vec<RenderBond>,
}

#[derive(Clone, Copy)]
enum AtomStereoDirection {
    Clockwise,
    CounterClockwise,
}

#[derive(Clone)]
struct SmilesAtom {
    element: String,
    atomic_number: u8,
    hydrogens: u8,
    charge: i8,
    stereo: Option<AtomStereoSpec>,
    stereo_annotation: Option<String>,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum SmilesBondKind {
    Single,
    Double,
    Triple,
    Aromatic,
}

#[derive(Clone)]
struct SmilesBond {
    atom1: usize,
    atom2: usize,
    kind: SmilesBondKind,
}

#[derive(Clone)]
struct SmilesGraph {
    atoms: Vec<SmilesAtom>,
    bonds: Vec<SmilesBond>,
    double_bond_stereo: Vec<DoubleBondStereoSpec>,
}

#[derive(Clone, Copy)]
struct AtomStereoSpec {
    looking_from: Option<usize>,
    atom1: Option<usize>,
    atom2: Option<usize>,
    direction: AtomStereoDirection,
}

#[derive(Clone)]
struct DoubleBondStereoSpec {
    bond_index: usize,
    atom1: usize,
    atom2: usize,
    is_z: bool,
}

struct CoordPayload {
    coords: Vec<(f32, f32)>,
    bond_styles: Vec<u8>,
}

#[wasm_func]
pub fn sdf_to_ast(sdf_data: &[u8], options: &[u8]) -> Result<Vec<u8>, String> {
    let sdf_str = std::str::from_utf8(sdf_data).map_err(|e| e.to_string())?;
    let mol = parse_sdf_string(sdf_str).map_err(|e| e.to_string())?;
    let stereo_map = extract_stereo(sdf_str);
    let render_mol = render_mol_from_sdf(&mol);
    ast_from_render_mol(&render_mol, parse_mode(options), &stereo_map)
}

#[wasm_func]
pub fn smiles_to_layout_input(smiles_data: &[u8]) -> Result<Vec<u8>, String> {
    let smiles = std::str::from_utf8(smiles_data).map_err(|e| e.to_string())?;
    let graph = parse_smiles_graph(smiles)?;
    Ok(encode_layout_input(&graph))
}

#[wasm_func]
pub fn smiles_to_full_layout_input(smiles_data: &[u8]) -> Result<Vec<u8>, String> {
    let smiles = std::str::from_utf8(smiles_data).map_err(|e| e.to_string())?;
    let graph = parse_smiles_graph(smiles)?;
    let expanded = expand_smiles_graph_hydrogens(&graph);
    Ok(encode_layout_input(&expanded))
}

#[wasm_func]
pub fn smiles_to_ast(smiles_data: &[u8], coords_data: &[u8], options: &[u8]) -> Result<Vec<u8>, String> {
    let smiles = std::str::from_utf8(smiles_data).map_err(|e| e.to_string())?;
    let mode = parse_mode(options);
    let mut graph = parse_smiles_graph(smiles)?;
    if mode == "full" {
        graph = expand_smiles_graph_hydrogens(&graph);
    }
    let coords = decode_coords(coords_data, graph.atoms.len(), graph.bonds.len())?;
    let render_mol = render_mol_from_smiles(&graph, &coords.coords);
    let stereo_map = build_smiles_stereo_map(&graph, &coords.bond_styles);
    ast_from_render_mol(&render_mol, mode, &stereo_map)
}

fn parse_mode(options: &[u8]) -> &str {
    std::str::from_utf8(options).unwrap_or("full")
}

fn render_mol_from_sdf(mol: &sdfrust::Molecule) -> RenderMol {
    let atoms = mol
        .atoms
        .iter()
        .map(|atom| RenderAtom {
            element: atom.element.clone(),
            x: atom.x,
            y: atom.y,
            hydrogens: 0,
            charge: 0,
            stereo_annotation: None,
        })
        .collect();

    let bonds = mol
        .bonds
        .iter()
        .map(|bond| RenderBond {
            atom1: bond.atom1,
            atom2: bond.atom2,
            kind: match bond.order {
                BondOrder::Double => BondKind::Double,
                BondOrder::Triple => BondKind::Triple,
                _ => BondKind::Single,
            },
        })
        .collect();

    RenderMol { atoms, bonds }
}

fn render_mol_from_smiles(graph: &SmilesGraph, coords: &[(f32, f32)]) -> RenderMol {
    let atoms = graph
        .atoms
        .iter()
        .zip(coords.iter())
        .map(|(atom, &(x, y))| RenderAtom {
            element: atom.element.clone(),
            x: x as f64,
            y: y as f64,
            hydrogens: atom.hydrogens,
            charge: atom.charge,
            stereo_annotation: atom.stereo_annotation.clone(),
        })
        .collect::<Vec<_>>();

    let bonds = graph
        .bonds
        .iter()
        .map(|bond| RenderBond {
            atom1: bond.atom1,
            atom2: bond.atom2,
            kind: match bond.kind {
                SmilesBondKind::Single => BondKind::Single,
                SmilesBondKind::Double => BondKind::Double,
                SmilesBondKind::Triple => BondKind::Triple,
                SmilesBondKind::Aromatic => BondKind::Single,
            },
        })
        .collect::<Vec<_>>();

    RenderMol { atoms, bonds }
}

fn expand_smiles_graph_hydrogens(graph: &SmilesGraph) -> SmilesGraph {
    let mut expanded = graph.clone();
    for atom in &mut expanded.atoms {
        atom.hydrogens = 0;
    }

    for atom_idx in 0..graph.atoms.len() {
        let hydrogen_count = graph.atoms[atom_idx].hydrogens;
        for _ in 0..hydrogen_count {
            let hydrogen_idx = expanded.atoms.len();
            expanded.atoms.push(SmilesAtom {
                element: "H".to_string(),
                atomic_number: 1,
                hydrogens: 0,
                charge: 0,
                stereo: None,
                stereo_annotation: None,
            });
            expanded.bonds.push(SmilesBond {
                atom1: atom_idx,
                atom2: hydrogen_idx,
                kind: SmilesBondKind::Single,
            });
        }
    }

    expanded
}

fn parse_smiles_graph(smiles: &str) -> Result<SmilesGraph, String> {
    let molecule = parse_smiles(smiles).map_err(|e| e.to_string())?;
    let mut incident_bonds = vec![Vec::<(usize, usize)>::new(); molecule.nodes().len()];
    let mut original_to_graph = vec![None; molecule.bonds().len()];
    let mut bonds = Vec::new();

    for (bond_idx, bond) in molecule.bonds().iter().enumerate() {
        let source = bond.source() as usize;
        let target = bond.target() as usize;
        if bond.kind() == SmilesBondType::Disconnected {
            continue;
        }

        incident_bonds[source].push((bond_idx, target));
        incident_bonds[target].push((bond_idx, source));

        let kind = match bond.kind() {
            SmilesBondType::Simple | SmilesBondType::Up | SmilesBondType::Down => {
                SmilesBondKind::Single
            }
            SmilesBondType::Double => SmilesBondKind::Double,
            SmilesBondType::Triple => SmilesBondKind::Triple,
            SmilesBondType::Aromatic => SmilesBondKind::Aromatic,
            SmilesBondType::Disconnected => continue,
            SmilesBondType::Quadruple => {
                return Err("Quadruple bonds are not supported yet".to_string())
            }
        };

        original_to_graph[bond_idx] = Some(bonds.len());
        bonds.push(SmilesBond {
            atom1: source,
            atom2: target,
            kind,
        });
    }

    let atoms = molecule
        .nodes()
        .iter()
        .enumerate()
        .map(|(atom_idx, node)| {
            let atomic_number = node.atom().element().atomic_number();
            if atomic_number == 0 {
                return Err(format!(
                    "Unsupported atom '{}' for SMILES layout",
                    node.atom().element()
                ));
            }

            Ok(SmilesAtom {
                element: node.atom().element().to_string(),
                atomic_number,
                hydrogens: node.hydrogens(),
                charge: node.atom().charge(),
                stereo: atom_stereo_spec(atom_idx, node.chirality(), &incident_bonds),
                stereo_annotation: chirality_annotation(node.chirality()),
            })
        })
        .collect::<Result<Vec<_>, _>>()?;

    let double_bond_stereo = build_double_bond_stereo_specs(&molecule, &original_to_graph);

    kekulize_aromatic_bonds(atoms.len(), &mut bonds);

    Ok(SmilesGraph {
        atoms,
        bonds,
        double_bond_stereo,
    })
}

fn atom_stereo_spec(
    atom_idx: usize,
    chirality: Option<opensmiles::Chirality>,
    incident_bonds: &[Vec<(usize, usize)>],
) -> Option<AtomStereoSpec> {
    let direction = match chirality {
        Some(opensmiles::Chirality::TH1) => AtomStereoDirection::CounterClockwise,
        Some(opensmiles::Chirality::TH2) => AtomStereoDirection::Clockwise,
        _ => return None,
    };

    let neighbors = incident_bonds
        .get(atom_idx)?
        .iter()
        .map(|&(_, neighbor)| neighbor)
        .collect::<Vec<_>>();

    if neighbors.len() < 3 || neighbors.len() > 4 {
        return None;
    }

    Some(AtomStereoSpec {
        looking_from: neighbors.first().copied(),
        atom1: neighbors.get(1).copied(),
        atom2: neighbors.get(2).copied(),
        direction,
    })
}

fn chirality_annotation(chirality: Option<opensmiles::Chirality>) -> Option<String> {
    match chirality {
        Some(opensmiles::Chirality::TH1 | opensmiles::Chirality::TH2) | None => None,
        Some(chirality) => Some(chirality.to_string()),
    }
}

fn build_double_bond_stereo_specs(
    molecule: &opensmiles::Molecule,
    original_to_graph: &[Option<usize>],
) -> Vec<DoubleBondStereoSpec> {
    let mut specs = Vec::new();

    for (bond_idx, bond) in molecule.bonds().iter().enumerate() {
        if bond.kind() != SmilesBondType::Double {
            continue;
        }

        let left_center = bond.source() as usize;
        let right_center = bond.target() as usize;
        let Some(left_marker) = directional_marker_for_center(molecule, left_center, right_center)
        else {
            continue;
        };
        let Some(right_marker) = directional_marker_for_center(molecule, right_center, left_center)
        else {
            continue;
        };

        let Some(graph_bond_idx) = original_to_graph[bond_idx] else {
            continue;
        };

        specs.push(DoubleBondStereoSpec {
            bond_index: graph_bond_idx,
            atom1: left_marker.0,
            atom2: right_marker.0,
            is_z: left_marker.1 == right_marker.1,
        });
    }

    specs
}

fn directional_marker_for_center(
    molecule: &opensmiles::Molecule,
    center: usize,
    other_center: usize,
) -> Option<(usize, i8)> {
    molecule
        .bonds()
        .iter()
        .find_map(|bond| match bond.kind() {
            SmilesBondType::Up | SmilesBondType::Down => {
                let source = bond.source() as usize;
                let target = bond.target() as usize;
                let sign = match bond.kind() {
                    SmilesBondType::Up => 1,
                    SmilesBondType::Down => -1,
                    _ => unreachable!(),
                };

                if source == center && target != other_center {
                    Some((target, sign))
                } else if target == center && source != other_center {
                    Some((source, -sign))
                } else {
                    None
                }
            }
            _ => None,
        })
}

fn kekulize_aromatic_bonds(atom_count: usize, bonds: &mut [SmilesBond]) {
    let mut aromatic_adj = vec![Vec::new(); atom_count];
    for (bond_idx, bond) in bonds.iter().enumerate() {
        if bond.kind == SmilesBondKind::Aromatic {
            aromatic_adj[bond.atom1].push(bond_idx);
            aromatic_adj[bond.atom2].push(bond_idx);
        }
    }

    let mut visited = vec![false; atom_count];
    for start in 0..atom_count {
        if visited[start] || aromatic_adj[start].is_empty() {
            continue;
        }

        let mut queue = VecDeque::from([start]);
        let mut component_atoms = Vec::new();
        let mut component_edges = HashSet::new();
        visited[start] = true;

        while let Some(atom) = queue.pop_front() {
            component_atoms.push(atom);
            for &bond_idx in &aromatic_adj[atom] {
                component_edges.insert(bond_idx);
                let bond = &bonds[bond_idx];
                let other = if bond.atom1 == atom { bond.atom2 } else { bond.atom1 };
                if !visited[other] {
                    visited[other] = true;
                    queue.push_back(other);
                }
            }
        }

        let mut edge_indices = component_edges.into_iter().collect::<Vec<_>>();
        edge_indices.sort_unstable_by_key(|&bond_idx| {
            let bond = &bonds[bond_idx];
            usize::MAX - (aromatic_adj[bond.atom1].len() + aromatic_adj[bond.atom2].len())
        });

        let mut used_atoms = vec![false; atom_count];
        let mut current = Vec::new();
        let mut best = Vec::new();
        search_max_matching(
            &edge_indices,
            bonds,
            0,
            &mut used_atoms,
            &mut current,
            &mut best,
        );

        let best: HashSet<_> = best.into_iter().collect();
        for bond_idx in edge_indices {
            bonds[bond_idx].kind = if best.contains(&bond_idx) {
                SmilesBondKind::Double
            } else {
                SmilesBondKind::Single
            };
        }
    }
}

fn search_max_matching(
    edge_indices: &[usize],
    bonds: &[SmilesBond],
    pos: usize,
    used_atoms: &mut [bool],
    current: &mut Vec<usize>,
    best: &mut Vec<usize>,
) {
    if pos == edge_indices.len() {
        if current.len() > best.len() {
            *best = current.clone();
        }
        return;
    }

    let remaining = edge_indices.len() - pos;
    if current.len() + remaining <= best.len() {
        return;
    }

    let bond_idx = edge_indices[pos];
    let bond = &bonds[bond_idx];

    if !used_atoms[bond.atom1] && !used_atoms[bond.atom2] {
        used_atoms[bond.atom1] = true;
        used_atoms[bond.atom2] = true;
        current.push(bond_idx);
        search_max_matching(edge_indices, bonds, pos + 1, used_atoms, current, best);
        current.pop();
        used_atoms[bond.atom1] = false;
        used_atoms[bond.atom2] = false;
    }

    search_max_matching(edge_indices, bonds, pos + 1, used_atoms, current, best);
}

fn encode_layout_input(graph: &SmilesGraph) -> Vec<u8> {
    let atom_stereo_count = graph.atoms.iter().filter(|atom| atom.stereo.is_some()).count();
    let mut buffer = Vec::with_capacity(
        20
            + graph.atoms.len()
            + graph.bonds.len() * 9
            + atom_stereo_count * 17
            + graph.double_bond_stereo.len() * 13,
    );
    buffer.extend_from_slice(LAYOUT_MAGIC);
    write_u32_le(&mut buffer, graph.atoms.len());
    write_u32_le(&mut buffer, graph.bonds.len());
    write_u32_le(&mut buffer, atom_stereo_count);
    write_u32_le(&mut buffer, graph.double_bond_stereo.len());

    for atom in &graph.atoms {
        buffer.push(atom.atomic_number);
    }

    for bond in &graph.bonds {
        write_u32_le(&mut buffer, bond.atom1);
        write_u32_le(&mut buffer, bond.atom2);
        buffer.push(match bond.kind {
            SmilesBondKind::Single => 1,
            SmilesBondKind::Double => 2,
            SmilesBondKind::Triple => 3,
            SmilesBondKind::Aromatic => 1,
        });
    }

    for (atom_idx, atom) in graph.atoms.iter().enumerate() {
        let Some(stereo) = atom.stereo else {
            continue;
        };
        write_u32_le(&mut buffer, atom_idx);
        write_u32_le(
            &mut buffer,
            stereo.looking_from.unwrap_or(usize::MAX),
        );
        write_u32_le(&mut buffer, stereo.atom1.unwrap_or(usize::MAX));
        write_u32_le(&mut buffer, stereo.atom2.unwrap_or(usize::MAX));
        buffer.push(match stereo.direction {
            AtomStereoDirection::Clockwise => 1,
            AtomStereoDirection::CounterClockwise => 2,
        });
    }

    for stereo in &graph.double_bond_stereo {
        write_u32_le(&mut buffer, stereo.bond_index);
        write_u32_le(&mut buffer, stereo.atom1);
        write_u32_le(&mut buffer, stereo.atom2);
        buffer.push(u8::from(stereo.is_z));
    }

    buffer
}

fn decode_coords(
    coords_data: &[u8],
    expected_atom_count: usize,
    expected_bond_count: usize,
) -> Result<CoordPayload, String> {
    if coords_data.len() < 12 {
        return Err("Coordinate payload is too short".to_string());
    }
    if &coords_data[..4] != COORD_MAGIC {
        return Err("Invalid coordinate payload header".to_string());
    }

    let mut offset = 4;
    let atom_count = read_u32_le(coords_data, &mut offset)? as usize;
    let bond_count = read_u32_le(coords_data, &mut offset)? as usize;
    if atom_count != expected_atom_count {
        return Err(format!(
            "Coordinate payload has {} atoms, expected {}",
            atom_count, expected_atom_count
        ));
    }
    if bond_count != expected_bond_count {
        return Err(format!(
            "Coordinate payload has {} bonds, expected {}",
            bond_count, expected_bond_count
        ));
    }

    let expected_len = 12 + atom_count * 8 + bond_count;
    if coords_data.len() != expected_len {
        return Err("Coordinate payload length does not match atom count".to_string());
    }

    let mut coords = Vec::with_capacity(atom_count);
    for _ in 0..atom_count {
        let x = read_f32_le(coords_data, &mut offset)?;
        let y = read_f32_le(coords_data, &mut offset)?;
        coords.push((x, y));
    }

    let bond_styles = coords_data[offset..offset + bond_count].to_vec();

    Ok(CoordPayload { coords, bond_styles })
}

fn build_smiles_stereo_map(
    graph: &SmilesGraph,
    bond_styles: &[u8],
) -> HashMap<(usize, usize), (u8, bool)> {
    let mut stereo_map = HashMap::new();

    for (bond, style) in graph.bonds.iter().zip(bond_styles.iter().copied()) {
        match style {
            1 => insert_stereo_edge(&mut stereo_map, bond.atom1, bond.atom2, 1),
            2 => insert_stereo_edge(&mut stereo_map, bond.atom2, bond.atom1, 1),
            3 => insert_stereo_edge(&mut stereo_map, bond.atom1, bond.atom2, 6),
            4 => insert_stereo_edge(&mut stereo_map, bond.atom2, bond.atom1, 6),
            _ => {}
        }
    }

    stereo_map
}

fn insert_stereo_edge(
    stereo_map: &mut HashMap<(usize, usize), (u8, bool)>,
    from: usize,
    to: usize,
    stereo: u8,
) {
    stereo_map.insert((from, to), (stereo, true));
    stereo_map.insert((to, from), (stereo, false));
}

fn ast_from_render_mol(
    mol: &RenderMol,
    mode_str: &str,
    stereo_map: &HashMap<(usize, usize), (u8, bool)>,
) -> Result<Vec<u8>, String> {
    let mut adj: Vec<Vec<(usize, BondKind, usize)>> = vec![Vec::new(); mol.atoms.len()];
    let mut total_len = 0.0;
    let mut bond_count = 0usize;

    for (i, bond) in mol.bonds.iter().enumerate() {
        adj[bond.atom1].push((bond.atom2, bond.kind, i));
        adj[bond.atom2].push((bond.atom1, bond.kind, i));

        let u = &mol.atoms[bond.atom1];
        let v = &mol.atoms[bond.atom2];
        let dx = u.x - v.x;
        let dy = u.y - v.y;
        total_len += (dx * dx + dy * dy).sqrt();
        bond_count += 1;
    }

    let avg_length = if bond_count > 0 {
        total_len / bond_count as f64
    } else {
        1.0
    };
    let avg_length = if avg_length < 1e-6 { 1.0 } else { avg_length };

    let mut visited_nodes = vec![false; mol.atoms.len()];
    let mut handled_bonds = vec![false; mol.bonds.len()];
    let labels = build_labels(mol, &adj, mode_str, &mut visited_nodes, &mut handled_bonds);
    let rings = find_rings(&adj, mol.atoms.len());

    let mut root_commands = Vec::new();
    for start_node in 0..mol.atoms.len() {
        if !visited_nodes[start_node] {
            dfs(
                start_node,
                &adj,
                mol,
                &labels,
                stereo_map,
                &rings,
                avg_length,
                &mut visited_nodes,
                &mut handled_bonds,
                &mut root_commands,
            );
        }
    }

    let mut buffer = Vec::new();
    ciborium::into_writer(&root_commands, &mut buffer).map_err(|e| e.to_string())?;
    Ok(buffer)
}

fn build_labels(
    mol: &RenderMol,
    adj: &[Vec<(usize, BondKind, usize)>],
    mode_str: &str,
    visited_nodes: &mut [bool],
    handled_bonds: &mut [bool],
) -> Vec<String> {
    let mut labels = vec![String::new(); mol.atoms.len()];

    for i in 0..mol.atoms.len() {
        let atom = &mol.atoms[i];

        if mode_str == "abbreviate" || mode_str == "skeletal" {
            if atom.element == "H" {
                visited_nodes[i] = true;
                continue;
            }

            let (explicit_h_count, explicit_h_bonds) = explicit_h_neighbors(i, mol, adj);
            for bond_idx in explicit_h_bonds {
                handled_bonds[bond_idx] = true;
            }

            let total_h = atom.hydrogens + explicit_h_count;
            if atom.element == "C" {
                if mode_str == "skeletal" {
                    labels[i] = if atom.stereo_annotation.is_some() {
                        atom.element.clone()
                    } else {
                        String::new()
                    };
                } else {
                    let heavy_atoms = adj[i].len().saturating_sub(explicit_h_count as usize);
                    labels[i] = if heavy_atoms <= 1 && total_h > 0 {
                        atom_label(&atom.element, total_h, atom.charge)
                    } else if atom.stereo_annotation.is_some() || atom.charge != 0 {
                        atom_label(&atom.element, 0, atom.charge)
                    } else {
                        String::new()
                    };
                }
            } else {
                labels[i] = atom_label(&atom.element, total_h, atom.charge);
            }
        } else {
            labels[i] = atom_label(&atom.element, 0, atom.charge);
        }
    }

    labels
}

fn explicit_h_neighbors(
    atom_idx: usize,
    mol: &RenderMol,
    adj: &[Vec<(usize, BondKind, usize)>],
) -> (u8, Vec<usize>) {
    let mut count = 0u8;
    let mut bond_indices = Vec::new();
    for &(neighbor, _, bond_idx) in &adj[atom_idx] {
        if mol.atoms[neighbor].element == "H" {
            count = count.saturating_add(1);
            bond_indices.push(bond_idx);
        }
    }
    (count, bond_indices)
}

fn atom_label(element: &str, hydrogen_count: u8, charge: i8) -> String {
    let mut label = if hydrogen_count == 0 {
        element.to_string()
    } else if hydrogen_count == 1 {
        format!("{element}H")
    } else {
        format!("{element}H_{hydrogen_count}")
    };
    label.push_str(&charge_suffix(charge));
    label
}

fn charge_suffix(charge: i8) -> String {
    match charge {
        0 => String::new(),
        1 => "^+".to_string(),
        -1 => "^-".to_string(),
        charge if charge > 1 => format!("^{}+", charge),
        charge => format!("^{}-", charge.abs()),
    }
}

fn extract_stereo(sdf_str: &str) -> HashMap<(usize, usize), (u8, bool)> {
    let mut map = HashMap::new();
    let mut lines = sdf_str.lines();

    lines.next();
    lines.next();
    lines.next();

    if let Some(counts_line) = lines.next() {
        if counts_line.len() >= 6 {
            let num_atoms = counts_line[0..3].trim().parse::<usize>().unwrap_or(0);
            let num_bonds = counts_line[3..6].trim().parse::<usize>().unwrap_or(0);

            for _ in 0..num_atoms {
                lines.next();
            }

            for _ in 0..num_bonds {
                if let Some(bond_line) = lines.next() {
                    if bond_line.len() >= 12 {
                        let a1 = bond_line[0..3].trim().parse::<usize>().unwrap_or(0);
                        let a2 = bond_line[3..6].trim().parse::<usize>().unwrap_or(0);
                        let stereo = bond_line[9..12].trim().parse::<u8>().unwrap_or(0);

                        if a1 > 0 && a2 > 0 && (stereo == 1 || stereo == 6) {
                            map.insert((a1 - 1, a2 - 1), (stereo, true));
                            map.insert((a2 - 1, a1 - 1), (stereo, false));
                        }
                    }
                }
            }
        }
    }
    map
}

fn find_rings(adj: &[Vec<(usize, BondKind, usize)>], num_atoms: usize) -> Vec<Vec<usize>> {
    let mut rings = Vec::new();
    let mut seen = HashSet::new();

    for u in 0..num_atoms {
        for &(v, _, bond_idx) in &adj[u] {
            if u >= v {
                continue;
            }

            let mut queue = VecDeque::new();
            let mut parent = vec![usize::MAX; num_atoms];
            queue.push_back(u);
            parent[u] = u;

            let mut found = false;
            while let Some(curr) = queue.pop_front() {
                if curr == v {
                    found = true;
                    break;
                }
                for &(next, _, b_idx) in &adj[curr] {
                    if b_idx == bond_idx {
                        continue;
                    }
                    if parent[next] == usize::MAX {
                        parent[next] = curr;
                        queue.push_back(next);
                    }
                }
            }

            if found {
                let mut cycle = Vec::new();
                let mut curr = v;
                while curr != u {
                    cycle.push(curr);
                    curr = parent[curr];
                }
                cycle.push(u);

                if (3..=8).contains(&cycle.len()) {
                    let min_idx = cycle
                        .iter()
                        .enumerate()
                        .min_by_key(|&(_, &val)| val)
                        .map(|(idx, _)| idx)
                        .unwrap_or(0);
                    let n = cycle.len();
                    let prev = cycle[(min_idx + n - 1) % n];
                    let next = cycle[(min_idx + 1) % n];
                    let mut normalized = Vec::with_capacity(n);

                    if next < prev {
                        for i in 0..n {
                            normalized.push(cycle[(min_idx + i) % n]);
                        }
                    } else {
                        for i in 0..n {
                            normalized.push(cycle[(min_idx + n - i) % n]);
                        }
                    }

                    if seen.insert(normalized.clone()) {
                        rings.push(normalized);
                    }
                }
            }
        }
    }
    rings
}

fn calculate_double_bond_offset(
    u: usize,
    v: usize,
    mol: &RenderMol,
    rings: &[Vec<usize>],
) -> Option<String> {
    let ring = rings.iter().find(|ring| ring.contains(&u) && ring.contains(&v))?;

    let mut cx = 0.0;
    let mut cy = 0.0;
    for &node in ring {
        cx += mol.atoms[node].x;
        cy += mol.atoms[node].y;
    }
    cx /= ring.len() as f64;
    cy /= ring.len() as f64;

    let a = &mol.atoms[u];
    let b = &mol.atoms[v];
    let dx = b.x - a.x;
    let dy = b.y - a.y;
    let cx_dir = cx - a.x;
    let cy_dir = cy - a.y;
    let cross = dx * cy_dir - dy * cx_dir;

    Some(if cross > 0.0 { "left" } else { "right" }.to_string())
}

fn calc_angle(u: &RenderAtom, v: &RenderAtom) -> f64 {
    let dx = v.x - u.x;
    let dy = v.y - u.y;
    dy.atan2(dx).to_degrees()
}

fn calc_length_scale(u: &RenderAtom, v: &RenderAtom, avg_len: f64) -> f64 {
    let dx = u.x - v.x;
    let dy = u.y - v.y;
    let len = (dx * dx + dy * dy).sqrt();
    len / avg_len
}

fn bond_func_name(
    kind: BondKind,
    u: usize,
    v: usize,
    stereo_map: &HashMap<(usize, usize), (u8, bool)>,
) -> &'static str {
    if kind == BondKind::Single {
        if let Some(&(stereo, is_forward)) = stereo_map.get(&(u, v)) {
            match stereo {
                1 => {
                    return if is_forward {
                        "cram-filled-left"
                    } else {
                        "cram-filled-right"
                    };
                }
                6 => {
                    return if is_forward {
                        "cram-dashed-left"
                    } else {
                        "cram-dashed-right"
                    };
                }
                _ => {}
            }
        }
    }

    match kind {
        BondKind::Double => "double",
        BondKind::Triple => "triple",
        BondKind::Single => "single",
    }
}

fn dfs(
    u: usize,
    adj: &[Vec<(usize, BondKind, usize)>],
    mol: &RenderMol,
    labels: &[String],
    stereo_map: &HashMap<(usize, usize), (u8, bool)>,
    rings: &[Vec<usize>],
    avg_length: f64,
    visited_nodes: &mut [bool],
    handled_bonds: &mut [bool],
    commands: &mut Vec<Command>,
) {
    visited_nodes[u] = true;
    let u_atom = &mol.atoms[u];

    let mut back_edges = Vec::new();
    for &(v, kind, bond_idx) in &adj[u] {
        if handled_bonds[bond_idx] {
            continue;
        }
        if visited_nodes[v] {
            back_edges.push((v, kind, bond_idx));
        }
    }

    let mut links = Vec::new();
    for &(v, kind, bond_idx) in &back_edges {
        handled_bonds[bond_idx] = true;
        let offset = if kind == BondKind::Double {
            calculate_double_bond_offset(u, v, mol, rings)
        } else {
            None
        };
        links.push(LinkData {
            target: format!("a{v}"),
            bond_type: bond_func_name(kind, u, v, stereo_map).to_string(),
            angle: calc_angle(u_atom, &mol.atoms[v]),
            offset,
            length_scale: calc_length_scale(u_atom, &mol.atoms[v], avg_length),
        });
    }

    commands.push(Command::Fragment {
        element: labels[u].clone(),
        name: format!("a{u}"),
        links,
        annotation: u_atom.stereo_annotation.clone(),
    });

    let mut forward_targets = Vec::new();
    for &(v, kind, bond_idx) in &adj[u] {
        if !handled_bonds[bond_idx] && !visited_nodes[v] {
            forward_targets.push((v, kind, bond_idx));
        }
    }

    for (idx, &(v, kind, bond_idx)) in forward_targets.iter().enumerate() {
        if visited_nodes[v] {
            continue;
        }
        handled_bonds[bond_idx] = true;

        let has_more = forward_targets[idx + 1..]
            .iter()
            .any(|(next_v, _, _)| !visited_nodes[*next_v]);

        let offset = if kind == BondKind::Double {
            calculate_double_bond_offset(u, v, mol, rings)
        } else {
            None
        };

        let bond_cmd = Command::Bond {
            bond_type: bond_func_name(kind, u, v, stereo_map).to_string(),
            angle: calc_angle(u_atom, &mol.atoms[v]),
            offset,
            length_scale: calc_length_scale(u_atom, &mol.atoms[v], avg_length),
        };

        if has_more {
            let mut branch_body = vec![bond_cmd];
            dfs(
                v,
                adj,
                mol,
                labels,
                stereo_map,
                rings,
                avg_length,
                visited_nodes,
                handled_bonds,
                &mut branch_body,
            );
            commands.push(Command::Branch { body: branch_body });
        } else {
            commands.push(bond_cmd);
            dfs(
                v,
                adj,
                mol,
                labels,
                stereo_map,
                rings,
                avg_length,
                visited_nodes,
                handled_bonds,
                commands,
            );
        }
    }
}

fn write_u32_le(buffer: &mut Vec<u8>, value: usize) {
    buffer.extend_from_slice(&(value as u32).to_le_bytes());
}

fn read_u32_le(bytes: &[u8], offset: &mut usize) -> Result<u32, String> {
    let end = *offset + 4;
    let chunk = bytes
        .get(*offset..end)
        .ok_or_else(|| "Unexpected end of payload".to_string())?;
    *offset = end;
    Ok(u32::from_le_bytes(chunk.try_into().unwrap()))
}

fn read_f32_le(bytes: &[u8], offset: &mut usize) -> Result<f32, String> {
    let end = *offset + 4;
    let chunk = bytes
        .get(*offset..end)
        .ok_or_else(|| "Unexpected end of payload".to_string())?;
    *offset = end;
    Ok(f32::from_le_bytes(chunk.try_into().unwrap()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use ciborium::Value;

    fn decode_layout_input(
        bytes: &[u8],
    ) -> (
        Vec<u8>,
        Vec<(u32, u32, u8)>,
        Vec<(u32, u32, u32, u32, u8)>,
        Vec<(u32, u32, u32, u8)>,
    ) {
        assert_eq!(&bytes[..4], LAYOUT_MAGIC);
        let mut offset = 4;
        let atom_count = read_u32_le(bytes, &mut offset).unwrap() as usize;
        let bond_count = read_u32_le(bytes, &mut offset).unwrap() as usize;
        let atom_stereo_count = read_u32_le(bytes, &mut offset).unwrap() as usize;
        let double_bond_stereo_count = read_u32_le(bytes, &mut offset).unwrap() as usize;

        let atoms = bytes[offset..offset + atom_count].to_vec();
        offset += atom_count;

        let mut bonds = Vec::new();
        for _ in 0..bond_count {
            let atom1 = read_u32_le(bytes, &mut offset).unwrap();
            let atom2 = read_u32_le(bytes, &mut offset).unwrap();
            let order = bytes[offset];
            offset += 1;
            bonds.push((atom1, atom2, order));
        }

        let mut atom_stereo = Vec::new();
        for _ in 0..atom_stereo_count {
            let atom = read_u32_le(bytes, &mut offset).unwrap();
            let looking_from = read_u32_le(bytes, &mut offset).unwrap();
            let atom1 = read_u32_le(bytes, &mut offset).unwrap();
            let atom2 = read_u32_le(bytes, &mut offset).unwrap();
            let direction = bytes[offset];
            offset += 1;
            atom_stereo.push((atom, looking_from, atom1, atom2, direction));
        }

        let mut double_bond_stereo = Vec::new();
        for _ in 0..double_bond_stereo_count {
            let bond = read_u32_le(bytes, &mut offset).unwrap();
            let atom1 = read_u32_le(bytes, &mut offset).unwrap();
            let atom2 = read_u32_le(bytes, &mut offset).unwrap();
            let is_z = bytes[offset];
            offset += 1;
            double_bond_stereo.push((bond, atom1, atom2, is_z));
        }

        (atoms, bonds, atom_stereo, double_bond_stereo)
    }

    fn collect_fragment_elements(commands: &[Value], output: &mut Vec<String>) {
        for command in commands {
            let Some(map) = command.as_map() else {
                continue;
            };
            let command_type = map
                .iter()
                .find(|(key, _)| key.as_text() == Some("type"))
                .and_then(|(_, value)| value.as_text());
            match command_type {
                Some("fragment") => {
                    if let Some(label) = map
                        .iter()
                        .find(|(key, _)| key.as_text() == Some("element"))
                        .and_then(|(_, value)| value.as_text())
                    {
                        if !label.is_empty() {
                            output.push(label.to_string());
                        }
                    }
                }
                Some("branch") => {
                    if let Some(body) = map
                        .iter()
                        .find(|(key, _)| key.as_text() == Some("body"))
                        .and_then(|(_, value)| value.as_array())
                    {
                        collect_fragment_elements(body, output);
                    }
                }
                _ => {}
            }
        }
    }

    #[test]
    fn smiles_payload_keeps_basic_connectivity() {
        let payload = smiles_to_layout_input(b"CCO").unwrap();
        let (atoms, bonds, atom_stereo, double_bond_stereo) = decode_layout_input(&payload);

        assert_eq!(atoms, vec![6, 6, 8]);
        assert_eq!(bonds, vec![(0, 1, 1), (1, 2, 1)]);
        assert!(atom_stereo.is_empty());
        assert!(double_bond_stereo.is_empty());
    }

    #[test]
    fn smiles_full_layout_input_expands_implicit_hydrogens() {
        let payload = smiles_to_full_layout_input(b"CCO").unwrap();
        let (atoms, bonds, atom_stereo, double_bond_stereo) = decode_layout_input(&payload);

        assert_eq!(atoms, vec![6, 6, 8, 1, 1, 1, 1, 1, 1]);
        assert_eq!(bonds.len(), 8);
        assert!(atom_stereo.is_empty());
        assert!(double_bond_stereo.is_empty());
    }

    #[test]
    fn aromatic_ring_is_kekulized_for_layout() {
        let payload = smiles_to_layout_input(b"c1ccccc1").unwrap();
        let (_, bonds, _, _) = decode_layout_input(&payload);
        let double_count = bonds.iter().filter(|(_, _, order)| *order == 2).count();
        let single_count = bonds.iter().filter(|(_, _, order)| *order == 1).count();

        assert_eq!(bonds.len(), 6);
        assert_eq!(double_count, 3);
        assert_eq!(single_count, 3);
    }

    #[test]
    fn smiles_abbreviate_mode_uses_folded_hydrogen_labels() {
        let mut coords = Vec::new();
        coords.extend_from_slice(COORD_MAGIC);
        coords.extend_from_slice(&(3u32).to_le_bytes());
        coords.extend_from_slice(&(2u32).to_le_bytes());
        for (x, y) in [(0.0f32, 0.0f32), (1.0, 0.0), (2.0, 0.0)] {
            coords.extend_from_slice(&x.to_le_bytes());
            coords.extend_from_slice(&y.to_le_bytes());
        }
        coords.extend_from_slice(&[0, 0]);

        let ast_bytes = smiles_to_ast(b"CCO", &coords, b"abbreviate").unwrap();
        let ast: Value = ciborium::from_reader(ast_bytes.as_slice()).unwrap();
        let ast = ast.as_array().unwrap();
        let mut labels = Vec::new();
        collect_fragment_elements(ast, &mut labels);

        assert_eq!(labels, vec!["CH_3".to_string(), "OH".to_string()]);
    }

    #[test]
    fn smiles_abbreviate_mode_preserves_formal_charge_labels() {
        let mut coords = Vec::new();
        coords.extend_from_slice(COORD_MAGIC);
        coords.extend_from_slice(&(2u32).to_le_bytes());
        coords.extend_from_slice(&(1u32).to_le_bytes());
        for (x, y) in [(0.0f32, 0.0f32), (1.0, 0.0)] {
            coords.extend_from_slice(&x.to_le_bytes());
            coords.extend_from_slice(&y.to_le_bytes());
        }
        coords.extend_from_slice(&[0]);

        let ast_bytes = smiles_to_ast(b"[n+]C", &coords, b"abbreviate").unwrap();
        let ast: Value = ciborium::from_reader(ast_bytes.as_slice()).unwrap();
        let ast = ast.as_array().unwrap();
        let mut labels = Vec::new();
        collect_fragment_elements(ast, &mut labels);

        assert_eq!(labels, vec!["N^+".to_string(), "CH_3".to_string()]);
    }

    #[test]
    fn smiles_full_mode_expands_implicit_hydrogens() {
        let mut coords = Vec::new();
        coords.extend_from_slice(COORD_MAGIC);
        coords.extend_from_slice(&(9u32).to_le_bytes());
        coords.extend_from_slice(&(8u32).to_le_bytes());
        for (x, y) in [
            (0.0f32, 0.0f32),
            (1.0, 0.0),
            (2.0, 0.0),
            (-0.5, 0.9),
            (-0.5, -0.9),
            (0.0, 1.1),
            (1.0, 1.0),
            (1.0, -1.0),
            (2.5, 0.8),
        ] {
            coords.extend_from_slice(&x.to_le_bytes());
            coords.extend_from_slice(&y.to_le_bytes());
        }
        coords.extend_from_slice(&[0; 8]);

        let ast_bytes = smiles_to_ast(b"CCO", &coords, b"full").unwrap();
        let ast: Value = ciborium::from_reader(ast_bytes.as_slice()).unwrap();
        let ast = ast.as_array().unwrap();
        let mut labels = Vec::new();
        collect_fragment_elements(ast, &mut labels);
        labels.sort();

        assert_eq!(
            labels,
            vec![
                "C".to_string(),
                "C".to_string(),
                "H".to_string(),
                "H".to_string(),
                "H".to_string(),
                "H".to_string(),
                "H".to_string(),
                "H".to_string(),
                "O".to_string(),
            ]
        );
    }

    #[test]
    fn tetrahedral_chirality_is_encoded_for_layout() {
        let payload = smiles_to_layout_input(b"N[C@@H](C)C(=O)O").unwrap();
        let (_, _, atom_stereo, _) = decode_layout_input(&payload);

        assert_eq!(atom_stereo.len(), 1);
        let (atom, looking_from, atom1, atom2, direction) = atom_stereo[0];
        assert_eq!(atom, 1);
        assert_eq!(looking_from, 0);
        assert_eq!(atom1, 2);
        assert_eq!(atom2, 3);
        assert_eq!(direction, 1);
    }

    #[test]
    fn ez_stereo_is_encoded_for_layout() {
        let trans = smiles_to_layout_input(br"F/C=C/F").unwrap();
        let cis = smiles_to_layout_input(br"F/C=C\F").unwrap();

        let (_, _, _, trans_db) = decode_layout_input(&trans);
        let (_, _, _, cis_db) = decode_layout_input(&cis);

        assert_eq!(trans_db, vec![(1, 0, 3, 0)]);
        assert_eq!(cis_db, vec![(1, 0, 3, 1)]);
    }

    #[test]
    fn complex_reported_smiles_is_accepted() {
        let smiles = b"CC[C@@H]([C@@H]1[C@H](C[C@@](O1)(CC)[C@H]2CC[C@@]([C@@H](O2)C)(CC)O)C)C(=O)[C@@H](C)[C@H]([C@H](C)CCC3=C(C=C(C(=C3C(=O)O)O)C)Br)O";
        let payload = smiles_to_layout_input(smiles).unwrap();
        let (atoms, bonds, atom_stereo, _) = decode_layout_input(&payload);

        assert!(!atoms.is_empty());
        assert!(!bonds.is_empty());
        assert!(!atom_stereo.is_empty());
    }

    #[test]
    fn non_tetrahedral_chirality_is_preserved_in_ast_annotations() {
        let mut coords = Vec::new();
        coords.extend_from_slice(COORD_MAGIC);
        coords.extend_from_slice(&(5u32).to_le_bytes());
        coords.extend_from_slice(&(4u32).to_le_bytes());
        for (x, y) in [(0.0f32, 0.0f32), (1.0, 0.0), (-1.0, 0.0), (0.0, 1.0), (0.0, -1.0)] {
            coords.extend_from_slice(&x.to_le_bytes());
            coords.extend_from_slice(&y.to_le_bytes());
        }
        coords.extend_from_slice(&[0, 0, 0, 0]);

        let ast_bytes = smiles_to_ast(b"[Pt@SP1](F)(Cl)(Br)I", &coords, b"skeletal").unwrap();
        let ast: ciborium::Value = ciborium::from_reader(ast_bytes.as_slice()).unwrap();
        let ast = ast.as_array().unwrap();
        let first = ast[0].as_map().unwrap();
        let first_label = first
            .iter()
            .find(|(key, _)| key.as_text() == Some("element"))
            .and_then(|(_, value)| value.as_text())
            .unwrap();
        let first_annotation = first
            .iter()
            .find(|(key, _)| key.as_text() == Some("annotation"))
            .and_then(|(_, value)| value.as_text())
            .unwrap();

        assert_eq!(first_label, "Pt");
        assert_eq!(first_annotation, "@SP1");
    }

    #[test]
    fn all_supported_chirality_classes_are_parsed() {
        let al = parse_smiles_graph("[C@AL1](F)(Cl)C").unwrap();
        let sp = parse_smiles_graph("[Pt@SP2](F)(Cl)(Br)I").unwrap();
        let tb = parse_smiles_graph("[As@TB5](F)(Cl)(Br)(N)S").unwrap();
        let oh = parse_smiles_graph("[Co@OH5](F)(Cl)(Br)(I)(N)S").unwrap();

        assert_eq!(al.atoms[0].stereo_annotation.as_deref(), Some("@AL1"));
        assert_eq!(sp.atoms[0].stereo_annotation.as_deref(), Some("@SP2"));
        assert_eq!(tb.atoms[0].stereo_annotation.as_deref(), Some("@TB5"));
        assert_eq!(oh.atoms[0].stereo_annotation.as_deref(), Some("@OH5"));
    }
}
