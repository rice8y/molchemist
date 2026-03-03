use serde::Serialize;
use wasm_minimal_protocol::*;
use sdfrust::{parse_sdf_string, Atom, BondOrder, Molecule};
use std::collections::{HashMap, HashSet, VecDeque};

initiate_protocol!();

#[derive(Serialize)]
#[serde(tag = "type")]
enum Command {
    #[serde(rename = "fragment")]
    Fragment {
        element: String,
        name: String,
        links: Vec<LinkData>,
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
    Branch {
        body: Vec<Command>,
    },
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

#[wasm_func]
pub fn sdf_to_ast(sdf_data: &[u8], options: &[u8]) -> Result<Vec<u8>, String> {
    let sdf_str = std::str::from_utf8(sdf_data).map_err(|e| e.to_string())?;
    let mol = parse_sdf_string(sdf_str).map_err(|e| e.to_string())?;
    let mode_str = std::str::from_utf8(options).unwrap_or("full");

    let stereo_map = extract_stereo(sdf_str);

    let mut adj: Vec<Vec<(usize, &BondOrder, usize)>> = vec![Vec::new(); mol.atom_count()];
    let mut total_len = 0.0;
    let mut bond_count = 0;

    for (i, bond) in mol.bonds.iter().enumerate() {
        adj[bond.atom1].push((bond.atom2, &bond.order, i));
        adj[bond.atom2].push((bond.atom1, &bond.order, i));
        
        let u = &mol.atoms[bond.atom1];
        let v = &mol.atoms[bond.atom2];
        let dx = u.x - v.x;
        let dy = u.y - v.y;
        total_len += (dx * dx + dy * dy).sqrt();
        bond_count += 1;
    }
    
    let avg_length = if bond_count > 0 { total_len / bond_count as f64 } else { 1.0 };
    let avg_length = if avg_length < 1e-6 { 1.0 } else { avg_length };

    let mut visited_nodes = vec![false; mol.atom_count()];
    let mut handled_bonds = vec![false; mol.bonds.len()];
    let mut labels = vec!["".to_string(); mol.atom_count()];

    for i in 0..mol.atom_count() {
        if mode_str == "abbreviate" || mode_str == "skeletal" {
            if mol.atoms[i].element == "H" {
                visited_nodes[i] = true;
            } else {
                let mut h_count = 0;
                for &(v, _, bond_idx) in &adj[i] {
                    if mol.atoms[v].element == "H" {
                        h_count += 1;
                        handled_bonds[bond_idx] = true;
                    }
                }
                if mol.atoms[i].element == "C" {
                    if mode_str == "skeletal" {
                        labels[i] = "".to_string();
                    } else {
                        let heavy_atoms = adj[i].len() - h_count;
                        if heavy_atoms <= 1 && h_count > 0 {
                            labels[i] = if h_count == 1 { "CH".to_string() } else { format!("CH_{}", h_count) };
                        } else {
                            labels[i] = "".to_string();
                        }
                    }
                } else {
                    if h_count == 0 {
                        labels[i] = mol.atoms[i].element.clone();
                    } else if h_count == 1 {
                        labels[i] = format!("{}H", mol.atoms[i].element);
                    } else {
                        labels[i] = format!("{}H_{}", mol.atoms[i].element, h_count);
                    }
                }
            }
        } else {
            labels[i] = mol.atoms[i].element.clone();
        }
    }

    let rings = find_rings(&adj, mol.atom_count());

    let mut root_commands = Vec::new();
    for start_node in 0..mol.atom_count() {
        if !visited_nodes[start_node] {
            dfs(
                start_node,
                &adj,
                &mol,
                &labels,
                &stereo_map,
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

fn extract_stereo(sdf_str: &str) -> HashMap<(usize, usize), (u8, bool)> {
    let mut map = HashMap::new();
    let mut lines = sdf_str.lines();
    
    lines.next(); lines.next(); lines.next();
    
    if let Some(counts_line) = lines.next() {
        if counts_line.len() >= 6 {
            let num_atoms = counts_line[0..3].trim().parse::<usize>().unwrap_or(0);
            let num_bonds = counts_line[3..6].trim().parse::<usize>().unwrap_or(0);
            
            for _ in 0..num_atoms { lines.next(); }
            
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

fn find_rings(adj: &[Vec<(usize, &BondOrder, usize)>], num_atoms: usize) -> Vec<Vec<usize>> {
    let mut rings = Vec::new();
    let mut seen = HashSet::new();

    for u in 0..num_atoms {
        for &(v, _, bond_idx) in &adj[u] {
            if u >= v { continue; } 
            
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
                    if b_idx == bond_idx { continue; } 
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

                if cycle.len() >= 3 && cycle.len() <= 8 {
                    let min_idx = cycle.iter().enumerate().min_by_key(|&(_, &val)| val).unwrap().0;
                    let mut normalized = Vec::new();
                    let n = cycle.len();
                    
                    let prev = cycle[(min_idx + n - 1) % n];
                    let next = cycle[(min_idx + 1) % n];
                    
                    if next < prev {
                        for i in 0..n { normalized.push(cycle[(min_idx + i) % n]); }
                    } else {
                        for i in 0..n { normalized.push(cycle[(min_idx + n - i) % n]); }
                    }
                    
                    if !seen.contains(&normalized) {
                        seen.insert(normalized.clone());
                        rings.push(normalized);
                    }
                }
            }
        }
    }
    rings
}

fn calculate_double_bond_offset(u: usize, v: usize, mol: &Molecule, rings: &[Vec<usize>]) -> Option<String> {
    let ring = rings.iter().find(|r| r.contains(&u) && r.contains(&v));
    if let Some(r) = ring {
        let mut cx = 0.0;
        let mut cy = 0.0;
        for &node in r {
            cx += mol.atoms[node].x;
            cy += mol.atoms[node].y;
        }
        cx /= r.len() as f64;
        cy /= r.len() as f64;

        let ax = mol.atoms[u].x;
        let ay = mol.atoms[u].y;
        let bx = mol.atoms[v].x;
        let by = mol.atoms[v].y;

        let dx = bx - ax;
        let dy = by - ay;
        let cx_dir = cx - ax;
        let cy_dir = cy - ay;

        let cross = dx * cy_dir - dy * cx_dir;
        
        if cross > 0.0 {
            return Some("left".to_string());
        } else {
            return Some("right".to_string());
        }
    }
    None
}

fn calc_angle(u: &Atom, v: &Atom) -> f64 {
    let dx = v.x - u.x;
    let dy = v.y - u.y;
    dy.atan2(dx).to_degrees()
}

fn calc_length_scale(u: &Atom, v: &Atom, avg_len: f64) -> f64 {
    let dx = u.x - v.x;
    let dy = u.y - v.y;
    let len = (dx * dx + dy * dy).sqrt();
    len / avg_len
}

fn bond_func_name(
    order: &BondOrder, 
    u: usize, 
    v: usize, 
    stereo_map: &HashMap<(usize, usize), (u8, bool)>
) -> &'static str {
    if order == &BondOrder::Single {
        if let Some(&(stereo, is_forward)) = stereo_map.get(&(u, v)) {
            match stereo {
                1 => return if is_forward { "cram-filled-left" } else { "cram-filled-right" },
                6 => return if is_forward { "cram-dashed-left" } else { "cram-dashed-right" },
                _ => {}
            }
        }
    }
    match order {
        BondOrder::Double => "double",
        BondOrder::Triple => "triple",
        _ => "single",
    }
}

fn dfs(
    u: usize,
    adj: &[Vec<(usize, &BondOrder, usize)>],
    mol: &Molecule,
    labels: &[String],
    stereo_map: &HashMap<(usize, usize), (u8, bool)>,
    rings: &[Vec<usize>],
    avg_length: f64,
    visited_nodes: &mut Vec<bool>,
    handled_bonds: &mut Vec<bool>,
    commands: &mut Vec<Command>,
) {
    visited_nodes[u] = true;
    let u_atom = &mol.atoms[u];

    let mut back_edges = Vec::new();
    for &(v, order, bond_idx) in &adj[u] {
        if handled_bonds[bond_idx] { continue; }
        if visited_nodes[v] {
            back_edges.push((v, order, bond_idx));
        }
    }

    let mut links = Vec::new();
    for &(v, order, bond_idx) in &back_edges {
        handled_bonds[bond_idx] = true;
        let offset = if order == &BondOrder::Double {
            calculate_double_bond_offset(u, v, mol, rings)
        } else {
            None
        };
        links.push(LinkData {
            target: format!("a{}", v),
            bond_type: bond_func_name(order, u, v, stereo_map).to_string(),
            angle: calc_angle(u_atom, &mol.atoms[v]),
            offset,
            length_scale: calc_length_scale(u_atom, &mol.atoms[v], avg_length),
        });
    }

    commands.push(Command::Fragment {
        element: labels[u].clone(),
        name: format!("a{}", u),
        links,
    });

    let mut forward_targets = Vec::new();
    for &(v, order, bond_idx) in &adj[u] {
        if !handled_bonds[bond_idx] && !visited_nodes[v] {
            forward_targets.push((v, order, bond_idx));
        }
    }

    for (idx, &(v, order, bond_idx)) in forward_targets.iter().enumerate() {
        if visited_nodes[v] { continue; }
        handled_bonds[bond_idx] = true;

        let mut has_more = false;
        for &(next_v, _, _) in &forward_targets[idx + 1..] {
            if !visited_nodes[next_v] {
                has_more = true;
                break;
            }
        }

        let offset = if order == &BondOrder::Double {
            calculate_double_bond_offset(u, v, mol, rings)
        } else {
            None
        };

        let bond_cmd = Command::Bond {
            bond_type: bond_func_name(order, u, v, stereo_map).to_string(),
            angle: calc_angle(u_atom, &mol.atoms[v]),
            offset,
            length_scale: calc_length_scale(u_atom, &mol.atoms[v], avg_length),
        };

        if has_more {
            let mut branch_body = Vec::new();
            branch_body.push(bond_cmd);
            dfs(v, adj, mol, labels, stereo_map, rings, avg_length, visited_nodes, handled_bonds, &mut branch_body);
            commands.push(Command::Branch { body: branch_body });
        } else {
            commands.push(bond_cmd);
            dfs(v, adj, mol, labels, stereo_map, rings, avg_length, visited_nodes, handled_bonds, commands);
        }
    }
}