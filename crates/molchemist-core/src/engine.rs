use crate::{
    parse as parse_smiles, AtomSymbol, BondType as SmilesBondType, Chirality as SmilesChirality,
    Molecule as SmilesMolecule,
};
use sdfrust::{parse_sdf_auto_string, BondOrder, BondStereo, SdfFormat};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet, VecDeque};
const LAYOUT_MAGIC: &[u8; 4] = b"MCG2";
const COORD_MAGIC: &[u8; 4] = b"MCC2";

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(tag = "type")]
pub enum Command {
    #[serde(rename = "fragment")]
    Fragment {
        element: String,
        name: String,
        links: Vec<LinkData>,
        #[serde(skip_serializing_if = "Option::is_none")]
        atom: Option<AtomLabel>,
        #[serde(skip_serializing_if = "Option::is_none")]
        annotation: Option<String>,
    },
    #[serde(rename = "bond")]
    Bond {
        name: String,
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

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct LinkData {
    pub target: String,
    pub name: String,
    #[serde(rename = "bondType")]
    pub bond_type: String,
    pub angle: f64,
    pub offset: Option<String>,
    #[serde(rename = "lengthScale")]
    pub length_scale: f64,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct AtomLabel {
    pub symbol: String,
    #[serde(rename = "hydrogenCount")]
    pub hydrogen_count: u8,
    pub charge: i8,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub isotope: Option<u16>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub radical: Option<u8>,
    #[serde(rename = "atomMap", skip_serializing_if = "Option::is_none")]
    pub atom_map: Option<u32>,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum RenderMode {
    #[default]
    Full,
    Abbreviate,
    Skeletal,
}

impl RenderMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Full => "full",
            Self::Abbreviate => "abbreviate",
            Self::Skeletal => "skeletal",
        }
    }

    pub fn parse(value: &str) -> Result<Self, String> {
        match value {
            "full" => Ok(Self::Full),
            "abbreviate" => Ok(Self::Abbreviate),
            "skeletal" => Ok(Self::Skeletal),
            _ => Err(format!(
                "Invalid rendering mode {value:?}; expected full, abbreviate, or skeletal"
            )),
        }
    }
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
    isotope: Option<u16>,
    radical: Option<u8>,
    atom_map: Option<u32>,
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

#[derive(Clone, Default)]
struct SdfAtomMetadata {
    isotope: Option<u16>,
    radical: Option<u8>,
    atom_map: Option<u32>,
}

#[derive(Clone, Default)]
struct RenderLabel {
    text: String,
    atom: Option<AtomLabel>,
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
    isotope: Option<u16>,
    atom_map: Option<u32>,
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

pub fn sdf_to_ast(sdf_data: &[u8], options: &[u8]) -> Result<Vec<u8>, String> {
    sdf_record_to_ast(sdf_data, options, 1)
}

pub fn sdf_record_to_ast(
    sdf_data: &[u8],
    options: &[u8],
    record: usize,
) -> Result<Vec<u8>, String> {
    let sdf_str = std::str::from_utf8(sdf_data).map_err(|e| e.to_string())?;
    let commands = sdf_record_to_commands(sdf_str, mode_from_options(options), record)?;
    commands_to_cbor(&commands)
}

pub fn sdf_to_commands(sdf: &str, mode: RenderMode) -> Result<Vec<Command>, String> {
    sdf_record_to_commands(sdf, mode, 1)
}

pub fn sdf_record_to_commands(
    sdf: &str,
    mode: RenderMode,
    record: usize,
) -> Result<Vec<Command>, String> {
    let record_sdf = select_sdf_record(sdf, record)?;
    let mol = parse_sdf_auto_string(record_sdf)
        .map_err(|error| format!("could not parse SDF record {record}: {error}"))?;
    validate_sdf_molecule(&mol, record)?;
    let stereo_map = extract_sdf_stereo(&mol);
    let render_mol = render_mol_from_sdf(&mol, record_sdf);
    Ok(ast_from_render_mol(&render_mol, mode.as_str(), &stereo_map))
}

fn select_sdf_record(sdf: &str, record: usize) -> Result<&str, String> {
    if record == 0 {
        return Err("SDF record numbers are one-based".to_string());
    }
    if sdf.trim().is_empty() {
        return Err("SDF input is empty".to_string());
    }

    let records = split_sdf_records(sdf);
    let selected = records.get(record - 1).copied().ok_or_else(|| {
        format!(
            "SDF record {record} does not exist; input contains {} record(s)",
            records.len()
        )
    })?;
    if selected.trim().is_empty() {
        return Err(format!("SDF record {record} is empty"));
    }
    Ok(selected)
}

fn split_sdf_records(sdf: &str) -> Vec<&str> {
    let mut records = Vec::new();
    let mut record_start = 0;
    let mut offset = 0;

    for line_with_ending in sdf.split_inclusive('\n') {
        let line = line_with_ending.trim_end_matches(&['\r', '\n'][..]);
        if line.trim() == "$$$$" {
            records.push(&sdf[record_start..offset]);
            record_start = offset + line_with_ending.len();
        }
        offset += line_with_ending.len();
    }

    if record_start < sdf.len() {
        let trailing = &sdf[record_start..];
        if !trailing.trim().is_empty() || records.is_empty() {
            records.push(trailing);
        }
    } else if records.is_empty() {
        records.push(sdf);
    }

    records
}

fn validate_sdf_molecule(mol: &sdfrust::Molecule, record: usize) -> Result<(), String> {
    if mol.atoms.is_empty() {
        return Err(format!("SDF record {record} contains no atoms"));
    }

    for (index, atom) in mol.atoms.iter().enumerate() {
        if atom.element.trim().is_empty() {
            return Err(format!(
                "SDF record {record} atom {} has no element symbol",
                index + 1
            ));
        }
        if !atom.x.is_finite() || !atom.y.is_finite() || !atom.z.is_finite() {
            return Err(format!(
                "SDF record {record} atom {} has non-finite coordinates",
                index + 1
            ));
        }
    }

    for (index, bond) in mol.bonds.iter().enumerate() {
        if bond.atom1 >= mol.atoms.len() || bond.atom2 >= mol.atoms.len() {
            return Err(format!(
                "SDF record {record} bond {} refers to an atom outside the record",
                index + 1
            ));
        }
    }

    Ok(())
}

pub fn smiles_to_layout_input(smiles_data: &[u8]) -> Result<Vec<u8>, String> {
    let smiles = std::str::from_utf8(smiles_data).map_err(|e| e.to_string())?;
    smiles_layout_input(smiles, RenderMode::Abbreviate)
}

pub fn smiles_to_full_layout_input(smiles_data: &[u8]) -> Result<Vec<u8>, String> {
    let smiles = std::str::from_utf8(smiles_data).map_err(|e| e.to_string())?;
    smiles_layout_input(smiles, RenderMode::Full)
}

pub fn smiles_layout_input(smiles: &str, mode: RenderMode) -> Result<Vec<u8>, String> {
    let graph = parse_smiles_graph(smiles)?;
    if mode == RenderMode::Full {
        Ok(encode_layout_input(&expand_smiles_graph_hydrogens(&graph)))
    } else {
        Ok(encode_layout_input(&graph))
    }
}

pub fn smiles_to_ast(
    smiles_data: &[u8],
    coords_data: &[u8],
    options: &[u8],
) -> Result<Vec<u8>, String> {
    let smiles = std::str::from_utf8(smiles_data).map_err(|e| e.to_string())?;
    let commands = smiles_to_commands_with_coords(smiles, coords_data, mode_from_options(options))?;
    commands_to_cbor(&commands)
}

pub fn smiles_to_commands_with_coords(
    smiles: &str,
    coords_data: &[u8],
    mode: RenderMode,
) -> Result<Vec<Command>, String> {
    let mut graph = parse_smiles_graph(smiles)?;
    if mode == RenderMode::Full {
        graph = expand_smiles_graph_hydrogens(&graph);
    }
    let coords = decode_coords(coords_data, graph.atoms.len(), graph.bonds.len())?;
    let render_mol = render_mol_from_smiles(&graph, &coords.coords);
    let stereo_map = build_smiles_stereo_map(&graph, &coords.bond_styles);
    Ok(ast_from_render_mol(&render_mol, mode.as_str(), &stereo_map))
}

fn mode_from_options(options: &[u8]) -> RenderMode {
    std::str::from_utf8(options)
        .ok()
        .and_then(|mode| RenderMode::parse(mode).ok())
        .unwrap_or_default()
}

fn commands_to_cbor(commands: &[Command]) -> Result<Vec<u8>, String> {
    let mut buffer = Vec::new();
    ciborium::into_writer(commands, &mut buffer).map_err(|e| e.to_string())?;
    Ok(buffer)
}

fn render_mol_from_sdf(mol: &sdfrust::Molecule, sdf: &str) -> RenderMol {
    let metadata = extract_sdf_atom_metadata(sdf, mol);
    let atoms = mol
        .atoms
        .iter()
        .zip(metadata)
        .map(|(atom, metadata)| RenderAtom {
            element: atom.element.clone(),
            x: atom.x,
            y: atom.y,
            hydrogens: 0,
            charge: atom.formal_charge,
            isotope: metadata.isotope,
            radical: atom.radical.or(metadata.radical),
            atom_map: atom.atom_atom_mapping.or(metadata.atom_map),
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

fn extract_sdf_atom_metadata(sdf: &str, mol: &sdfrust::Molecule) -> Vec<SdfAtomMetadata> {
    let mut metadata = vec![SdfAtomMetadata::default(); mol.atoms.len()];
    if mol.format_version == SdfFormat::V3000 {
        for (metadata, atom) in metadata.iter_mut().zip(&mol.atoms) {
            metadata.isotope = u16::try_from(atom.mass_difference)
                .ok()
                .filter(|&mass| mass > 0);
        }
        extract_v3000_atom_metadata(sdf, &mut metadata);
        return metadata;
    }

    let lines = sdf.lines().collect::<Vec<_>>();
    let Some(counts_line) = lines.get(3) else {
        return metadata;
    };
    let atom_count = counts_line
        .get(0..3)
        .and_then(|value| value.trim().parse::<usize>().ok())
        .unwrap_or(0)
        .min(metadata.len());

    for (atom_idx, (metadata, atom)) in metadata
        .iter_mut()
        .zip(&mol.atoms)
        .take(atom_count)
        .enumerate()
    {
        let Some(line) = lines.get(4 + atom_idx) else {
            break;
        };

        let mass_difference = line
            .get(34..36)
            .and_then(|value| value.trim().parse::<i16>().ok())
            .unwrap_or(0);
        if mass_difference != 0 {
            metadata.isotope = nominal_isotope(&atom.element)
                .and_then(|mass| mass.checked_add(mass_difference))
                .and_then(|mass| u16::try_from(mass).ok());
        }

        let charge_code = line
            .get(36..39)
            .and_then(|value| value.trim().parse::<u8>().ok())
            .unwrap_or(0);
        if charge_code == 4 {
            metadata.radical = Some(2);
        }

        metadata.atom_map = line
            .get(60..63)
            .and_then(|value| value.trim().parse::<u32>().ok())
            .filter(|&value| value > 0);
    }

    for line in &lines {
        if line.starts_with("M  END") || line.starts_with("$$$$") {
            break;
        } else if line.starts_with("M  ISO") {
            for (atom_idx, value) in parse_mdl_property_pairs(line) {
                if let (Some(atom), Ok(isotope)) =
                    (metadata.get_mut(atom_idx), u16::try_from(value))
                {
                    atom.isotope = (isotope > 0).then_some(isotope);
                }
            }
        } else if line.starts_with("M  RAD") {
            for (atom_idx, value) in parse_mdl_property_pairs(line) {
                if let (Some(atom), Ok(radical)) = (metadata.get_mut(atom_idx), u8::try_from(value))
                {
                    atom.radical = (radical > 0).then_some(radical);
                }
            }
        }
    }

    metadata
}

fn extract_v3000_atom_metadata(sdf: &str, metadata: &mut [SdfAtomMetadata]) {
    let mut in_atom_block = false;
    let mut atom_index = 0;

    for line in v3000_logical_lines(sdf) {
        if line == "BEGIN ATOM" {
            in_atom_block = true;
            continue;
        }
        if line == "END ATOM" {
            break;
        }
        if !in_atom_block {
            continue;
        }

        let Some(atom_metadata) = metadata.get_mut(atom_index) else {
            break;
        };
        atom_index += 1;
        let parts = line.split_whitespace().collect::<Vec<_>>();
        atom_metadata.atom_map = parts
            .get(5)
            .and_then(|value| value.parse::<u32>().ok())
            .filter(|&value| value > 0);

        for part in parts.iter().skip(6) {
            let Some((key, value)) = part.split_once('=') else {
                continue;
            };
            match key {
                "MASS" => {
                    atom_metadata.isotope =
                        value.parse::<u16>().ok().filter(|&isotope| isotope > 0);
                }
                "RAD" => {
                    atom_metadata.radical = value.parse::<u8>().ok().filter(|&radical| radical > 0);
                }
                _ => {}
            }
        }
    }
}

fn v3000_logical_lines(sdf: &str) -> Vec<String> {
    let mut lines = Vec::new();
    let mut current = String::new();

    for line in sdf.lines() {
        let Some(content) = line.strip_prefix("M  V30 ") else {
            continue;
        };
        if let Some(content) = content.strip_suffix('-') {
            current.push_str(content);
        } else {
            current.push_str(content);
            lines.push(std::mem::take(&mut current));
        }
    }

    if !current.is_empty() {
        lines.push(current);
    }
    lines
}

fn nominal_isotope(element: &str) -> Option<i16> {
    let element = element.parse::<AtomSymbol>().ok()?;
    let mass = element.standard_mass();
    (mass > 0.0 && mass.is_finite()).then_some(mass.round() as i16)
}

fn parse_mdl_property_pairs(line: &str) -> Vec<(usize, i32)> {
    let count = line
        .get(6..9)
        .and_then(|value| value.trim().parse::<usize>().ok())
        .unwrap_or(0);
    let mut pairs = Vec::with_capacity(count);

    for index in 0..count {
        let start = 9 + index * 8;
        let Some(atom_number) = line
            .get(start..start + 4)
            .and_then(|value| value.trim().parse::<usize>().ok())
        else {
            break;
        };
        let Some(value) = line
            .get(start + 4..start + 8)
            .and_then(|value| value.trim().parse::<i32>().ok())
        else {
            break;
        };
        if atom_number > 0 {
            pairs.push((atom_number - 1, value));
        }
    }

    pairs
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
            isotope: atom.isotope,
            radical: None,
            atom_map: atom.atom_map,
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
                isotope: None,
                atom_map: None,
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
            // Coordgen does not have a chemical atomic number for the SMILES
            // wildcard. Carbon provides stable neutral geometry while the
            // original `*` symbol remains in the rendering model.
            let layout_atomic_number = if atomic_number == 0 { 6 } else { atomic_number };

            SmilesAtom {
                element: node.atom().element().to_string(),
                atomic_number: layout_atomic_number,
                hydrogens: node.hydrogens(),
                charge: node.atom().charge(),
                isotope: node.atom().isotope(),
                atom_map: node.class().map(u32::from),
                stereo: atom_stereo_spec(atom_idx, node.chirality(), &incident_bonds),
                stereo_annotation: chirality_annotation(node.chirality()),
            }
        })
        .collect::<Vec<_>>();

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
    chirality: Option<SmilesChirality>,
    incident_bonds: &[Vec<(usize, usize)>],
) -> Option<AtomStereoSpec> {
    let direction = match chirality {
        Some(SmilesChirality::TH1) => AtomStereoDirection::CounterClockwise,
        Some(SmilesChirality::TH2) => AtomStereoDirection::Clockwise,
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

fn chirality_annotation(chirality: Option<SmilesChirality>) -> Option<String> {
    match chirality {
        Some(SmilesChirality::TH1 | SmilesChirality::TH2) | None => None,
        Some(chirality) => Some(chirality.to_string()),
    }
}

fn build_double_bond_stereo_specs(
    molecule: &SmilesMolecule,
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
    molecule: &SmilesMolecule,
    center: usize,
    other_center: usize,
) -> Option<(usize, i8)> {
    molecule.bonds().iter().find_map(|bond| match bond.kind() {
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
                let other = if bond.atom1 == atom {
                    bond.atom2
                } else {
                    bond.atom1
                };
                if !visited[other] {
                    visited[other] = true;
                    queue.push_back(other);
                }
            }
        }

        let mut edge_indices = component_edges.into_iter().collect::<Vec<_>>();
        edge_indices.sort_unstable_by_key(|&bond_idx| {
            let bond = &bonds[bond_idx];
            (
                usize::MAX - (aromatic_adj[bond.atom1].len() + aromatic_adj[bond.atom2].len()),
                bond_idx,
            )
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
    let atom_stereo_count = graph
        .atoms
        .iter()
        .filter(|atom| atom.stereo.is_some())
        .count();
    let mut buffer = Vec::with_capacity(
        20 + graph.atoms.len()
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
        write_u32_le(&mut buffer, stereo.looking_from.unwrap_or(usize::MAX));
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

    Ok(CoordPayload {
        coords,
        bond_styles,
    })
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
) -> Vec<Command> {
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
    let context = DfsContext {
        adj: &adj,
        mol,
        labels: &labels,
        stereo_map,
        rings: &rings,
        avg_length,
    };

    let mut root_commands = Vec::new();
    for start_node in 0..mol.atoms.len() {
        if !visited_nodes[start_node] {
            dfs(
                start_node,
                &context,
                &mut visited_nodes,
                &mut handled_bonds,
                &mut root_commands,
            );
        }
    }

    root_commands
}

fn build_labels(
    mol: &RenderMol,
    adj: &[Vec<(usize, BondKind, usize)>],
    mode_str: &str,
    visited_nodes: &mut [bool],
    handled_bonds: &mut [bool],
) -> Vec<RenderLabel> {
    let mut labels = vec![RenderLabel::default(); mol.atoms.len()];

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
                    if atom.stereo_annotation.is_some()
                        || atom.charge != 0
                        || has_visible_atom_metadata(atom)
                    {
                        labels[i] = render_label(atom, total_h);
                    }
                } else {
                    let heavy_atoms = adj[i].len().saturating_sub(explicit_h_count as usize);
                    labels[i] = if (heavy_atoms <= 1 && total_h > 0)
                        || atom.stereo_annotation.is_some()
                        || atom.charge != 0
                        || has_visible_atom_metadata(atom)
                    {
                        render_label(atom, total_h)
                    } else {
                        RenderLabel::default()
                    };
                }
            } else {
                labels[i] = render_label(atom, total_h);
            }
        } else {
            labels[i] = render_label(atom, 0);
        }
    }

    labels
}

fn has_visible_atom_metadata(atom: &RenderAtom) -> bool {
    atom.element == "*"
        || atom.isotope.is_some()
        || atom.radical.is_some()
        || atom.atom_map.is_some()
}

fn render_label(atom: &RenderAtom, hydrogen_count: u8) -> RenderLabel {
    let text = atom_label(&atom.element, hydrogen_count, atom.charge);
    let structured = has_visible_atom_metadata(atom).then(|| AtomLabel {
        symbol: atom.element.clone(),
        hydrogen_count,
        charge: atom.charge,
        isotope: atom.isotope,
        radical: atom.radical,
        atom_map: atom.atom_map,
    });

    RenderLabel {
        text,
        atom: structured,
    }
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

fn extract_sdf_stereo(mol: &sdfrust::Molecule) -> HashMap<(usize, usize), (u8, bool)> {
    let mut map = HashMap::new();
    for bond in &mol.bonds {
        let stereo = match bond.stereo {
            BondStereo::Up => 1,
            BondStereo::Down => 6,
            BondStereo::None | BondStereo::Either => continue,
        };
        if bond.atom1 < mol.atoms.len() && bond.atom2 < mol.atoms.len() {
            map.insert((bond.atom1, bond.atom2), (stereo, true));
            map.insert((bond.atom2, bond.atom1), (stereo, false));
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
    let ring = rings
        .iter()
        .filter(|ring| ring_contains_edge(ring, u, v))
        .min_by_key(|ring| ring.len())?;

    let (cx, cy) = ring_centroid(ring, mol);

    let a = &mol.atoms[u];
    let b = &mol.atoms[v];
    let dx = b.x - a.x;
    let dy = b.y - a.y;
    let cx_dir = cx - a.x;
    let cy_dir = cy - a.y;
    let cross = dx * cy_dir - dy * cx_dir;

    Some(if cross > 0.0 { "left" } else { "right" }.to_string())
}

fn ring_contains_edge(ring: &[usize], u: usize, v: usize) -> bool {
    ring.iter().enumerate().any(|(idx, &node)| {
        let next = ring[(idx + 1) % ring.len()];
        (node == u && next == v) || (node == v && next == u)
    })
}

fn ring_centroid(ring: &[usize], mol: &RenderMol) -> (f64, f64) {
    if ring.len() < 3 {
        return average_ring_center(ring, mol);
    }

    let mut twice_area = 0.0;
    let mut cx = 0.0;
    let mut cy = 0.0;

    for idx in 0..ring.len() {
        let a = &mol.atoms[ring[idx]];
        let b = &mol.atoms[ring[(idx + 1) % ring.len()]];
        let cross = a.x * b.y - b.x * a.y;
        twice_area += cross;
        cx += (a.x + b.x) * cross;
        cy += (a.y + b.y) * cross;
    }

    if twice_area.abs() < 1e-9 {
        return average_ring_center(ring, mol);
    }

    (cx / (3.0 * twice_area), cy / (3.0 * twice_area))
}

fn average_ring_center(ring: &[usize], mol: &RenderMol) -> (f64, f64) {
    let mut cx = 0.0;
    let mut cy = 0.0;
    for &node in ring {
        cx += mol.atoms[node].x;
        cy += mol.atoms[node].y;
    }
    (cx / ring.len() as f64, cy / ring.len() as f64)
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

struct DfsContext<'a> {
    adj: &'a [Vec<(usize, BondKind, usize)>],
    mol: &'a RenderMol,
    labels: &'a [RenderLabel],
    stereo_map: &'a HashMap<(usize, usize), (u8, bool)>,
    rings: &'a [Vec<usize>],
    avg_length: f64,
}

fn dfs(
    u: usize,
    context: &DfsContext<'_>,
    visited_nodes: &mut [bool],
    handled_bonds: &mut [bool],
    commands: &mut Vec<Command>,
) {
    visited_nodes[u] = true;
    let u_atom = &context.mol.atoms[u];

    let mut back_edges = Vec::new();
    for &(v, kind, bond_idx) in &context.adj[u] {
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
            calculate_double_bond_offset(u, v, context.mol, context.rings)
        } else {
            None
        };
        links.push(LinkData {
            target: format!("a{v}"),
            name: format!("b{bond_idx}"),
            bond_type: bond_func_name(kind, u, v, context.stereo_map).to_string(),
            angle: calc_angle(u_atom, &context.mol.atoms[v]),
            offset,
            length_scale: calc_length_scale(u_atom, &context.mol.atoms[v], context.avg_length),
        });
    }

    commands.push(Command::Fragment {
        element: context.labels[u].text.clone(),
        name: format!("a{u}"),
        links,
        atom: context.labels[u].atom.clone(),
        annotation: u_atom.stereo_annotation.clone(),
    });

    let mut forward_targets = Vec::new();
    for &(v, kind, bond_idx) in &context.adj[u] {
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
            calculate_double_bond_offset(u, v, context.mol, context.rings)
        } else {
            None
        };

        let bond_cmd = Command::Bond {
            name: format!("b{bond_idx}"),
            bond_type: bond_func_name(kind, u, v, context.stereo_map).to_string(),
            angle: calc_angle(u_atom, &context.mol.atoms[v]),
            offset,
            length_scale: calc_length_scale(u_atom, &context.mol.atoms[v], context.avg_length),
        };

        if has_more {
            let mut branch_body = vec![bond_cmd];
            dfs(v, context, visited_nodes, handled_bonds, &mut branch_body);
            commands.push(Command::Branch { body: branch_body });
        } else {
            commands.push(bond_cmd);
            dfs(v, context, visited_nodes, handled_bonds, commands);
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

    type EncodedBond = (u32, u32, u8);
    type EncodedAtomStereo = (u32, u32, u32, u32, u8);
    type EncodedDoubleBondStereo = (u32, u32, u32, u8);
    type LayoutInput = (
        Vec<u8>,
        Vec<EncodedBond>,
        Vec<EncodedAtomStereo>,
        Vec<EncodedDoubleBondStereo>,
    );

    fn decode_layout_input(bytes: &[u8]) -> LayoutInput {
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

    fn coordinate_payload(points: &[(f32, f32)], bond_count: usize) -> Vec<u8> {
        let mut coords = Vec::new();
        coords.extend_from_slice(COORD_MAGIC);
        coords.extend_from_slice(&(points.len() as u32).to_le_bytes());
        coords.extend_from_slice(&(bond_count as u32).to_le_bytes());
        for &(x, y) in points {
            coords.extend_from_slice(&x.to_le_bytes());
            coords.extend_from_slice(&y.to_le_bytes());
        }
        coords.extend(std::iter::repeat_n(0, bond_count));
        coords
    }

    fn first_fragment(commands: &[Command]) -> (&str, Option<&AtomLabel>) {
        match commands.first().unwrap() {
            Command::Fragment { element, atom, .. } => (element, atom.as_ref()),
            command => panic!("expected fragment, got {command:?}"),
        }
    }

    const CARBON_V2000: &str = concat!(
        "carbon\n",
        "  molchemist\n",
        "\n",
        "  1  0  0  0  0  0  0  0  0  0999 V2000\n",
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
        "M  END\n",
    );

    const OXYGEN_V2000: &str = concat!(
        "oxygen\n",
        "  molchemist\n",
        "\n",
        "  1  0  0  0  0  0  0  0  0  0999 V2000\n",
        "    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n",
        "M  END\n",
    );

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
    fn sdf_preserves_charge_isotope_radical_and_atom_map() {
        let sdf = concat!(
            "metadata\n",
            "  molchemist\n",
            "\n",
            "  1  0  0  0  0  0  0  0  0  0999 V2000\n",
            "    0.0000    0.0000    0.0000 C   0  3  0  0  0  0  0  0  0  7  0  0\n",
            "M  ISO  1   1  13\n",
            "M  RAD  1   1   2\n",
            "M  END\n",
        );

        let commands = sdf_to_commands(sdf, RenderMode::Full).unwrap();
        let (label, atom) = first_fragment(&commands);
        let atom = atom.unwrap();

        assert_eq!(label, "C^+");
        assert_eq!(atom.symbol, "C");
        assert_eq!(atom.charge, 1);
        assert_eq!(atom.isotope, Some(13));
        assert_eq!(atom.radical, Some(2));
        assert_eq!(atom.atom_map, Some(7));
    }

    #[test]
    fn sdf_atom_block_mass_difference_becomes_an_isotope() {
        let sdf = concat!(
            "isotope\n",
            "  molchemist\n",
            "\n",
            "  1  0  0  0  0  0  0  0  0  0999 V2000\n",
            "    0.0000    0.0000    0.0000 C   1  0  0  0  0  0  0  0  0  0  0  0\n",
            "M  END\n",
        );

        let commands = sdf_to_commands(sdf, RenderMode::Full).unwrap();
        let (_, atom) = first_fragment(&commands);

        assert_eq!(atom.unwrap().isotope, Some(13));
    }

    #[test]
    fn sdf_auto_detects_v3000_and_preserves_atom_metadata() {
        let sdf = concat!(
            "v3000 metadata\n",
            "  molchemist\n",
            "\n",
            "  0  0  0     0  0            999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 2 1 0 0 0\n",
            "M  V30 BEGIN ATOM\n",
            "M  V30 1 U 0.0000 0.0000 0.0000 7 CHG=1 MASS=238 RAD=2\n",
            "M  V30 2 O 1.5000 0.0000 0.0000 0 CHG=-1\n",
            "M  V30 END ATOM\n",
            "M  V30 BEGIN BOND\n",
            "M  V30 1 1 1 2 CFG=1\n",
            "M  V30 END BOND\n",
            "M  V30 END CTAB\n",
            "M  END\n",
        );

        let commands = sdf_to_commands(sdf, RenderMode::Full).unwrap();
        let (label, atom) = first_fragment(&commands);
        let atom = atom.unwrap();

        assert_eq!(label, "U^+");
        assert_eq!(atom.symbol, "U");
        assert_eq!(atom.charge, 1);
        assert_eq!(atom.isotope, Some(238));
        assert_eq!(atom.radical, Some(2));
        assert_eq!(atom.atom_map, Some(7));
        assert!(commands.iter().any(|command| matches!(
            command,
            Command::Bond { bond_type, .. } if bond_type.starts_with("cram-filled")
        )));
    }

    #[test]
    fn sdf_record_selection_is_one_based_and_format_independent() {
        let sdf = format!("{CARBON_V2000}$$$$\r\n{OXYGEN_V2000}$$$$\r\n\r\n");

        let commands = sdf_record_to_commands(&sdf, RenderMode::Full, 2).unwrap();
        assert_eq!(first_fragment(&commands).0, "O");
        assert_eq!(
            sdf_record_to_commands(&sdf, RenderMode::Full, 0).unwrap_err(),
            "SDF record numbers are one-based"
        );
        assert_eq!(
            sdf_record_to_commands(&sdf, RenderMode::Full, 3).unwrap_err(),
            "SDF record 3 does not exist; input contains 2 record(s)"
        );
    }

    #[test]
    fn sdf_rejects_empty_structures_and_non_finite_coordinates() {
        let empty = concat!(
            "empty\n",
            "  molchemist\n",
            "\n",
            "  0  0  0  0  0  0  0  0  0  0999 V2000\n",
            "M  END\n",
        );
        assert_eq!(
            sdf_to_commands(empty, RenderMode::Full).unwrap_err(),
            "SDF record 1 contains no atoms"
        );

        let non_finite = CARBON_V2000.replacen("    0.0000", "       NaN", 1);
        assert_eq!(
            sdf_to_commands(&non_finite, RenderMode::Full).unwrap_err(),
            "SDF record 1 atom 1 has non-finite coordinates"
        );
    }

    #[test]
    fn smiles_preserves_isotope_and_atom_class() {
        let coords = coordinate_payload(&[(0.0, 0.0), (1.0, 0.0)], 1);
        let commands =
            smiles_to_commands_with_coords("[13CH3:7]C", &coords, RenderMode::Abbreviate).unwrap();
        let (label, atom) = first_fragment(&commands);
        let atom = atom.unwrap();

        assert_eq!(label, "CH_3");
        assert_eq!(atom.symbol, "C");
        assert_eq!(atom.hydrogen_count, 3);
        assert_eq!(atom.isotope, Some(13));
        assert_eq!(atom.atom_map, Some(7));
    }

    #[test]
    fn skeletal_mode_keeps_charged_carbon_visible() {
        let coords = coordinate_payload(&[(0.0, 0.0), (1.0, 0.0)], 1);
        let commands =
            smiles_to_commands_with_coords("[CH2-]C", &coords, RenderMode::Skeletal).unwrap();
        let (label, atom) = first_fragment(&commands);

        assert_eq!(label, "CH_2^-");
        assert!(atom.is_none());
    }

    #[test]
    fn smiles_wildcard_uses_neutral_layout_surrogate_and_visible_label() {
        let payload = smiles_to_layout_input(b"*").unwrap();
        let (atoms, bonds, _, _) = decode_layout_input(&payload);
        assert_eq!(atoms, vec![6]);
        assert!(bonds.is_empty());

        let coords = coordinate_payload(&[(0.0, 0.0)], 0);
        let commands = smiles_to_commands_with_coords("*", &coords, RenderMode::Skeletal).unwrap();
        let (label, atom) = first_fragment(&commands);

        assert_eq!(label, "*");
        assert_eq!(atom.unwrap().symbol, "*");
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
        assert_eq!(
            bonds.iter().map(|(_, _, order)| *order).collect::<Vec<_>>(),
            vec![2, 1, 2, 1, 2, 1]
        );
    }

    #[test]
    fn branch_leading_bond_type_does_not_leak_into_smiles_layout() {
        let payload = smiles_to_layout_input(b"C(=CC=O)C=C(C(=O)O)N").unwrap();
        let (_, bonds, _, _) = decode_layout_input(&payload);

        assert_eq!(bonds[0], (0, 1, 2));
        assert_eq!(bonds[1], (1, 2, 1));
        assert_eq!(bonds[2], (2, 3, 2));
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

    fn test_atom(element: &str, x: f64, y: f64) -> RenderAtom {
        RenderAtom {
            element: element.to_string(),
            x,
            y,
            hydrogens: 0,
            charge: 0,
            isotope: None,
            radical: None,
            atom_map: None,
            stereo_annotation: None,
        }
    }

    #[test]
    fn double_bond_offset_uses_ring_where_bond_is_an_edge() {
        let mol = RenderMol {
            atoms: vec![
                test_atom("C", 0.0, 0.0),
                test_atom("C", 1.0, 0.0),
                test_atom("C", 1.0, 1.0),
                test_atom("C", 0.0, 1.0),
                test_atom("C", 0.5, -1.0),
                test_atom("C", -0.5, -1.0),
            ],
            bonds: Vec::new(),
        };
        let rings = vec![vec![0, 4, 1, 5], vec![0, 1, 2, 3]];

        assert_eq!(
            calculate_double_bond_offset(0, 1, &mol, &rings).as_deref(),
            Some("left")
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
        for (x, y) in [
            (0.0f32, 0.0f32),
            (1.0, 0.0),
            (-1.0, 0.0),
            (0.0, 1.0),
            (0.0, -1.0),
        ] {
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
