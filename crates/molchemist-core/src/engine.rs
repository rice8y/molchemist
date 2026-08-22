use crate::semantic::{contractible_superatom_indexes, SdfAtomMetadata};
use crate::{
    parse as parse_smiles, AtomSymbol, BondType as SmilesBondType, ChemicalBondOrder,
    ChemicalRecord, Chirality as SmilesChirality, Molecule as SmilesMolecule,
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
    /// Reset placement before rendering the next disconnected component.
    #[serde(rename = "component-break")]
    ComponentBreak,
    /// CTfile display metadata rendered as an overlay after the molecular skeleton.
    #[serde(rename = "ctfile")]
    Ctfile {
        sgroups: Vec<CtfileSGroupDepiction>,
        highlights: Vec<CtfileHighlightDepiction>,
        #[serde(rename = "atomQueries")]
        atom_queries: Vec<CtfileAtomQueryDepiction>,
        #[serde(rename = "bondQueries")]
        bond_queries: Vec<CtfileBondQueryDepiction>,
        #[serde(rename = "atomAnnotations")]
        atom_annotations: Vec<CtfileAtomAnnotationDepiction>,
        #[serde(rename = "variableAttachments")]
        variable_attachments: Vec<CtfileVariableAttachmentDepiction>,
    },
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
    #[serde(rename = "rgroupLabel", skip_serializing_if = "Option::is_none")]
    pub rgroup_label: Option<u8>,
    #[serde(rename = "querySymbol", default)]
    pub query_symbol: bool,
    #[serde(rename = "queryNegated", default)]
    pub query_negated: bool,
    #[serde(default)]
    pub hidden: bool,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct CtfilePoint {
    pub atom_index: usize,
    pub offset: [f64; 2],
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct CtfileSGroupDepiction {
    pub id: u32,
    pub kind: String,
    pub left_top: CtfilePoint,
    pub left_bottom: CtfilePoint,
    pub right_top: CtfilePoint,
    pub right_bottom: CtfilePoint,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct CtfileHighlightBond {
    pub bond_index: usize,
    pub atom1_index: usize,
    pub atom2_index: usize,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct CtfileHighlightDepiction {
    pub atom_indexes: Vec<usize>,
    pub bonds: Vec<CtfileHighlightBond>,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct CtfileBondQueryDepiction {
    pub bond_index: usize,
    pub label: String,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct CtfileAtomQueryDepiction {
    pub atom_index: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub symbol: Option<String>,
    pub label: String,
    pub compact_label: String,
    pub constraint_lines: Vec<String>,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct CtfileAtomAnnotationDepiction {
    pub atom_index: usize,
    pub kind: String,
    pub label: String,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct CtfileVariableAttachmentDepiction {
    pub bond_index: usize,
    pub endpoint_atom_indexes: Vec<usize>,
    pub mode: String,
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
    Quadruple,
    Aromatic,
    SingleOrDouble,
    SingleOrAromatic,
    DoubleOrAromatic,
    Any,
    Coordination,
    Hydrogen,
}

impl BondKind {
    fn has_parallel_line(self) -> bool {
        matches!(
            self,
            Self::Double | Self::Aromatic | Self::SingleOrDouble | Self::DoubleOrAromatic
        )
    }

    fn contributes_to_average_length(self) -> bool {
        self != Self::Hydrogen
    }
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
    rgroup_label: Option<u8>,
    query_symbol: bool,
    query_negated: bool,
    stereo_annotation: Option<String>,
    force_label: bool,
    hide_label: bool,
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

struct RecordRenderProjection {
    mol: RenderMol,
    /// Original semantic atom index to projected depiction atom index.
    atom_indexes: Vec<usize>,
    /// Original semantic bond index to projected depiction bond index. Internal
    /// bonds of a contracted superatom deliberately have no projected bond.
    bond_indexes: Vec<Option<usize>>,
    contracted_sgroup_ids: HashSet<u32>,
    contracted_atom_indexes: HashSet<usize>,
    coordinates_preserved: bool,
}

#[derive(Clone, Default)]
struct RenderLabel {
    text: String,
    atom: Option<AtomLabel>,
}

#[derive(Clone, Copy, PartialEq, Eq)]
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
    extended_stereo: Option<ExtendedAtomStereoSpec>,
    stereo_annotation: Option<String>,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum SmilesBondKind {
    Single,
    Double,
    Triple,
    Quadruple,
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
    neighbor_order: Vec<Vec<usize>>,
    hydrogen_order_positions: Vec<usize>,
    double_bond_stereo: Vec<DoubleBondStereoSpec>,
}

#[derive(Clone, Copy)]
struct AtomStereoSpec {
    looking_from: Option<usize>,
    atom1: Option<usize>,
    atom2: Option<usize>,
    direction: AtomStereoDirection,
}

#[derive(Clone, Copy)]
enum PlanarShape {
    U,
    Four,
    Z,
}

#[derive(Clone)]
enum ExtendedAtomStereoSpec {
    Allene {
        direction: AtomStereoDirection,
    },
    SquarePlanar {
        ligands: [usize; 4],
        shape: PlanarShape,
    },
    TrigonalBipyramidal {
        axis_from: usize,
        axis_to: usize,
        equatorial: [usize; 3],
        direction: AtomStereoDirection,
    },
    Octahedral {
        axis_from: usize,
        axis_to: usize,
        equatorial: [usize; 4],
        shape: PlanarShape,
        direction: AtomStereoDirection,
    },
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

pub fn inspect_sdf_record(sdf: &str, record: usize) -> Result<ChemicalRecord, String> {
    Ok(parse_sdf_record(sdf, record)?.chemical)
}

pub fn sdf_record_to_inspection_cbor(sdf_data: &[u8], record: usize) -> Result<Vec<u8>, String> {
    let sdf = std::str::from_utf8(sdf_data).map_err(|error| error.to_string())?;
    let inspection = inspect_sdf_record(sdf, record)?;
    let mut buffer = Vec::new();
    ciborium::into_writer(&inspection, &mut buffer).map_err(|error| error.to_string())?;
    Ok(buffer)
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

pub fn sdf_record_to_layout_input(sdf_data: &[u8], record: usize) -> Result<Vec<u8>, String> {
    let sdf = std::str::from_utf8(sdf_data).map_err(|error| error.to_string())?;
    let parsed = parse_sdf_record(sdf, record)?;
    if sdf_requires_layout(&parsed.molecule) {
        Ok(encode_layout_input(&sdf_layout_graph(&parsed.molecule)))
    } else {
        Ok(Vec::new())
    }
}

pub fn sdf_record_to_ast_with_coords(
    sdf_data: &[u8],
    coords_data: &[u8],
    options: &[u8],
    record: usize,
) -> Result<Vec<u8>, String> {
    let sdf = std::str::from_utf8(sdf_data).map_err(|error| error.to_string())?;
    let commands =
        sdf_record_to_commands_with_coords(sdf, coords_data, mode_from_options(options), record)?;
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
    let parsed = parse_sdf_record(sdf, record)?;
    let stereo_map = extract_sdf_stereo(&parsed.molecule, parsed.raw_record);
    let stereo_annotations = extract_sdf_stereo_annotations(parsed.raw_record, &parsed.molecule);
    let render_mol = render_mol_from_record(&parsed.chemical, stereo_annotations);
    let projection = project_record_render(&parsed.chemical, render_mol);
    Ok(ast_from_chemical_record(
        &parsed.chemical,
        &projection,
        mode.as_str(),
        &stereo_map,
    ))
}

pub fn sdf_record_to_commands_with_coords(
    sdf: &str,
    coords_data: &[u8],
    mode: RenderMode,
    record: usize,
) -> Result<Vec<Command>, String> {
    let parsed = parse_sdf_record(sdf, record)?;
    let coords = decode_coords(
        coords_data,
        parsed.molecule.atoms.len(),
        parsed.molecule.bonds.len(),
    )?;
    let stereo_map = extract_sdf_stereo(&parsed.molecule, parsed.raw_record);
    let stereo_annotations = extract_sdf_stereo_annotations(parsed.raw_record, &parsed.molecule);
    let mut render_mol = render_mol_from_record(&parsed.chemical, stereo_annotations);
    for (atom, &(x, y)) in render_mol.atoms.iter_mut().zip(&coords.coords) {
        atom.x = f64::from(x);
        atom.y = f64::from(y);
    }
    let projection = project_record_render(&parsed.chemical, render_mol);
    Ok(ast_from_chemical_record(
        &parsed.chemical,
        &projection,
        mode.as_str(),
        &stereo_map,
    ))
}

struct ParsedSdfRecord<'a> {
    molecule: sdfrust::Molecule,
    raw_record: &'a str,
    chemical: ChemicalRecord,
}

fn parse_sdf_record(sdf: &str, record: usize) -> Result<ParsedSdfRecord<'_>, String> {
    let raw_record = select_sdf_record(sdf, record)?;
    let parser_record = normalize_v3000_atom_lists(raw_record);
    let molecule = parse_sdf_auto_string(&parser_record)
        .map_err(|error| format!("could not parse SDF record {record}: {error}"))?;
    validate_sdf_molecule(&molecule, record)?;
    let metadata = extract_sdf_atom_metadata(raw_record, &molecule);
    let chemical = ChemicalRecord::from_sdf(&molecule, raw_record, record, &metadata);
    Ok(ParsedSdfRecord {
        molecule,
        raw_record,
        chemical,
    })
}

fn normalize_v3000_atom_lists(raw_record: &str) -> String {
    let mut output = String::with_capacity(raw_record.len());
    let mut in_atom_block = false;
    for line_with_ending in raw_record.split_inclusive('\n') {
        let line = line_with_ending.trim_end_matches(['\r', '\n']);
        let ending = &line_with_ending[line.len()..];
        let Some(content) = line.strip_prefix("M  V30 ") else {
            output.push_str(line_with_ending);
            continue;
        };
        if content == "BEGIN ATOM" {
            in_atom_block = true;
            output.push_str(line_with_ending);
            continue;
        }
        if content == "END ATOM" {
            in_atom_block = false;
            output.push_str(line_with_ending);
            continue;
        }
        if !in_atom_block {
            output.push_str(line_with_ending);
            continue;
        }

        let Some(id_end) = content.find(char::is_whitespace) else {
            output.push_str(line_with_ending);
            continue;
        };
        let after_id = content[id_end..].trim_start();
        let list_start = after_id
            .strip_prefix('"')
            .unwrap_or(after_id)
            .strip_prefix("NOT ")
            .unwrap_or_else(|| after_id.strip_prefix('"').unwrap_or(after_id));
        if !list_start.starts_with('[') {
            output.push_str(line_with_ending);
            continue;
        }
        let Some(list_end) = list_start.find(']') else {
            output.push_str(line_with_ending);
            continue;
        };
        let suffix = list_start[list_end + 1..]
            .strip_prefix('"')
            .unwrap_or(&list_start[list_end + 1..]);
        output.push_str("M  V30 ");
        output.push_str(&content[..id_end]);
        output.push_str(" L");
        output.push_str(suffix);
        output.push_str(ending);
    }
    if !raw_record.ends_with('\n') && output.ends_with('\n') {
        output.pop();
    }
    output
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

fn sdf_requires_layout(mol: &sdfrust::Molecule) -> bool {
    if mol.bonds.is_empty() {
        return false;
    }

    let coordinate_scale = mol.atoms.iter().fold(1.0f64, |scale, atom| {
        scale.max(atom.x.abs()).max(atom.y.abs())
    });
    let minimum_length = coordinate_scale * f64::EPSILON * 64.0;
    let mut lengths = Vec::with_capacity(mol.bonds.len());

    for bond in &mol.bonds {
        let atom1 = &mol.atoms[bond.atom1];
        let atom2 = &mol.atoms[bond.atom2];
        let length = (atom2.x - atom1.x).hypot(atom2.y - atom1.y);
        if !length.is_finite() || length <= minimum_length {
            return true;
        }
        if bond.order != BondOrder::Hydrogen {
            lengths.push(length);
        }
    }

    if lengths.is_empty() {
        return false;
    }
    let typical_length = median(&mut lengths).unwrap_or(1.0);
    lengths
        .iter()
        .any(|length| *length < typical_length * 0.05 || *length > typical_length * 20.0)
}

fn sdf_layout_graph(mol: &sdfrust::Molecule) -> SmilesGraph {
    let atoms = mol
        .atoms
        .iter()
        .map(|atom| {
            let atomic_number = atom
                .element
                .parse::<AtomSymbol>()
                .map(|element| element.atomic_number())
                .ok()
                .filter(|&number| number > 0)
                .unwrap_or(6);
            SmilesAtom {
                element: atom.element.clone(),
                atomic_number,
                hydrogens: 0,
                charge: atom.formal_charge,
                isotope: None,
                atom_map: None,
                stereo: None,
                extended_stereo: None,
                stereo_annotation: None,
            }
        })
        .collect::<Vec<_>>();
    let bonds = mol
        .bonds
        .iter()
        .map(|bond| SmilesBond {
            atom1: bond.atom1,
            atom2: bond.atom2,
            kind: match bond.order {
                BondOrder::Double | BondOrder::DoubleOrAromatic => SmilesBondKind::Double,
                BondOrder::Triple => SmilesBondKind::Triple,
                BondOrder::Aromatic => SmilesBondKind::Aromatic,
                _ => SmilesBondKind::Single,
            },
        })
        .collect::<Vec<_>>();
    let mut neighbor_order = vec![Vec::new(); atoms.len()];
    for bond in &bonds {
        neighbor_order[bond.atom1].push(bond.atom2);
        neighbor_order[bond.atom2].push(bond.atom1);
    }

    SmilesGraph {
        atoms,
        bonds,
        neighbor_order,
        hydrogen_order_positions: vec![0; mol.atoms.len()],
        double_bond_stereo: Vec::new(),
    }
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
    let mut coords = decode_coords(coords_data, graph.atoms.len(), graph.bonds.len())?;
    let mut stereo_map = build_smiles_stereo_map(&graph, &coords.bond_styles);
    let depicted_centers =
        apply_extended_stereo_depictions(&graph, &mut coords.coords, &mut stereo_map);
    let render_mol = render_mol_from_smiles(&graph, &coords.coords, &depicted_centers);
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

fn render_mol_from_record(
    record: &ChemicalRecord,
    stereo_annotations: Vec<Option<String>>,
) -> RenderMol {
    let mut atoms = record
        .atoms
        .iter()
        .zip(stereo_annotations)
        .map(|(atom, stereo_annotation)| RenderAtom {
            element: atom.element.clone(),
            x: atom.coordinates[0],
            y: atom.coordinates[1],
            hydrogens: 0,
            charge: atom.formal_charge,
            isotope: atom.isotope,
            radical: atom.radical,
            atom_map: atom.atom_map,
            rgroup_label: None,
            query_symbol: false,
            query_negated: false,
            stereo_annotation,
            force_label: false,
            hide_label: false,
        })
        .collect::<Vec<_>>();

    for atom in &record.atoms {
        let mut label = atom.element.clone();
        if let Some(query) = &atom.query {
            if let Some(symbol) = atom_query_glyph_symbol(query) {
                label = symbol;
                atoms[atom.index].query_symbol = true;
                atoms[atom.index].query_negated = query.is_not_list;
            } else if matches!(atom.element.as_str(), "A" | "Q" | "*") {
                label = atom.element.clone();
            }
        }
        if label == "R#" {
            if let Some(rgroup) = atom.rgroup_labels.first().or(atom.rgroup_label.as_ref()) {
                label = "R".to_string();
                atoms[atom.index].rgroup_label = Some(*rgroup);
            }
        }
        atoms[atom.index].element = label;
        atoms[atom.index].force_label =
            atom.query.is_some() || atom.element == "R#" || !atom.rgroup_labels.is_empty();
        atoms[atom.index].hide_label = false;
    }

    for group in &record.sgroups {
        if group.kind == "superatom" && group.atom_source_ids.len() == 1 {
            let Some(label) = group.label.as_ref().or(group.subscript.as_ref()) else {
                continue;
            };
            let source_id = group.atom_source_ids[0];
            if let Some(atom) = record.atoms.iter().find(|atom| atom.source_id == source_id) {
                atoms[atom.index].element = typographic_superatom_label(label);
                atoms[atom.index].force_label = true;
            }
        }
    }

    let bonds = record
        .bonds
        .iter()
        .map(|bond| RenderBond {
            atom1: bond.atom1_index,
            atom2: bond.atom2_index,
            kind: match bond.order {
                ChemicalBondOrder::Single => BondKind::Single,
                ChemicalBondOrder::Double => BondKind::Double,
                ChemicalBondOrder::Triple => BondKind::Triple,
                ChemicalBondOrder::Aromatic => BondKind::Aromatic,
                ChemicalBondOrder::SingleOrDouble => BondKind::SingleOrDouble,
                ChemicalBondOrder::SingleOrAromatic => BondKind::SingleOrAromatic,
                ChemicalBondOrder::DoubleOrAromatic => BondKind::DoubleOrAromatic,
                ChemicalBondOrder::Any => BondKind::Any,
                ChemicalBondOrder::Coordination => BondKind::Coordination,
                ChemicalBondOrder::Hydrogen => BondKind::Hydrogen,
            },
        })
        .collect();

    RenderMol { atoms, bonds }
}

fn project_record_render(record: &ChemicalRecord, mol: RenderMol) -> RecordRenderProjection {
    let coordinates_preserved = record.atoms.iter().all(|atom| {
        mol.atoms.get(atom.index).is_some_and(|rendered| {
            (rendered.x - atom.coordinates[0]).abs() <= 1e-8
                && (rendered.y - atom.coordinates[1]).abs() <= 1e-8
        })
    });
    let source_atom_indexes = record
        .atoms
        .iter()
        .map(|atom| (atom.source_id, atom.index))
        .collect::<HashMap<_, _>>();
    let contractible_groups = contractible_superatom_indexes(&record.sgroups);
    let mut atom_group_indexes = vec![None; record.atoms.len()];
    for &group_index in &contractible_groups {
        for source_id in &record.sgroups[group_index].atom_source_ids {
            if let Some(&atom_index) = source_atom_indexes.get(source_id) {
                atom_group_indexes[atom_index] = Some(group_index);
            }
        }
    }

    let mut atoms = Vec::with_capacity(record.atoms.len());
    let mut atom_indexes = vec![usize::MAX; record.atoms.len()];
    let mut emitted_groups = HashSet::new();
    let mut contracted_atom_indexes = HashSet::new();
    let mut contracted_sgroup_ids = HashSet::new();

    for original_index in 0..record.atoms.len() {
        let Some(group_index) = atom_group_indexes[original_index] else {
            atom_indexes[original_index] = atoms.len();
            atoms.push(mol.atoms[original_index].clone());
            continue;
        };
        if !emitted_groups.insert(group_index) {
            continue;
        }

        let group = &record.sgroups[group_index];
        let member_indexes = group
            .atom_source_ids
            .iter()
            .filter_map(|source_id| source_atom_indexes.get(source_id).copied())
            .collect::<HashSet<_>>();
        let projected_index = atoms.len();
        for &member_index in &member_indexes {
            atom_indexes[member_index] = projected_index;
            contracted_atom_indexes.insert(member_index);
        }

        let (x, y) = superatom_position(record, &mol, group, &member_indexes);
        let label = group
            .label
            .as_ref()
            .or(group.subscript.as_ref())
            .expect("contractible superatoms have labels");
        atoms.push(RenderAtom {
            element: typographic_superatom_label(label),
            x,
            y,
            hydrogens: 0,
            charge: 0,
            isotope: None,
            radical: None,
            atom_map: None,
            rgroup_label: None,
            query_symbol: false,
            query_negated: false,
            stereo_annotation: None,
            force_label: true,
            hide_label: false,
        });
        contracted_sgroup_ids.insert(group.id);
    }

    debug_assert!(atom_indexes.iter().all(|&index| index != usize::MAX));
    let mut bonds = Vec::with_capacity(record.bonds.len());
    let mut bond_indexes = vec![None; record.bonds.len()];
    for bond in &record.bonds {
        let atom1 = atom_indexes[bond.atom1_index];
        let atom2 = atom_indexes[bond.atom2_index];
        if atom1 == atom2 {
            continue;
        }
        bond_indexes[bond.index] = Some(bonds.len());
        bonds.push(RenderBond {
            atom1,
            atom2,
            kind: mol.bonds[bond.index].kind,
        });
    }

    RecordRenderProjection {
        mol: RenderMol { atoms, bonds },
        atom_indexes,
        bond_indexes,
        contracted_sgroup_ids,
        contracted_atom_indexes,
        coordinates_preserved,
    }
}

fn superatom_position(
    record: &ChemicalRecord,
    mol: &RenderMol,
    group: &crate::ChemicalSGroup,
    member_indexes: &HashSet<usize>,
) -> (f64, f64) {
    let source_atom_indexes = record
        .atoms
        .iter()
        .map(|atom| (atom.source_id, atom.index))
        .collect::<HashMap<_, _>>();
    let explicit_attachment_indexes = group
        .attachment_points
        .iter()
        .filter_map(|point| source_atom_indexes.get(&point.atom_source_id).copied())
        .filter(|index| member_indexes.contains(index))
        .collect::<HashSet<_>>();
    let mut boundary_indexes = Vec::new();
    let mut seen = HashSet::new();
    for bond in &record.bonds {
        let atom1_inside = member_indexes.contains(&bond.atom1_index);
        let atom2_inside = member_indexes.contains(&bond.atom2_index);
        if atom1_inside == atom2_inside {
            continue;
        }
        let inside_index = if atom1_inside {
            bond.atom1_index
        } else {
            bond.atom2_index
        };
        if seen.insert(inside_index) {
            boundary_indexes.push(inside_index);
        }
    }
    let position_indexes = if !explicit_attachment_indexes.is_empty() {
        explicit_attachment_indexes.into_iter().collect::<Vec<_>>()
    } else if boundary_indexes.is_empty() {
        member_indexes.iter().copied().collect::<Vec<_>>()
    } else {
        boundary_indexes
    };
    let count = position_indexes.len() as f64;
    let x = position_indexes
        .iter()
        .map(|&index| mol.atoms[index].x)
        .sum::<f64>()
        / count;
    let y = position_indexes
        .iter()
        .map(|&index| mol.atoms[index].y)
        .sum::<f64>()
        / count;
    (x, y)
}

fn typographic_superatom_label(label: &str) -> String {
    let mut output = String::with_capacity(label.len());
    let mut previous_is_letter = false;
    let mut in_subscript_run = false;
    for character in label.chars() {
        let should_subscript =
            character.is_ascii_digit() && (previous_is_letter || in_subscript_run);
        let rendered = if should_subscript {
            match character {
                '0' => '₀',
                '1' => '₁',
                '2' => '₂',
                '3' => '₃',
                '4' => '₄',
                '5' => '₅',
                '6' => '₆',
                '7' => '₇',
                '8' => '₈',
                '9' => '₉',
                _ => character,
            }
        } else {
            character
        };
        output.push(rendered);
        in_subscript_run = should_subscript;
        previous_is_letter = character.is_alphabetic();
    }
    output
}

fn query_count(value: i8) -> String {
    match value {
        -2 => "*".to_string(),
        -1 => "0".to_string(),
        value => value.to_string(),
    }
}

fn extract_sdf_stereo_annotations(sdf: &str, mol: &sdfrust::Molecule) -> Vec<Option<String>> {
    let mut annotations = vec![Vec::<String>::new(); mol.atoms.len()];

    for (atom_idx, atom) in mol.atoms.iter().enumerate() {
        if let Some(parity @ 1..=3) = atom.stereo_parity {
            annotations[atom_idx].push(format!("CFG={parity}"));
        }
    }

    for group in &mol.stereogroups {
        let label = match group.group_type {
            sdfrust::StereoGroupType::Absolute => "ABS".to_string(),
            sdfrust::StereoGroupType::Or => format!("OR{}", group.group_number),
            sdfrust::StereoGroupType::And => format!("AND{}", group.group_number),
        };
        for &atom_idx in &group.atoms {
            push_stereo_annotation(&mut annotations, atom_idx, &label);
        }
    }

    if mol.format_version == SdfFormat::V3000 {
        let atom_ids = mol
            .atoms
            .iter()
            .enumerate()
            .map(|(atom_idx, atom)| (atom.v3000_id.unwrap_or((atom_idx + 1) as u32), atom_idx))
            .collect::<HashMap<_, _>>();

        for line in v3000_logical_lines(sdf) {
            let Some(label) = v3000_stereo_group_label(&line) else {
                continue;
            };
            let Some(start) = line.find("ATOMS=(") else {
                continue;
            };
            let values = &line[start + "ATOMS=(".len()..];
            let Some(end) = values.find(')') else {
                continue;
            };
            for atom_id in values[..end]
                .split_whitespace()
                .skip(1)
                .filter_map(|value| value.parse::<u32>().ok())
            {
                if let Some(&atom_idx) = atom_ids.get(&atom_id) {
                    push_stereo_annotation(&mut annotations, atom_idx, &label);
                }
            }
        }
    }

    annotations
        .into_iter()
        .map(|annotations| (!annotations.is_empty()).then(|| annotations.join("; ")))
        .collect()
}

fn push_stereo_annotation(annotations: &mut [Vec<String>], atom_idx: usize, label: &str) {
    let Some(atom_annotations) = annotations.get_mut(atom_idx) else {
        return;
    };
    if !atom_annotations
        .iter()
        .any(|annotation| annotation == label)
    {
        atom_annotations.push(label.to_string());
    }
}

fn v3000_stereo_group_label(line: &str) -> Option<String> {
    let tag = line.split_whitespace().next()?.strip_prefix("MDLV30/")?;
    let label = if tag == "STEABS" {
        "ABS".to_string()
    } else if let Some(group) = tag.strip_prefix("STEREL") {
        format!("OR{}", if group.is_empty() { "1" } else { group })
    } else if let Some(group) = tag.strip_prefix("STERAC") {
        format!("AND{}", if group.is_empty() { "1" } else { group })
    } else {
        return None;
    };
    Some(label)
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

fn render_mol_from_smiles(
    graph: &SmilesGraph,
    coords: &[(f32, f32)],
    depicted_centers: &HashSet<usize>,
) -> RenderMol {
    let atoms = graph
        .atoms
        .iter()
        .zip(coords.iter())
        .enumerate()
        .map(|(atom_idx, (atom, &(x, y)))| RenderAtom {
            element: atom.element.clone(),
            x: x as f64,
            y: y as f64,
            hydrogens: atom.hydrogens,
            charge: atom.charge,
            isotope: atom.isotope,
            radical: None,
            atom_map: atom.atom_map,
            rgroup_label: None,
            query_symbol: false,
            query_negated: false,
            stereo_annotation: (!depicted_centers.contains(&atom_idx))
                .then(|| atom.stereo_annotation.clone())
                .flatten(),
            force_label: false,
            hide_label: false,
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
                SmilesBondKind::Quadruple => BondKind::Quadruple,
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
        for hydrogen_offset in 0..hydrogen_count {
            let hydrogen_idx = expanded.atoms.len();
            if let Some(stereo) = &mut expanded.atoms[atom_idx].stereo {
                if stereo.looking_from.is_none() {
                    stereo.looking_from = Some(hydrogen_idx);
                } else if stereo.atom1.is_none() {
                    stereo.atom1 = Some(hydrogen_idx);
                } else if stereo.atom2.is_none() {
                    stereo.atom2 = Some(hydrogen_idx);
                }
            }
            expanded.atoms.push(SmilesAtom {
                element: "H".to_string(),
                atomic_number: 1,
                hydrogens: 0,
                charge: 0,
                isotope: None,
                atom_map: None,
                stereo: None,
                extended_stereo: None,
                stereo_annotation: None,
            });
            expanded.bonds.push(SmilesBond {
                atom1: atom_idx,
                atom2: hydrogen_idx,
                kind: SmilesBondKind::Single,
            });
            // A bracket hydrogen follows the incoming neighbor, when present,
            // and precedes ring bonds, branches, and the continuing chain.
            let insertion_position = graph.hydrogen_order_positions[atom_idx]
                .saturating_add(usize::from(hydrogen_offset))
                .min(expanded.neighbor_order[atom_idx].len());
            expanded.neighbor_order[atom_idx].insert(insertion_position, hydrogen_idx);
            expanded.neighbor_order.push(vec![atom_idx]);
            expanded.hydrogen_order_positions.push(1);
        }
    }

    expanded
}

fn parse_smiles_graph(smiles: &str) -> Result<SmilesGraph, String> {
    let molecule = parse_smiles(smiles).map_err(|e| e.to_string())?;
    let aromatic_double_bonds = crate::ast::aromaticity::aromatic_kekule_bonds(&molecule)
        .map_err(|error| error.to_string())?;
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
            SmilesBondType::Quadruple => SmilesBondKind::Quadruple,
            SmilesBondType::Aromatic => {
                if aromatic_double_bonds.contains(&bond_idx) {
                    SmilesBondKind::Double
                } else {
                    SmilesBondKind::Single
                }
            }
            SmilesBondType::Disconnected => continue,
        };

        original_to_graph[bond_idx] = Some(bonds.len());
        bonds.push(SmilesBond {
            atom1: source,
            atom2: target,
            kind,
        });
    }

    let (chiral_neighbor_order, hydrogen_order_positions) =
        smiles_neighbor_order(smiles, molecule.nodes().len()).unwrap_or_else(|| {
            (
                incident_bonds
                    .iter()
                    .map(|bonds| bonds.iter().map(|&(_, neighbor)| neighbor).collect())
                    .collect(),
                vec![0; molecule.nodes().len()],
            )
        });

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
                stereo: atom_stereo_spec(
                    atom_idx,
                    node.chirality(),
                    node.hydrogens(),
                    &chiral_neighbor_order,
                    hydrogen_order_positions
                        .get(atom_idx)
                        .copied()
                        .unwrap_or_default(),
                ),
                extended_stereo: extended_atom_stereo_spec(
                    node.chirality(),
                    chiral_neighbor_order.get(atom_idx).map(Vec::as_slice),
                ),
                stereo_annotation: chirality_annotation(node.chirality()),
            }
        })
        .collect::<Vec<_>>();

    let double_bond_stereo = build_double_bond_stereo_specs(&molecule, &original_to_graph);

    Ok(SmilesGraph {
        atoms,
        bonds,
        neighbor_order: chiral_neighbor_order,
        hydrogen_order_positions,
        double_bond_stereo,
    })
}

fn smiles_neighbor_order(smiles: &str, atom_count: usize) -> Option<(Vec<Vec<usize>>, Vec<usize>)> {
    let bytes = smiles.as_bytes();
    let mut neighbors = vec![Vec::<Option<usize>>::new(); atom_count];
    let mut hydrogen_order_positions = vec![0; atom_count];
    let mut current: Option<usize> = None;
    let mut branches = Vec::new();
    let mut rings = HashMap::<u16, (usize, usize)>::new();
    let mut next_atom = 0usize;
    let mut cursor = 0usize;

    while cursor < bytes.len() {
        let byte = bytes[cursor];
        let atom_width = match byte {
            b'[' => {
                let end = bytes[cursor + 1..].iter().position(|&byte| byte == b']')?;
                end + 2
            }
            b'*' => 1,
            b'C' if bytes.get(cursor + 1) == Some(&b'l') => 2,
            b'B' if bytes.get(cursor + 1) == Some(&b'r') => 2,
            byte if byte.is_ascii_alphabetic() => 1,
            _ => 0,
        };

        if atom_width > 0 {
            if next_atom >= atom_count {
                return None;
            }
            let atom = next_atom;
            next_atom += 1;
            if let Some(previous) = current {
                neighbors[previous].push(Some(atom));
                neighbors[atom].push(Some(previous));
            }
            // Bracket properties occur after the incoming connection and
            // before any ring bonds, branches, or continuing-chain neighbor.
            hydrogen_order_positions[atom] = neighbors[atom].len();
            current = Some(atom);
            cursor += atom_width;
            continue;
        }

        match byte {
            b'(' => branches.push(current),
            b')' => current = branches.pop()?,
            b'.' => current = None,
            b'%' => {
                let first = *bytes.get(cursor + 1)?;
                let second = *bytes.get(cursor + 2)?;
                if !first.is_ascii_digit() || !second.is_ascii_digit() {
                    return None;
                }
                let ring = u16::from(first - b'0') * 10 + u16::from(second - b'0');
                record_ring_neighbor(&mut neighbors, &mut rings, current?, ring)?;
                cursor += 3;
                continue;
            }
            byte if byte.is_ascii_digit() => {
                record_ring_neighbor(&mut neighbors, &mut rings, current?, u16::from(byte - b'0'))?;
            }
            b' ' | b'\t' | b'\r' | b'\n' => break,
            _ => {}
        }
        cursor += 1;
    }

    if next_atom != atom_count || !branches.is_empty() || !rings.is_empty() {
        return None;
    }
    let neighbors = neighbors
        .into_iter()
        .map(|neighbors| neighbors.into_iter().collect::<Option<Vec<_>>>())
        .collect::<Option<Vec<_>>>()?;
    Some((neighbors, hydrogen_order_positions))
}

fn record_ring_neighbor(
    neighbors: &mut [Vec<Option<usize>>],
    rings: &mut HashMap<u16, (usize, usize)>,
    atom: usize,
    ring: u16,
) -> Option<()> {
    if let Some((opening_atom, opening_position)) = rings.remove(&ring) {
        *neighbors.get_mut(opening_atom)?.get_mut(opening_position)? = Some(atom);
        neighbors.get_mut(atom)?.push(Some(opening_atom));
    } else {
        let position = neighbors.get(atom)?.len();
        neighbors.get_mut(atom)?.push(None);
        rings.insert(ring, (atom, position));
    }
    Some(())
}

fn atom_stereo_spec(
    atom_idx: usize,
    chirality: Option<SmilesChirality>,
    hydrogen_count: u8,
    neighbor_order: &[Vec<usize>],
    hydrogen_order_position: usize,
) -> Option<AtomStereoSpec> {
    let direction = match chirality {
        Some(SmilesChirality::TH1) => AtomStereoDirection::CounterClockwise,
        Some(SmilesChirality::TH2) => AtomStereoDirection::Clockwise,
        _ => return None,
    };

    let mut neighbors = neighbor_order
        .get(atom_idx)?
        .iter()
        .map(|&neighbor| Some(neighbor))
        .collect::<Vec<_>>();
    match hydrogen_count {
        0 => {}
        1 => {
            // The bracket H occupies its lexical position after the incoming
            // neighbor (if any), rather than preceding every graph neighbor.
            if hydrogen_order_position > neighbors.len() {
                return None;
            }
            neighbors.insert(hydrogen_order_position, None);
        }
        _ => return None,
    }

    if neighbors.len() == 3 {
        // A three-coordinate tetrahedral center may use an implicit lone pair
        // as its fourth neighbor. Coordgen represents that position as null.
        neighbors.push(None);
    }
    if neighbors.len() != 4 {
        return None;
    }

    Some(AtomStereoSpec {
        looking_from: neighbors[0],
        atom1: neighbors[1],
        atom2: neighbors[2],
        direction,
    })
}

fn extended_atom_stereo_spec(
    chirality: Option<SmilesChirality>,
    neighbor_order: Option<&[usize]>,
) -> Option<ExtendedAtomStereoSpec> {
    let chirality = chirality?;
    match chirality {
        SmilesChirality::AL1 => Some(ExtendedAtomStereoSpec::Allene {
            direction: AtomStereoDirection::CounterClockwise,
        }),
        SmilesChirality::AL2 => Some(ExtendedAtomStereoSpec::Allene {
            direction: AtomStereoDirection::Clockwise,
        }),
        SmilesChirality::SP1 | SmilesChirality::SP2 | SmilesChirality::SP3 => {
            let ligands = neighbor_order?.try_into().ok()?;
            let shape = match chirality {
                SmilesChirality::SP1 => PlanarShape::U,
                SmilesChirality::SP2 => PlanarShape::Four,
                SmilesChirality::SP3 => PlanarShape::Z,
                _ => unreachable!(),
            };
            Some(ExtendedAtomStereoSpec::SquarePlanar { ligands, shape })
        }
        chirality
            if (SmilesChirality::TB1 as u8..=SmilesChirality::TB20 as u8)
                .contains(&(chirality as u8)) =>
        {
            let neighbors = neighbor_order?;
            let number = usize::from(chirality as u8 - SmilesChirality::TB1 as u8);
            let (from_position, to_position, direction) = TB_PERMUTATIONS[number];
            let equatorial = neighbors
                .iter()
                .enumerate()
                .filter(|(position, _)| *position != from_position && *position != to_position)
                .map(|(_, &neighbor)| neighbor)
                .collect::<Vec<_>>()
                .try_into()
                .ok()?;
            Some(ExtendedAtomStereoSpec::TrigonalBipyramidal {
                axis_from: *neighbors.get(from_position)?,
                axis_to: *neighbors.get(to_position)?,
                equatorial,
                direction,
            })
        }
        chirality
            if (SmilesChirality::OH1 as u8..=SmilesChirality::OH30 as u8)
                .contains(&(chirality as u8)) =>
        {
            let neighbors = neighbor_order?;
            let number = usize::from(chirality as u8 - SmilesChirality::OH1 as u8);
            let (to_position, shape, direction) = OH_PERMUTATIONS[number];
            let equatorial = neighbors
                .iter()
                .enumerate()
                .filter(|(position, _)| *position != 0 && *position != to_position)
                .map(|(_, &neighbor)| neighbor)
                .collect::<Vec<_>>()
                .try_into()
                .ok()?;
            Some(ExtendedAtomStereoSpec::Octahedral {
                axis_from: *neighbors.first()?,
                axis_to: *neighbors.get(to_position)?,
                equatorial,
                shape,
                direction,
            })
        }
        _ => None,
    }
}

const TB_PERMUTATIONS: [(usize, usize, AtomStereoDirection); 20] = [
    (0, 4, AtomStereoDirection::CounterClockwise),
    (0, 4, AtomStereoDirection::Clockwise),
    (0, 3, AtomStereoDirection::CounterClockwise),
    (0, 3, AtomStereoDirection::Clockwise),
    (0, 2, AtomStereoDirection::CounterClockwise),
    (0, 2, AtomStereoDirection::Clockwise),
    (0, 1, AtomStereoDirection::CounterClockwise),
    (0, 1, AtomStereoDirection::Clockwise),
    (1, 4, AtomStereoDirection::CounterClockwise),
    (1, 3, AtomStereoDirection::CounterClockwise),
    (1, 4, AtomStereoDirection::Clockwise),
    (1, 3, AtomStereoDirection::Clockwise),
    (1, 2, AtomStereoDirection::CounterClockwise),
    (1, 2, AtomStereoDirection::Clockwise),
    (2, 4, AtomStereoDirection::CounterClockwise),
    (2, 3, AtomStereoDirection::CounterClockwise),
    (3, 4, AtomStereoDirection::CounterClockwise),
    (3, 4, AtomStereoDirection::Clockwise),
    (2, 3, AtomStereoDirection::Clockwise),
    (2, 4, AtomStereoDirection::Clockwise),
];

const OH_PERMUTATIONS: [(usize, PlanarShape, AtomStereoDirection); 30] = [
    (5, PlanarShape::U, AtomStereoDirection::CounterClockwise),
    (5, PlanarShape::U, AtomStereoDirection::Clockwise),
    (4, PlanarShape::U, AtomStereoDirection::CounterClockwise),
    (5, PlanarShape::Z, AtomStereoDirection::CounterClockwise),
    (4, PlanarShape::Z, AtomStereoDirection::CounterClockwise),
    (3, PlanarShape::U, AtomStereoDirection::CounterClockwise),
    (3, PlanarShape::Z, AtomStereoDirection::CounterClockwise),
    (5, PlanarShape::Four, AtomStereoDirection::Clockwise),
    (4, PlanarShape::Four, AtomStereoDirection::Clockwise),
    (5, PlanarShape::Four, AtomStereoDirection::CounterClockwise),
    (4, PlanarShape::Four, AtomStereoDirection::CounterClockwise),
    (3, PlanarShape::Four, AtomStereoDirection::Clockwise),
    (3, PlanarShape::Four, AtomStereoDirection::CounterClockwise),
    (5, PlanarShape::Z, AtomStereoDirection::Clockwise),
    (4, PlanarShape::Z, AtomStereoDirection::Clockwise),
    (4, PlanarShape::U, AtomStereoDirection::Clockwise),
    (3, PlanarShape::Z, AtomStereoDirection::Clockwise),
    (3, PlanarShape::U, AtomStereoDirection::Clockwise),
    (2, PlanarShape::U, AtomStereoDirection::CounterClockwise),
    (2, PlanarShape::Z, AtomStereoDirection::CounterClockwise),
    (2, PlanarShape::Four, AtomStereoDirection::Clockwise),
    (2, PlanarShape::Four, AtomStereoDirection::CounterClockwise),
    (2, PlanarShape::Z, AtomStereoDirection::Clockwise),
    (2, PlanarShape::U, AtomStereoDirection::Clockwise),
    (1, PlanarShape::U, AtomStereoDirection::CounterClockwise),
    (1, PlanarShape::Z, AtomStereoDirection::CounterClockwise),
    (1, PlanarShape::Four, AtomStereoDirection::Clockwise),
    (1, PlanarShape::Four, AtomStereoDirection::CounterClockwise),
    (1, PlanarShape::Z, AtomStereoDirection::Clockwise),
    (1, PlanarShape::U, AtomStereoDirection::Clockwise),
];

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
            SmilesBondKind::Quadruple => 4,
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
    for atom in 0..atom_count {
        let x = read_f32_le(coords_data, &mut offset)?;
        let y = read_f32_le(coords_data, &mut offset)?;
        if !x.is_finite() || !y.is_finite() {
            return Err(format!(
                "Coordinate payload atom {} has non-finite coordinates",
                atom + 1
            ));
        }
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

fn apply_extended_stereo_depictions(
    graph: &SmilesGraph,
    coords: &mut [(f32, f32)],
    stereo_map: &mut HashMap<(usize, usize), (u8, bool)>,
) -> HashSet<usize> {
    let mut depicted = HashSet::new();

    for center in 0..graph.atoms.len() {
        let Some(spec) = graph.atoms[center].extended_stereo.clone() else {
            continue;
        };

        let was_depicted = match spec {
            ExtendedAtomStereoSpec::Allene { direction } => {
                depict_allene_stereo(graph, center, direction, stereo_map)
            }
            ExtendedAtomStereoSpec::SquarePlanar { ligands, shape } => {
                let targets = planar_target_angles(shape, None);
                layout_coordination_center(graph, coords, center, &ligands, &targets)
            }
            ExtendedAtomStereoSpec::TrigonalBipyramidal {
                axis_from,
                axis_to,
                equatorial,
                direction,
            } => {
                let equatorial_targets = match direction {
                    AtomStereoDirection::CounterClockwise => [90.0, 210.0, 330.0],
                    AtomStereoDirection::Clockwise => [90.0, -30.0, -150.0],
                };
                let ligands = [
                    axis_from,
                    axis_to,
                    equatorial[0],
                    equatorial[1],
                    equatorial[2],
                ];
                let targets = [
                    0.0,
                    180.0,
                    equatorial_targets[0],
                    equatorial_targets[1],
                    equatorial_targets[2],
                ];
                if has_single_bond(graph, center, axis_from)
                    && has_single_bond(graph, center, axis_to)
                    && layout_coordination_center(graph, coords, center, &ligands, &targets)
                {
                    add_axis_stereo(stereo_map, center, axis_from, axis_to);
                    true
                } else {
                    false
                }
            }
            ExtendedAtomStereoSpec::Octahedral {
                axis_from,
                axis_to,
                equatorial,
                shape,
                direction,
            } => {
                let planar_targets = planar_target_angles(shape, Some(direction));
                let ligands = [
                    axis_from,
                    axis_to,
                    equatorial[0],
                    equatorial[1],
                    equatorial[2],
                    equatorial[3],
                ];
                let targets = [
                    45.0,
                    225.0,
                    planar_targets[0],
                    planar_targets[1],
                    planar_targets[2],
                    planar_targets[3],
                ];
                if has_single_bond(graph, center, axis_from)
                    && has_single_bond(graph, center, axis_to)
                    && layout_coordination_center(graph, coords, center, &ligands, &targets)
                {
                    add_axis_stereo(stereo_map, center, axis_from, axis_to);
                    true
                } else {
                    false
                }
            }
        };

        if was_depicted {
            depicted.insert(center);
        }
    }

    depicted
}

fn planar_target_angles(shape: PlanarShape, direction: Option<AtomStereoDirection>) -> [f64; 4] {
    // The slots progress anticlockwise. The permutations are the three
    // Hamiltonian path shapes defined by OpenSMILES: U, 4, and Z.
    let slots = match shape {
        PlanarShape::U => [0, 1, 2, 3],
        PlanarShape::Four => [0, 2, 3, 1],
        PlanarShape::Z => [0, 1, 3, 2],
    };
    let slot_angles = [90.0, 180.0, 270.0, 0.0];
    let mut targets = slots.map(|slot| slot_angles[slot]);
    if direction == Some(AtomStereoDirection::Clockwise) {
        targets = targets.map(|target| 180.0 - target);
    }
    targets
}

fn layout_coordination_center<const N: usize>(
    graph: &SmilesGraph,
    coords: &mut [(f32, f32)],
    center: usize,
    ligands: &[usize; N],
    target_angles: &[f64; N],
) -> bool {
    if center >= coords.len() || ligands.iter().any(|&ligand| ligand >= coords.len()) {
        return false;
    }

    let Some(components) = ligand_branch_components(graph, center, ligands) else {
        return false;
    };
    let (center_x, center_y) = coords[center];
    let rotations = ligands
        .iter()
        .zip(target_angles.iter())
        .map(|(&ligand, &target_angle)| {
            let (ligand_x, ligand_y) = coords[ligand];
            let dx = ligand_x - center_x;
            let dy = ligand_y - center_y;
            if dx.abs() < f32::EPSILON && dy.abs() < f32::EPSILON {
                return None;
            }
            let current_angle = f64::from(dy).atan2(f64::from(dx));
            current_angle
                .is_finite()
                .then_some(target_angle.to_radians() - current_angle)
        })
        .collect::<Option<Vec<_>>>();
    let Some(rotations) = rotations else {
        return false;
    };

    for (component, rotation) in components.iter().zip(rotations) {
        let (sin_rotation, cos_rotation) = rotation.sin_cos();
        for &atom in component {
            let dx = f64::from(coords[atom].0 - center_x);
            let dy = f64::from(coords[atom].1 - center_y);
            coords[atom] = (
                center_x + (dx * cos_rotation - dy * sin_rotation) as f32,
                center_y + (dx * sin_rotation + dy * cos_rotation) as f32,
            );
        }
    }

    true
}

fn ligand_branch_components<const N: usize>(
    graph: &SmilesGraph,
    center: usize,
    ligands: &[usize; N],
) -> Option<Vec<Vec<usize>>> {
    let mut adjacency = vec![Vec::new(); graph.atoms.len()];
    for bond in &graph.bonds {
        adjacency.get_mut(bond.atom1)?.push(bond.atom2);
        adjacency.get_mut(bond.atom2)?.push(bond.atom1);
    }

    let mut claimed = HashSet::new();
    let mut components = Vec::with_capacity(N);
    for &ligand in ligands {
        if !adjacency.get(center)?.contains(&ligand) || claimed.contains(&ligand) {
            return None;
        }

        let mut component = Vec::new();
        let mut queue = VecDeque::from([ligand]);
        let mut visited = HashSet::from([center]);
        while let Some(atom) = queue.pop_front() {
            if !visited.insert(atom) {
                continue;
            }
            component.push(atom);
            for &neighbor in adjacency.get(atom)? {
                if !visited.contains(&neighbor) {
                    queue.push_back(neighbor);
                }
            }
        }
        if component.iter().any(|atom| claimed.contains(atom)) {
            return None;
        }
        claimed.extend(component.iter().copied());
        components.push(component);
    }

    Some(components)
}

fn add_axis_stereo(
    stereo_map: &mut HashMap<(usize, usize), (u8, bool)>,
    center: usize,
    axis_from: usize,
    axis_to: usize,
) {
    insert_stereo_edge(stereo_map, center, axis_from, 1);
    insert_stereo_edge(stereo_map, center, axis_to, 6);
}

fn depict_allene_stereo(
    graph: &SmilesGraph,
    center: usize,
    direction: AtomStereoDirection,
    stereo_map: &mut HashMap<(usize, usize), (u8, bool)>,
) -> bool {
    let double_neighbors = graph
        .neighbor_order
        .get(center)
        .into_iter()
        .flatten()
        .copied()
        .filter(|&neighbor| has_bond_kind(graph, center, neighbor, SmilesBondKind::Double))
        .collect::<Vec<_>>();
    let Ok([first_side, second_side]) = <[usize; 2]>::try_from(double_neighbors) else {
        return false;
    };
    let Some((first_terminal, first_previous, first_length)) =
        allene_terminal(graph, center, first_side)
    else {
        return false;
    };
    let Some((second_terminal, second_previous, second_length)) =
        allene_terminal(graph, center, second_side)
    else {
        return false;
    };
    if (first_length + second_length) % 2 != 0 || first_terminal == second_terminal {
        return false;
    }

    let first_substituents = allene_terminal_substituents(graph, first_terminal, first_previous);
    let second_substituents = allene_terminal_substituents(graph, second_terminal, second_previous);
    if first_substituents.is_empty()
        || second_substituents.is_empty()
        || first_substituents.len() + usize::from(graph.atoms[first_terminal].hydrogens) < 2
        || second_substituents.len() + usize::from(graph.atoms[second_terminal].hydrogens) < 2
    {
        return false;
    }

    // OpenSMILES depicts the second terminal plane in perspective. For AL1
    // its first substituent points away and the second points towards the
    // viewer; AL2 reverses those two bonds.
    let (first_style, second_style) = match direction {
        AtomStereoDirection::CounterClockwise => (6, 1),
        AtomStereoDirection::Clockwise => (1, 6),
    };
    let styles = [first_style, second_style];
    let implicit_hydrogen_offset = usize::from(graph.atoms[second_terminal].hydrogens.min(1));
    for (position, &substituent) in second_substituents.iter().enumerate() {
        let Some(&style) = styles.get(implicit_hydrogen_offset + position) else {
            break;
        };
        insert_stereo_edge(stereo_map, second_terminal, substituent, style);
    }
    true
}

fn allene_terminal(
    graph: &SmilesGraph,
    center: usize,
    first: usize,
) -> Option<(usize, usize, usize)> {
    let mut previous = center;
    let mut current = first;
    let mut length = 1usize;
    let mut visited = HashSet::from([center]);

    loop {
        if !visited.insert(current) {
            return None;
        }
        let continuations = graph
            .neighbor_order
            .get(current)?
            .iter()
            .copied()
            .filter(|&neighbor| {
                neighbor != previous
                    && has_bond_kind(graph, current, neighbor, SmilesBondKind::Double)
            })
            .collect::<Vec<_>>();
        match continuations.as_slice() {
            [] => return Some((current, previous, length)),
            &[next] => {
                previous = current;
                current = next;
                length += 1;
            }
            _ => return None,
        }
    }
}

fn allene_terminal_substituents(
    graph: &SmilesGraph,
    terminal: usize,
    previous: usize,
) -> Vec<usize> {
    graph
        .neighbor_order
        .get(terminal)
        .into_iter()
        .flatten()
        .copied()
        .filter(|&neighbor| {
            neighbor != previous && has_bond_kind(graph, terminal, neighbor, SmilesBondKind::Single)
        })
        .collect()
}

fn has_single_bond(graph: &SmilesGraph, atom1: usize, atom2: usize) -> bool {
    has_bond_kind(graph, atom1, atom2, SmilesBondKind::Single)
}

fn has_bond_kind(graph: &SmilesGraph, atom1: usize, atom2: usize, kind: SmilesBondKind) -> bool {
    graph.bonds.iter().any(|bond| {
        ((bond.atom1 == atom1 && bond.atom2 == atom2)
            || (bond.atom1 == atom2 && bond.atom2 == atom1))
            && bond.kind == kind
    })
}

fn ast_from_chemical_record(
    record: &ChemicalRecord,
    projection: &RecordRenderProjection,
    mode_str: &str,
    stereo_map: &HashMap<(usize, usize), (u8, bool)>,
) -> Vec<Command> {
    let stereo_map = project_stereo_map(stereo_map, projection);
    let mut commands = ast_from_render_mol(&projection.mol, mode_str, &stereo_map);
    let depictions = ctfile_depictions(record, projection);
    if !depictions.sgroups.is_empty()
        || !depictions.highlights.is_empty()
        || !depictions.atom_queries.is_empty()
        || !depictions.bond_queries.is_empty()
        || !depictions.atom_annotations.is_empty()
        || !depictions.variable_attachments.is_empty()
    {
        commands.push(Command::Ctfile {
            sgroups: depictions.sgroups,
            highlights: depictions.highlights,
            atom_queries: depictions.atom_queries,
            bond_queries: depictions.bond_queries,
            atom_annotations: depictions.atom_annotations,
            variable_attachments: depictions.variable_attachments,
        });
    }
    commands
}

fn project_stereo_map(
    stereo_map: &HashMap<(usize, usize), (u8, bool)>,
    projection: &RecordRenderProjection,
) -> HashMap<(usize, usize), (u8, bool)> {
    stereo_map
        .iter()
        .filter_map(|(&(atom1, atom2), &stereo)| {
            let atom1 = *projection.atom_indexes.get(atom1)?;
            let atom2 = *projection.atom_indexes.get(atom2)?;
            (atom1 != atom2).then_some(((atom1, atom2), stereo))
        })
        .collect()
}

struct CtfileDepictions {
    sgroups: Vec<CtfileSGroupDepiction>,
    highlights: Vec<CtfileHighlightDepiction>,
    atom_queries: Vec<CtfileAtomQueryDepiction>,
    bond_queries: Vec<CtfileBondQueryDepiction>,
    atom_annotations: Vec<CtfileAtomAnnotationDepiction>,
    variable_attachments: Vec<CtfileVariableAttachmentDepiction>,
}

fn ctfile_depictions(
    record: &ChemicalRecord,
    projection: &RecordRenderProjection,
) -> CtfileDepictions {
    let mol = &projection.mol;
    let atom_indexes = record
        .atoms
        .iter()
        .map(|atom| (atom.source_id, projection.atom_indexes[atom.index]))
        .collect::<HashMap<_, _>>();
    let bond_indexes = record
        .bonds
        .iter()
        .filter_map(|bond| projection.bond_indexes[bond.index].map(|index| (bond.source_id, index)))
        .collect::<HashMap<_, _>>();
    let average_length = render_average_bond_length(mol);

    let sgroups = record
        .sgroups
        .iter()
        .filter(|group| group.kind != "unknown")
        .filter(|group| !(group.kind == "superatom" && group.atom_source_ids.len() == 1))
        .filter(|group| !projection.contracted_sgroup_ids.contains(&group.id))
        .filter_map(|group| {
            let mut seen = HashSet::new();
            let indexes = group
                .atom_source_ids
                .iter()
                .filter_map(|id| atom_indexes.get(id).copied())
                .filter(|index| seen.insert(*index))
                .collect::<Vec<_>>();
            if indexes.is_empty() {
                return None;
            }

            let (left, right) = if projection.coordinates_preserved && group.brackets.len() >= 2 {
                let mut brackets = group.brackets.clone();
                brackets.sort_by(|left, right| {
                    ((left[0] + left[2]) / 2.0).total_cmp(&((right[0] + right[2]) / 2.0))
                });
                (brackets[0], *brackets.last().expect("two brackets exist"))
            } else {
                let min_x = indexes
                    .iter()
                    .map(|&index| mol.atoms[index].x)
                    .fold(f64::INFINITY, f64::min);
                let max_x = indexes
                    .iter()
                    .map(|&index| mol.atoms[index].x)
                    .fold(f64::NEG_INFINITY, f64::max);
                let min_y = indexes
                    .iter()
                    .map(|&index| mol.atoms[index].y)
                    .fold(f64::INFINITY, f64::min);
                let max_y = indexes
                    .iter()
                    .map(|&index| mol.atoms[index].y)
                    .fold(f64::NEG_INFINITY, f64::max);
                let padding = average_length * 0.35;
                (
                    [
                        min_x - padding,
                        max_y + padding,
                        min_x - padding,
                        min_y - padding,
                    ],
                    [
                        max_x + padding,
                        max_y + padding,
                        max_x + padding,
                        min_y - padding,
                    ],
                )
            };
            let label = sgroup_depiction_label(group);
            Some(CtfileSGroupDepiction {
                id: group.id,
                kind: group.kind.clone(),
                left_top: ctfile_point(left[0], left[1], &indexes, mol, average_length),
                left_bottom: ctfile_point(left[2], left[3], &indexes, mol, average_length),
                right_top: ctfile_point(right[0], right[1], &indexes, mol, average_length),
                right_bottom: ctfile_point(right[2], right[3], &indexes, mol, average_length),
                label,
            })
        })
        .collect();

    let highlights = record
        .collections
        .iter()
        .filter(|collection| collection.kind == "highlight")
        .map(|collection| {
            let mut seen_atoms = HashSet::new();
            let mut highlighted_atoms = collection
                .atom_source_ids
                .iter()
                .filter_map(|id| atom_indexes.get(id).copied())
                .filter(|index| seen_atoms.insert(*index))
                .collect::<Vec<_>>();
            let source_bond_indexes = record
                .bonds
                .iter()
                .map(|bond| (bond.source_id, bond.index))
                .collect::<HashMap<_, _>>();
            let mut bonds = Vec::new();
            for source_id in &collection.bond_source_ids {
                let Some(&source_bond_index) = source_bond_indexes.get(source_id) else {
                    continue;
                };
                if let Some(&bond_index) = bond_indexes.get(source_id) {
                    let bond = &mol.bonds[bond_index];
                    bonds.push(CtfileHighlightBond {
                        bond_index,
                        atom1_index: bond.atom1,
                        atom2_index: bond.atom2,
                    });
                } else {
                    let source_bond = &record.bonds[source_bond_index];
                    let atom1_index = projection.atom_indexes[source_bond.atom1_index];
                    let atom2_index = projection.atom_indexes[source_bond.atom2_index];
                    if atom1_index == atom2_index && seen_atoms.insert(atom1_index) {
                        highlighted_atoms.push(atom1_index);
                    }
                }
            }
            CtfileHighlightDepiction {
                atom_indexes: highlighted_atoms,
                bonds,
            }
        })
        .filter(|highlight| !highlight.atom_indexes.is_empty() || !highlight.bonds.is_empty())
        .collect();

    let bond_queries = record
        .bonds
        .iter()
        .filter_map(|bond| {
            let bond_index = projection.bond_indexes[bond.index]?;
            bond.topology.map(|topology| CtfileBondQueryDepiction {
                bond_index,
                label: match topology {
                    1 => "rn",
                    2 => "ch",
                    _ => "rn/ch",
                }
                .to_string(),
            })
        })
        .collect();

    let atom_annotations = record
        .link_nodes
        .iter()
        .filter(|node| node.min_repeat == 1 && link_node_atoms_are_renderable(node, &atom_indexes))
        .filter_map(|node| {
            let source_id = node.connections.first()?.atom_source_id;
            let atom_index = *atom_indexes.get(&source_id)?;
            Some(CtfileAtomAnnotationDepiction {
                atom_index,
                kind: "link-node".to_string(),
                label: if node.min_repeat == node.max_repeat {
                    node.min_repeat.to_string()
                } else {
                    format!("{}–{}", node.min_repeat, node.max_repeat)
                },
            })
        })
        .collect::<Vec<_>>();

    let variable_attachments = record
        .bonds
        .iter()
        .filter(|bond| {
            bond.endpoint_source_ids.len() >= 2
                && matches!(bond.attachment_mode.as_deref(), Some("ALL" | "ANY"))
                && bond
                    .endpoint_source_ids
                    .iter()
                    .all(|source_id| atom_indexes.contains_key(source_id))
        })
        .filter_map(|bond| {
            let bond_index = projection.bond_indexes[bond.index]?;
            let mut seen = HashSet::new();
            let endpoint_atom_indexes = bond
                .endpoint_source_ids
                .iter()
                .filter_map(|source_id| atom_indexes.get(source_id).copied())
                .filter(|index| seen.insert(*index))
                .collect::<Vec<_>>();
            (!endpoint_atom_indexes.is_empty()).then_some(CtfileVariableAttachmentDepiction {
                bond_index,
                endpoint_atom_indexes,
                mode: bond
                    .attachment_mode
                    .clone()
                    .expect("validated attachment mode"),
            })
        })
        .collect();

    let atom_queries = record
        .atoms
        .iter()
        .filter(|atom| !projection.contracted_atom_indexes.contains(&atom.index))
        .filter_map(|atom| {
            let query = atom.query.as_ref()?;
            let label = atom_query_constraint_label(query);
            let compact_label = atom_query_compact_label(query);
            let constraint_lines = atom_query_constraint_lines(query);
            let symbol = atom_query_symbol(query);
            (symbol.is_some() || !label.is_empty()).then_some(CtfileAtomQueryDepiction {
                atom_index: projection.atom_indexes[atom.index],
                symbol,
                label,
                compact_label,
                constraint_lines,
            })
        })
        .collect();

    CtfileDepictions {
        sgroups,
        highlights,
        atom_queries,
        bond_queries,
        atom_annotations,
        variable_attachments,
    }
}

fn link_node_atoms_are_renderable(
    node: &crate::ChemicalLinkNode,
    atom_indexes: &HashMap<u32, usize>,
) -> bool {
    node.max_repeat > 0
        && !node.connections.is_empty()
        && node.connections.iter().all(|connection| {
            atom_indexes.contains_key(&connection.atom_source_id)
                && atom_indexes.contains_key(&connection.neighbor_source_id)
        })
}

fn atom_query_constraint_label(query: &crate::ChemicalAtomQuery) -> String {
    atom_query_constraint_terms(query).join(" · ")
}

fn atom_query_compact_label(query: &crate::ChemicalAtomQuery) -> String {
    let mut constraints = Vec::new();
    if let Some(value) = query.hydrogen_count {
        constraints.push(format!("(H{})", query_count(value)));
    }
    if let Some(value) = query.substitution_count {
        constraints.push(format!("(s{})", query_count(value)));
    }
    if query.unsaturated == Some(true) {
        constraints.push("(u)".to_string());
    }
    if let Some(value) = query.ring_bond_count {
        constraints.push(format!("(r{})", query_count(value)));
    }
    if let Some(value) = query.valence {
        constraints.push(format!("(v{})", query_count(value)));
    }
    constraints.join(" ")
}

fn atom_query_constraint_lines(query: &crate::ChemicalAtomQuery) -> Vec<String> {
    atom_query_constraint_terms(query)
        .chunks(2)
        .map(|line| line.join(" · "))
        .collect()
}

fn atom_query_constraint_terms(query: &crate::ChemicalAtomQuery) -> Vec<String> {
    let mut constraints = Vec::new();
    if let Some(value) = query.hydrogen_count {
        constraints.push(match value {
            0 | -1 => "no implicit H".to_string(),
            value if value > 0 => format!("implicit H >= {value}"),
            value => format!("implicit H: {}", query_count(value)),
        });
    }
    if let Some(value) = query.substitution_count {
        constraints.push(format!(
            "substitution count: {}",
            query_count_description(value)
        ));
    }
    if query.unsaturated == Some(true) {
        constraints.push("unsaturated".to_string());
    }
    if let Some(value) = query.ring_bond_count {
        constraints.push(format!(
            "ring-bond count: {}",
            query_count_description(value)
        ));
    }
    if let Some(value) = query.valence {
        constraints.push(format!("valence: {}", query_count(value)));
    }
    constraints
}

fn query_count_description(value: i8) -> String {
    match value {
        -2 => "as drawn".to_string(),
        -1 => "0".to_string(),
        value => value.to_string(),
    }
}

fn atom_query_symbol(query: &crate::ChemicalAtomQuery) -> Option<String> {
    query.elements.as_ref().map(|elements| {
        format!(
            "{}[{}]",
            if query.is_not_list { "NOT " } else { "" },
            elements.join(",")
        )
    })
}

fn atom_query_glyph_symbol(query: &crate::ChemicalAtomQuery) -> Option<String> {
    query
        .elements
        .as_ref()
        .map(|elements| format!("[{}]", elements.join(",")))
}

fn render_average_bond_length(mol: &RenderMol) -> f64 {
    let mut lengths = mol
        .bonds
        .iter()
        .filter(|bond| bond.kind.contributes_to_average_length())
        .filter_map(|bond| {
            let atom1 = mol.atoms.get(bond.atom1)?;
            let atom2 = mol.atoms.get(bond.atom2)?;
            let length = (atom2.x - atom1.x).hypot(atom2.y - atom1.y);
            (length.is_finite() && length > f64::EPSILON).then_some(length)
        })
        .collect::<Vec<_>>();
    median(&mut lengths).unwrap_or(1.0)
}

fn ctfile_point(
    x: f64,
    y: f64,
    atom_indexes: &[usize],
    mol: &RenderMol,
    average_length: f64,
) -> CtfilePoint {
    let atom_index = atom_indexes
        .iter()
        .copied()
        .min_by(|&left, &right| {
            let left_distance = (mol.atoms[left].x - x).hypot(mol.atoms[left].y - y);
            let right_distance = (mol.atoms[right].x - x).hypot(mol.atoms[right].y - y);
            left_distance.total_cmp(&right_distance)
        })
        .unwrap_or(0);
    CtfilePoint {
        atom_index,
        offset: [
            (x - mol.atoms[atom_index].x) / average_length,
            (y - mol.atoms[atom_index].y) / average_length,
        ],
    }
}

fn sgroup_depiction_label(group: &crate::ChemicalSGroup) -> Option<String> {
    if group.kind == "data" {
        let value = group
            .field_values
            .first()
            .or(group.field_value.as_ref())
            .cloned();
        return match (&group.field_name, value) {
            (Some(name), Some(value)) => Some(format!("{name}: {value}")),
            (None, Some(value)) => Some(value),
            (Some(name), None) => Some(name.clone()),
            (None, None) => group.label.clone(),
        };
    }
    group.label.clone().or_else(|| group.subscript.clone())
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
    let mut primary_lengths = Vec::new();
    let mut fallback_lengths = Vec::new();

    for (i, bond) in mol.bonds.iter().enumerate() {
        adj[bond.atom1].push((bond.atom2, bond.kind, i));
        adj[bond.atom2].push((bond.atom1, bond.kind, i));

        let u = &mol.atoms[bond.atom1];
        let v = &mol.atoms[bond.atom2];
        let dx = u.x - v.x;
        let dy = u.y - v.y;
        let length = dx.hypot(dy);
        if length.is_finite() && length > f64::EPSILON {
            fallback_lengths.push(length);
            if bond.kind.contributes_to_average_length() {
                primary_lengths.push(length);
            }
        }
    }

    let avg_length = median(&mut primary_lengths)
        .or_else(|| median(&mut fallback_lengths))
        .unwrap_or(1.0);

    let mut visited_nodes = vec![false; mol.atoms.len()];
    let mut handled_bonds = vec![false; mol.bonds.len()];
    let components = connected_components(&adj);
    let labels = build_labels(
        mol,
        &adj,
        mode_str,
        stereo_map,
        &mut visited_nodes,
        &mut handled_bonds,
    );
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
    let mut has_component = false;
    for component in components {
        let mut emitted_component = false;
        for start_node in component {
            if !visited_nodes[start_node] {
                if !emitted_component && has_component {
                    root_commands.push(Command::ComponentBreak);
                }
                dfs(
                    start_node,
                    &context,
                    &mut visited_nodes,
                    &mut handled_bonds,
                    &mut root_commands,
                );
                emitted_component = true;
            }
        }
        if emitted_component {
            has_component = true;
        }
    }

    root_commands
}

fn median(values: &mut [f64]) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    values.sort_by(f64::total_cmp);
    let middle = values.len() / 2;
    if values.len() % 2 == 0 {
        Some(values[middle - 1] / 2.0 + values[middle] / 2.0)
    } else {
        Some(values[middle])
    }
}

fn build_labels(
    mol: &RenderMol,
    adj: &[Vec<(usize, BondKind, usize)>],
    mode_str: &str,
    stereo_map: &HashMap<(usize, usize), (u8, bool)>,
    visited_nodes: &mut [bool],
    handled_bonds: &mut [bool],
) -> Vec<RenderLabel> {
    let mut labels = vec![RenderLabel::default(); mol.atoms.len()];

    for i in 0..mol.atoms.len() {
        let atom = &mol.atoms[i];

        if mode_str == "abbreviate" || mode_str == "skeletal" {
            if atom.element == "H" {
                if is_foldable_hydrogen(i, mol, adj, stereo_map) {
                    visited_nodes[i] = true;
                } else {
                    labels[i] = render_label(atom, atom.hydrogens);
                }
                continue;
            }

            let (explicit_h_count, explicit_h_bonds) =
                explicit_h_neighbors(i, mol, adj, stereo_map);
            for bond_idx in explicit_h_bonds {
                handled_bonds[bond_idx] = true;
            }

            let total_h = atom.hydrogens + explicit_h_count;
            if atom.element == "C" {
                let heavy_atoms = adj[i].len().saturating_sub(explicit_h_count as usize);
                if mode_str == "skeletal" {
                    if heavy_atoms == 0
                        || atom.stereo_annotation.is_some()
                        || atom.charge != 0
                        || has_visible_atom_metadata(atom)
                    {
                        labels[i] = render_label(atom, total_h);
                    }
                } else {
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

fn connected_components(adj: &[Vec<(usize, BondKind, usize)>]) -> Vec<Vec<usize>> {
    let mut components = Vec::new();
    let mut visited = vec![false; adj.len()];

    for start in 0..adj.len() {
        if visited[start] {
            continue;
        }

        let mut component = Vec::new();
        let mut queue = VecDeque::from([start]);
        visited[start] = true;
        while let Some(atom) = queue.pop_front() {
            component.push(atom);
            for &(neighbor, _, _) in &adj[atom] {
                if !visited[neighbor] {
                    visited[neighbor] = true;
                    queue.push_back(neighbor);
                }
            }
        }
        components.push(component);
    }

    components
}

fn has_visible_atom_metadata(atom: &RenderAtom) -> bool {
    atom.force_label
        || atom.element == "*"
        || atom.isotope.is_some()
        || atom.radical.is_some()
        || atom.atom_map.is_some()
        || atom.rgroup_label.is_some()
}

fn render_label(atom: &RenderAtom, hydrogen_count: u8) -> RenderLabel {
    let text = atom_label(&atom.element, hydrogen_count, atom.charge);
    let structured = (hydrogen_count > 1 || atom.charge != 0 || has_visible_atom_metadata(atom))
        .then(|| AtomLabel {
            symbol: atom.element.clone(),
            hydrogen_count,
            charge: atom.charge,
            isotope: atom.isotope,
            radical: atom.radical,
            atom_map: atom.atom_map,
            rgroup_label: atom.rgroup_label,
            query_symbol: atom.query_symbol,
            query_negated: atom.query_negated,
            hidden: atom.hide_label,
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
    stereo_map: &HashMap<(usize, usize), (u8, bool)>,
) -> (u8, Vec<usize>) {
    let mut count = 0u8;
    let mut bond_indices = Vec::new();
    for &(neighbor, _, bond_idx) in &adj[atom_idx] {
        if is_foldable_hydrogen(neighbor, mol, adj, stereo_map) {
            count = count.saturating_add(1);
            bond_indices.push(bond_idx);
        }
    }
    (count, bond_indices)
}

fn is_foldable_hydrogen(
    atom_idx: usize,
    mol: &RenderMol,
    adj: &[Vec<(usize, BondKind, usize)>],
    stereo_map: &HashMap<(usize, usize), (u8, bool)>,
) -> bool {
    mol.atoms[atom_idx].element == "H"
        && matches!(
            adj[atom_idx].as_slice(),
            [(neighbor, _, _)]
                if mol.atoms[*neighbor].element != "H"
                    && !stereo_map.contains_key(&(atom_idx, *neighbor))
        )
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

fn extract_sdf_stereo(mol: &sdfrust::Molecule, sdf: &str) -> HashMap<(usize, usize), (u8, bool)> {
    let mut map = HashMap::new();
    for bond in &mol.bonds {
        let stereo = match bond.stereo {
            BondStereo::Up => 1,
            BondStereo::Either => 4,
            BondStereo::Down => 6,
            BondStereo::None => continue,
        };
        if bond.atom1 < mol.atoms.len() && bond.atom2 < mol.atoms.len() {
            map.insert((bond.atom1, bond.atom2), (stereo, true));
            map.insert((bond.atom2, bond.atom1), (stereo, false));
        }
    }

    if mol.format_version == SdfFormat::V2000 {
        let lines = sdf.lines().collect::<Vec<_>>();
        let atom_count = lines
            .get(3)
            .and_then(|line| line.get(0..3))
            .and_then(|value| value.trim().parse::<usize>().ok())
            .unwrap_or(0);

        for (bond_idx, bond) in mol.bonds.iter().enumerate() {
            if bond.order != BondOrder::Double {
                continue;
            }
            let stereo = lines
                .get(4 + atom_count + bond_idx)
                .and_then(|line| line.get(9..12))
                .and_then(|value| value.trim().parse::<u8>().ok());
            if stereo == Some(3) && bond.atom1 < mol.atoms.len() && bond.atom2 < mol.atoms.len() {
                map.insert((bond.atom1, bond.atom2), (3, true));
                map.insert((bond.atom2, bond.atom1), (3, false));
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
    let angle = dy.atan2(dx).to_degrees();
    if angle.is_finite() {
        angle
    } else {
        0.0
    }
}

fn calc_length_scale(u: &RenderAtom, v: &RenderAtom, avg_len: f64) -> f64 {
    let dx = u.x - v.x;
    let dy = u.y - v.y;
    let len = dx.hypot(dy);
    let scale = len / avg_len;
    if scale.is_finite() && scale > f64::EPSILON {
        scale
    } else {
        1.0
    }
}

fn bond_func_name(
    kind: BondKind,
    u: usize,
    v: usize,
    coordination_forward: bool,
    stereo_map: &HashMap<(usize, usize), (u8, bool)>,
) -> &'static str {
    if kind == BondKind::Double && matches!(stereo_map.get(&(u, v)), Some(&(3 | 4, _))) {
        return "crossed-double";
    }

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
                4 => return "either",
                _ => {}
            }
        }
    }

    match kind {
        BondKind::Double => "double",
        BondKind::Triple => "triple",
        BondKind::Quadruple => "quadruple",
        BondKind::Single => "single",
        BondKind::Aromatic => "aromatic",
        BondKind::SingleOrDouble => "single-or-double",
        BondKind::SingleOrAromatic => "single-or-aromatic",
        BondKind::DoubleOrAromatic => "double-or-aromatic",
        BondKind::Any => "any",
        BondKind::Coordination if coordination_forward => "coordination-right",
        BondKind::Coordination => "coordination-left",
        BondKind::Hydrogen => "hydrogen",
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
        let offset = if kind.has_parallel_line() {
            calculate_double_bond_offset(u, v, context.mol, context.rings)
        } else {
            None
        };
        links.push(LinkData {
            target: format!("a{v}"),
            name: format!("b{bond_idx}"),
            bond_type: bond_func_name(
                kind,
                u,
                v,
                context.mol.bonds[bond_idx].atom1 == u,
                context.stereo_map,
            )
            .to_string(),
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

        let offset = if kind.has_parallel_line() {
            calculate_double_bond_offset(u, v, context.mol, context.rings)
        } else {
            None
        };

        let bond_cmd = Command::Bond {
            name: format!("b{bond_idx}"),
            bond_type: bond_func_name(
                kind,
                u,
                v,
                context.mol.bonds[bond_idx].atom1 == u,
                context.stereo_map,
            )
            .to_string(),
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

    fn atom_angle(coords: &[(f32, f32)], center: usize, atom: usize) -> f64 {
        let dx = f64::from(coords[atom].0 - coords[center].0);
        let dy = f64::from(coords[atom].1 - coords[center].1);
        dy.atan2(dx).to_degrees().rem_euclid(360.0)
    }

    fn assert_angle(actual: f64, expected: f64) {
        let difference = (actual - expected.rem_euclid(360.0)).abs();
        let circular_difference = difference.min(360.0 - difference);
        assert!(
            circular_difference < 1e-4,
            "expected {expected} degrees, got {actual}"
        );
    }

    fn first_fragment(commands: &[Command]) -> (&str, Option<&AtomLabel>) {
        match commands.first().unwrap() {
            Command::Fragment { element, atom, .. } => (element, atom.as_ref()),
            command => panic!("expected fragment, got {command:?}"),
        }
    }

    fn contains_fragment(commands: &[Command], expected: &str) -> bool {
        commands.iter().any(|command| match command {
            Command::Fragment { element, .. } => element == expected,
            Command::Branch { body } => contains_fragment(body, expected),
            Command::Bond { .. } | Command::ComponentBreak | Command::Ctfile { .. } => false,
        })
    }

    fn collect_bonds(commands: &[Command], output: &mut Vec<(String, f64, Option<String>)>) {
        for command in commands {
            match command {
                Command::Fragment { links, .. } => {
                    output.extend(links.iter().map(|link| {
                        (
                            link.bond_type.clone(),
                            link.length_scale,
                            link.offset.clone(),
                        )
                    }));
                }
                Command::Bond {
                    bond_type,
                    length_scale,
                    offset,
                    ..
                } => output.push((bond_type.clone(), *length_scale, offset.clone())),
                Command::Branch { body } => collect_bonds(body, output),
                Command::ComponentBreak | Command::Ctfile { .. } => {}
            }
        }
    }

    fn bond_data(commands: &[Command]) -> Vec<(String, f64, Option<String>)> {
        let mut output = Vec::new();
        collect_bonds(commands, &mut output);
        output
    }

    fn collect_fragments(commands: &[Command], output: &mut Vec<(String, String, Option<String>)>) {
        for command in commands {
            match command {
                Command::Fragment {
                    element,
                    name,
                    annotation,
                    ..
                } => output.push((element.clone(), name.clone(), annotation.clone())),
                Command::Branch { body } => collect_fragments(body, output),
                _ => {}
            }
        }
    }

    fn fragment_data(commands: &[Command]) -> Vec<(String, String, Option<String>)> {
        let mut output = Vec::new();
        collect_fragments(commands, &mut output);
        output
    }

    fn v2000_bond(order: u8, stereo: u8) -> String {
        format!(
            concat!(
                "bond semantics\n",
                "  molchemist\n",
                "\n",
                "  2  1  0  0  0  0  0  0  0  0999 V2000\n",
                "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
                "    1.5000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n",
                "  1  2{order:>3}{stereo:>3}  0  0  0\n",
                "M  END\n",
            ),
            order = order,
            stereo = stereo,
        )
    }

    fn degenerate_v2000(second_z: f64, stereo: u8) -> String {
        format!(
            concat!(
                "degenerate layout\n",
                "  molchemist\n",
                "\n",
                "  2  1  0  0  0  0  0  0  0  0999 V2000\n",
                "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
                "    0.0000    0.0000{second_z:>10.4} N   0  0  0  0  0  0  0  0  0  0  0  0\n",
                "  1  2  1{stereo:>3}  0  0  0\n",
                "M  END\n",
            ),
            second_z = second_z,
            stereo = stereo,
        )
    }

    const UNEVEN_V2000: &str = concat!(
        "uneven layout\n",
        "  molchemist\n",
        "\n",
        "  3  2  0  0  0  0  0  0  0  0999 V2000\n",
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
        "    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
        "  101.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n",
        "  1  2  1  0  0  0  0\n",
        "  2  3  1  0  0  0  0\n",
        "M  END\n",
    );

    fn v3000_bond(order: u8, atom1: usize, atom2: usize) -> String {
        format!(
            concat!(
                "bond semantics\n",
                "  molchemist\n",
                "\n",
                "  0  0  0     0  0            999 V3000\n",
                "M  V30 BEGIN CTAB\n",
                "M  V30 COUNTS 2 1 0 0 0\n",
                "M  V30 BEGIN ATOM\n",
                "M  V30 1 C 0.0000 0.0000 0.0000 0\n",
                "M  V30 2 N 1.5000 0.0000 0.0000 0\n",
                "M  V30 END ATOM\n",
                "M  V30 BEGIN BOND\n",
                "M  V30 1 {order} {atom1} {atom2}\n",
                "M  V30 END BOND\n",
                "M  V30 END CTAB\n",
                "M  END\n",
            ),
            order = order,
            atom1 = atom1,
            atom2 = atom2,
        )
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
    fn sdf_preserves_v2000_aromatic_and_query_bond_orders() {
        let cases = [
            (4, "aromatic"),
            (5, "single-or-double"),
            (6, "single-or-aromatic"),
            (7, "double-or-aromatic"),
            (8, "any"),
        ];

        for (order, expected) in cases {
            let commands = sdf_to_commands(&v2000_bond(order, 0), RenderMode::Full).unwrap();
            assert_eq!(bond_data(&commands)[0].0, expected, "bond order {order}");
        }
    }

    #[test]
    fn sdf_preserves_v3000_coordination_direction_and_hydrogen_bonds() {
        let forward = sdf_to_commands(&v3000_bond(9, 1, 2), RenderMode::Full).unwrap();
        let reverse = sdf_to_commands(&v3000_bond(9, 2, 1), RenderMode::Full).unwrap();
        let hydrogen = sdf_to_commands(&v3000_bond(10, 1, 2), RenderMode::Full).unwrap();

        assert_eq!(bond_data(&forward)[0].0, "coordination-right");
        assert_eq!(bond_data(&reverse)[0].0, "coordination-left");
        assert_eq!(bond_data(&hydrogen)[0].0, "hydrogen");
        assert!((bond_data(&hydrogen)[0].1 - 1.0).abs() < 1e-9);
    }

    #[test]
    fn sdf_preserves_either_stereochemistry_as_a_wavy_bond() {
        let commands = sdf_to_commands(&v2000_bond(1, 4), RenderMode::Full).unwrap();

        assert_eq!(bond_data(&commands)[0].0, "either");
    }

    #[test]
    fn sdf_preserves_undefined_double_bond_stereo_as_a_crossed_double() {
        let v2000 = sdf_to_commands(&v2000_bond(2, 3), RenderMode::Full).unwrap();
        let v3000 = concat!(
            "undefined double\n",
            "  molchemist\n",
            "\n",
            "  0  0  0     0  0            999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 2 1 0 0 0\n",
            "M  V30 BEGIN ATOM\n",
            "M  V30 1 C 0.0000 0.0000 0.0000 0\n",
            "M  V30 2 N 1.5000 0.0000 0.0000 0\n",
            "M  V30 END ATOM\n",
            "M  V30 BEGIN BOND\n",
            "M  V30 1 2 1 2 CFG=2\n",
            "M  V30 END BOND\n",
            "M  V30 END CTAB\n",
            "M  END\n",
        );
        let v3000 = sdf_to_commands(v3000, RenderMode::Full).unwrap();

        assert_eq!(bond_data(&v2000)[0].0, "crossed-double");
        assert_eq!(bond_data(&v3000)[0].0, "crossed-double");
    }

    #[test]
    fn abbreviated_modes_do_not_fold_a_stereochemical_hydrogen_bond() {
        let sdf = concat!(
            "stereochemical hydrogen\n",
            "  molchemist\n",
            "\n",
            "  5  4  0  0  0  0  0  0  0  0999 V2000\n",
            "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
            "    1.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n",
            "   -1.0000    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n",
            "    0.0000    1.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n",
            "    0.0000   -1.0000    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0\n",
            "  1  2  1  1  0  0  0\n",
            "  1  3  1  0  0  0  0\n",
            "  1  4  1  0  0  0  0\n",
            "  1  5  1  0  0  0  0\n",
            "M  END\n",
        );

        for mode in [RenderMode::Abbreviate, RenderMode::Skeletal] {
            let commands = sdf_to_commands(sdf, mode).unwrap();
            assert!(fragment_data(&commands)
                .iter()
                .any(|(element, name, _)| element == "H" && name == "a1"));
            assert!(bond_data(&commands)
                .iter()
                .any(|(bond_type, _, _)| bond_type.starts_with("cram-filled")));
        }
    }

    #[test]
    fn sdf_preserves_atom_parity_and_enhanced_stereo_groups_as_annotations() {
        let sdf = concat!(
            "enhanced stereo\n",
            "  molchemist\n",
            "\n",
            "  0  0  0     0  0            999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 2 1 0 0 0\n",
            "M  V30 BEGIN ATOM\n",
            "M  V30 10 C 0.0000 0.0000 0.0000 0 CFG=1\n",
            "M  V30 20 F 1.5000 0.0000 0.0000 0\n",
            "M  V30 END ATOM\n",
            "M  V30 BEGIN BOND\n",
            "M  V30 1 1 10 20\n",
            "M  V30 END BOND\n",
            "M  V30 BEGIN COLLECTION\n",
            "M  V30 MDLV30/STEREL7 ATOMS=(1 10)\n",
            "M  V30 END COLLECTION\n",
            "M  V30 END CTAB\n",
            "M  END\n",
        );

        let commands = sdf_to_commands(sdf, RenderMode::Skeletal).unwrap();
        let carbon = fragment_data(&commands)
            .into_iter()
            .find(|(_, name, _)| name == "a0")
            .unwrap();

        assert_eq!(carbon.0, "C");
        assert_eq!(carbon.2.as_deref(), Some("CFG=1; OR7"));
    }

    #[test]
    fn sdf_inspection_preserves_source_ids_properties_and_unsupported_features() {
        let sdf = concat!(
            "foundation\n",
            "  molchemist\n",
            "semantic record\n",
            "  0  0  0     0  0            999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 2 1 1 0 0\n",
            "M  V30 BEGIN ATOM\n",
            "M  V30 10 C 0.0000 0.0000 0.0000 7 MASS=13\n",
            "M  V30 20 O 1.5000 0.0000 0.0000 0\n",
            "M  V30 END ATOM\n",
            "M  V30 BEGIN BOND\n",
            "M  V30 7 1 10 20 TOPO=1 RXCTR=1\n",
            "M  V30 END BOND\n",
            "M  V30 BEGIN SGROUP\n",
            "M  V30 3 SUP ATOMS=(1 10) LABEL=\"Me\"\n",
            "M  V30 END SGROUP\n",
            "M  V30 BEGIN COLLECTION\n",
            "M  V30 MDLV30/HILITE ATOMS=(2 10 20) BONDS=(1 7)\n",
            "M  V30 END COLLECTION\n",
            "M  V30 END CTAB\n",
            "M  END\n",
            "> <NAME>\nfirst\nline two\n\n",
            "> <NAME>\nsecond\n\n",
            "$$$$\n",
        );

        let record = inspect_sdf_record(sdf, 1).unwrap();
        assert_eq!(record.schema_version, 1);
        assert_eq!(record.name, "foundation");
        assert_eq!(record.atoms[0].source_id, 10);
        assert_eq!(record.atoms[0].atom_map, Some(7));
        assert_eq!(record.atoms[0].isotope, Some(13));
        assert_eq!(record.bonds[0].source_id, 7);
        assert_eq!(record.bonds[0].atom1_source_id, 10);
        assert_eq!(record.bonds[0].atom2_source_id, 20);
        assert_eq!(record.sgroups[0].atom_source_ids, vec![10]);
        assert_eq!(record.collections[0].kind, "highlight");
        assert_eq!(record.collections[0].atom_source_ids, vec![10, 20]);
        assert_eq!(record.collections[0].bond_source_ids, vec![7]);
        assert_eq!(record.properties.len(), 2);
        assert_eq!(record.properties[0].value, "first\nline two");
        assert_eq!(record.properties[1].name, "NAME");
        assert_eq!(record.properties[1].value, "second");
        assert_eq!(
            record
                .diagnostics
                .iter()
                .map(|diagnostic| diagnostic.code.as_str())
                .collect::<Vec<_>>(),
            vec!["reaction-center-not-depicted",]
        );
        assert!(record.raw_record.starts_with("foundation\n"));
        assert!(!record.raw_record.contains("$$$$"));
    }

    #[test]
    fn ctfile_queries_sgroups_and_highlights_reach_the_depiction_ast() {
        let sdf = concat!(
            "ctfile depiction\n",
            "  molchemist\n",
            "\n",
            "  0  0  0     0  0            999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 4 2 1 0 0\n",
            "M  V30 BEGIN ATOM\n",
            "M  V30 10 \"NOT [C,N]\" 0.0000 0.0000 0.0000 0 HCOUNT=1 SUBST=3 UNSAT=1 RBCNT=2\n",
            "M  V30 20 R# 1.5000 0.0000 0.0000 0 RGROUPS=(1 7)\n",
            "M  V30 30 C 3.0000 0.0000 0.0000 0\n",
            "M  V30 40 C 4.5000 0.0000 0.0000 0\n",
            "M  V30 END ATOM\n",
            "M  V30 BEGIN BOND\n",
            "M  V30 5 1 10 20 TOPO=1\n",
            "M  V30 6 1 30 40\n",
            "M  V30 END BOND\n",
            "M  V30 BEGIN SGROUP\n",
            "M  V30 3 SRU 0 ATOMS=(2 30 40) LABEL=n BRKXYZ=(9 2.6 0.5 0 2.6 -0.5 0 0 0 0) BRKXYZ=(9 4.9 0.5 0 4.9 -0.5 0 0 0 0)\n",
            "M  V30 END SGROUP\n",
            "M  V30 BEGIN COLLECTION\n",
            "M  V30 MDLV30/HILITE ATOMS=(2 10 20) BONDS=(1 5)\n",
            "M  V30 END COLLECTION\n",
            "M  V30 END CTAB\n",
            "M  END\n",
        );

        let record = inspect_sdf_record(sdf, 1).unwrap();
        let query = record.atoms[0].query.as_ref().unwrap();
        assert_eq!(
            query.elements.as_deref(),
            Some(&["C".to_string(), "N".to_string()][..])
        );
        assert_eq!(query.hydrogen_count, Some(1));
        assert!(query.is_not_list);
        assert_eq!(query.substitution_count, Some(3));
        assert_eq!(query.unsaturated, Some(true));
        assert_eq!(query.ring_bond_count, Some(2));
        assert_eq!(record.atoms[1].rgroup_labels, vec![7]);
        assert_eq!(record.sgroups[0].brackets.len(), 2);
        assert!(record.diagnostics.is_empty());

        let commands = sdf_to_commands(sdf, RenderMode::Full).unwrap();
        assert!(commands.iter().any(|command| matches!(
            command,
            Command::Fragment {
                element,
                name,
                atom: Some(atom),
                ..
            } if element == "[C,N]"
                && name == "a0"
                && atom.query_symbol
                && atom.query_negated
        )));
        let rgroup_label = commands.iter().find_map(|command| match command {
            Command::Fragment {
                atom: Some(atom), ..
            } if atom.symbol == "R" && atom.rgroup_label == Some(7) => Some(atom.symbol.as_str()),
            _ => None,
        });
        assert_eq!(rgroup_label, Some("R"));
        match commands.last().unwrap() {
            Command::Ctfile {
                sgroups,
                highlights,
                atom_queries,
                bond_queries,
                ..
            } => {
                assert_eq!(sgroups.len(), 1);
                assert_eq!(sgroups[0].label.as_deref(), Some("n"));
                assert_eq!(highlights[0].atom_indexes, vec![0, 1]);
                assert_eq!(highlights[0].bonds[0].bond_index, 0);
                assert_eq!(atom_queries[0].symbol.as_deref(), Some("NOT [C,N]"));
                assert_eq!(
                    atom_queries[0].label,
                    "implicit H >= 1 · substitution count: 3 · unsaturated · ring-bond count: 2"
                );
                assert_eq!(atom_queries[0].compact_label, "(H1) (s3) (u) (r2)");
                assert_eq!(
                    atom_queries[0].constraint_lines,
                    [
                        "implicit H >= 1 · substitution count: 3",
                        "unsaturated · ring-bond count: 2",
                    ]
                );
                assert_eq!(bond_queries[0].label, "rn");
            }
            command => panic!("expected CTfile overlay, got {command:?}"),
        }
    }

    #[test]
    fn real_v2000_link_node_and_sap_metadata_are_normalized_and_depicted() {
        let link = include_str!("../../molchemist-cli/tests/fixtures/rdkit/Sgroups_Link_01.mol");
        let record = inspect_sdf_record(link, 1).unwrap();
        assert_eq!(record.link_nodes.len(), 1);
        assert_eq!(record.link_nodes[0].min_repeat, 1);
        assert_eq!(record.link_nodes[0].max_repeat, 3);
        assert_eq!(record.link_nodes[0].connections.len(), 2);
        assert!(record.diagnostics.is_empty());
        let commands = sdf_to_commands(link, RenderMode::Skeletal).unwrap();
        assert!(commands.iter().any(|command| matches!(
            command,
            Command::Ctfile { atom_annotations, .. }
                if atom_annotations.iter().any(|annotation|
                    annotation.kind == "link-node" && annotation.label == "1–3")
        )));

        let sap = include_str!("../../molchemist-cli/tests/fixtures/rdkit/sgroup_ap_bug.mol");
        let record = inspect_sdf_record(sap, 1).unwrap();
        assert_eq!(record.sgroups[0].attachment_points[0].atom_source_id, 1);
        assert_eq!(
            record.sgroups[0].attachment_points[0].leaving_atom_source_id,
            Some(5)
        );
        assert_eq!(record.sgroups[2].attachment_points.len(), 2);
        assert_eq!(record.sgroups[2].attachment_points[0].id, "Al");
        assert_eq!(record.sgroups[2].attachment_points[1].id, "Br");
        assert!(record.diagnostics.is_empty());
    }

    #[test]
    fn v3000_attachment_points_variable_bonds_and_link_nodes_reach_the_ast() {
        let sdf = concat!(
            "attachments\n  molchemist\n\n",
            "  0  0  0     0  0            999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 4 3 0 0 0\n",
            "M  V30 BEGIN ATOM\n",
            "M  V30 10 C 0 0 0 0 ATTCHPT=-1\n",
            "M  V30 20 C 1.5 0 0 0\n",
            "M  V30 30 C 2.25 1.25 0 0\n",
            "M  V30 40 C 2.25 -1.25 0 0\n",
            "M  V30 END ATOM\n",
            "M  V30 BEGIN BOND\n",
            "M  V30 7 1 10 20 ENDPTS=(3 20 30 40) ATTACH=ANY\n",
            "M  V30 8 1 20 30\n",
            "M  V30 9 1 20 40\n",
            "M  V30 END BOND\n",
            "M  V30 LINKNODE 1 4 2 20 30 20 40\n",
            "M  V30 END CTAB\n",
            "M  END\n",
        );
        let record = inspect_sdf_record(sdf, 1).unwrap();
        assert_eq!(record.atoms[0].attachment_points, vec![-1]);
        assert_eq!(record.bonds[0].endpoint_source_ids, vec![20, 30, 40]);
        assert_eq!(record.bonds[0].attachment_mode.as_deref(), Some("ANY"));
        assert_eq!(record.link_nodes[0].max_repeat, 4);
        assert_eq!(record.diagnostics.len(), 1);
        assert_eq!(
            record.diagnostics[0].code,
            "rgroup-attachment-context-not-supported"
        );

        let commands = sdf_to_commands(sdf, RenderMode::Skeletal).unwrap();
        assert!(commands.iter().any(|command| matches!(
            command,
            Command::Ctfile {
                atom_annotations,
                variable_attachments,
                ..
            } if atom_annotations.iter().all(|annotation| annotation.kind == "link-node")
                && atom_annotations.iter().any(|annotation| annotation.label == "1–4")
                && variable_attachments.iter().any(|attachment|
                    attachment.mode == "ANY" && attachment.endpoint_atom_indexes.len() == 3)
        )));
    }

    #[test]
    fn invalid_multicenter_bonds_and_link_nodes_are_diagnostic() {
        let sdf = concat!(
            "invalid attachment semantics\n  molchemist\n\n",
            "  0  0  0     0  0            999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 3 1 0 0 0\n",
            "M  V30 BEGIN ATOM\n",
            "M  V30 1 C 0 0 0 0\n",
            "M  V30 2 * 1.5 0 0 0\n",
            "M  V30 3 F 3 0 0 0\n",
            "M  V30 END ATOM\n",
            "M  V30 BEGIN BOND\n",
            "M  V30 1 1 2 3 ENDPTS=(1 1) ATTACH=MAYBE\n",
            "M  V30 END BOND\n",
            "M  V30 LINKNODE 2 4 1 2 3\n",
            "M  V30 END CTAB\n",
            "M  END\n",
        );

        let record = inspect_sdf_record(sdf, 1).unwrap();
        assert_eq!(record.bonds[0].endpoint_source_ids, vec![1]);
        assert_eq!(record.bonds[0].attachment_mode.as_deref(), Some("MAYBE"));
        assert_eq!(record.link_nodes[0].min_repeat, 2);
        assert_eq!(
            record
                .diagnostics
                .iter()
                .map(|diagnostic| diagnostic.code.as_str())
                .collect::<Vec<_>>(),
            vec![
                "invalid-multicenter-bond-semantics",
                "invalid-link-node-semantics"
            ]
        );

        let commands = sdf_to_commands(sdf, RenderMode::Skeletal).unwrap();
        assert!(commands.iter().all(|command| !matches!(
            command,
            Command::Ctfile {
                variable_attachments,
                ..
            } if !variable_attachments.is_empty()
        )));
    }

    #[test]
    fn standard_atom_attributes_drive_glyphs_and_user_collections_are_preserved() {
        let sdf = concat!(
            "collections\n  molchemist\n\n",
            "  0  0  0     0  0            999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 2 1 0 0 0\n",
            "M  V30 BEGIN ATOM\n",
            "M  V30 10 \"NOT [C,N]\" 0 0 0 0\n",
            "M  V30 20 R# 1.5 0 0 0 RGROUPS=(1 7)\n",
            "M  V30 END ATOM\n",
            "M  V30 BEGIN BOND\n",
            "M  V30 1 1 10 20\n",
            "M  V30 END BOND\n",
            "M  V30 BEGIN COLLECTION\n",
            "M  V30 MOLCHEMIST/SELECTION ATOMS=(2 10 20)\n",
            "M  V30 END COLLECTION\n",
            "M  V30 END CTAB\n",
            "M  END\n",
        );
        let record = inspect_sdf_record(sdf, 1).unwrap();
        assert_eq!(record.collections.len(), 1);
        assert_eq!(record.collections[0].name, "MOLCHEMIST/SELECTION");
        assert_eq!(record.collections[0].kind, "user-defined");
        assert_eq!(record.collections[0].atom_source_ids, vec![10, 20]);
        assert_eq!(record.diagnostics.len(), 1);
        assert_eq!(
            record.diagnostics[0].code,
            "user-defined-collection-not-depicted"
        );
        let commands = sdf_to_commands(sdf, RenderMode::Skeletal).unwrap();
        assert!(commands.iter().any(|command| matches!(
            command,
            Command::Fragment { atom: Some(atom), .. }
                if atom.symbol == "[C,N]" && atom.query_negated
        )));
        assert!(commands.iter().any(|command| matches!(
            command,
            Command::Fragment { atom: Some(atom), .. }
                if atom.symbol == "R" && atom.rgroup_label == Some(7)
        )));
    }

    #[test]
    fn uninterpreted_sgroups_and_non_graph_highlights_are_diagnostic() {
        let sdf = concat!(
            "future metadata\n  molchemist\n\n",
            "  0  0  0     0  0            999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 1 0 1 0 0\n",
            "M  V30 BEGIN ATOM\n",
            "M  V30 1 C 0 0 0 0\n",
            "M  V30 END ATOM\n",
            "M  V30 BEGIN SGROUP\n",
            "M  V30 3 ZZZ 0 ATOMS=(1 1)\n",
            "M  V30 END SGROUP\n",
            "M  V30 BEGIN COLLECTION\n",
            "M  V30 MDLV30/HILITE SGROUPS=(1 3)\n",
            "M  V30 END COLLECTION\n",
            "M  V30 END CTAB\n",
            "M  END\n",
        );

        let record = inspect_sdf_record(sdf, 1).unwrap();
        assert_eq!(record.sgroups[0].type_code, "ZZZ");
        assert_eq!(record.sgroups[0].kind, "unknown");
        assert_eq!(record.collections[0].sgroup_source_ids, vec![3]);
        assert_eq!(
            record
                .diagnostics
                .iter()
                .map(|diagnostic| diagnostic.code.as_str())
                .collect::<Vec<_>>(),
            vec!["unknown-sgroup-semantics", "highlight-members-not-depicted"]
        );
    }

    #[test]
    fn real_v2000_multi_atom_superatoms_contract_to_labelled_graph_nodes() {
        let sdf =
            include_str!("../../molchemist-cli/tests/fixtures/rdkit/Sgroups_Abbreviations.mol");
        let record = inspect_sdf_record(sdf, 1).unwrap();
        assert!(record.diagnostics.is_empty());

        let commands = sdf_to_commands(sdf, RenderMode::Full).unwrap();
        let fragments = fragment_data(&commands);
        assert_eq!(fragments.len(), 8);
        assert!(fragments.iter().any(|(element, _, _)| element == "NO₂"));
        assert!(fragments.iter().any(|(element, _, _)| element == "COOH"));
        assert_eq!(bond_data(&commands).len(), 8);
        assert!(!commands
            .iter()
            .any(|command| matches!(command, Command::Ctfile { .. })));
    }

    #[test]
    fn contracted_superatom_remaps_crossing_and_internal_bond_highlights() {
        let sdf = concat!(
            "contracted highlight\n",
            "  molchemist\n",
            "\n",
            "  0  0  0     0  0            999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 3 2 1 0 0\n",
            "M  V30 BEGIN ATOM\n",
            "M  V30 10 C 0.0000 0.0000 0.0000 0\n",
            "M  V30 20 C 1.5000 0.0000 0.0000 0\n",
            "M  V30 30 O 3.0000 0.0000 0.0000 0\n",
            "M  V30 END ATOM\n",
            "M  V30 BEGIN BOND\n",
            "M  V30 1 1 10 20\n",
            "M  V30 2 1 20 30\n",
            "M  V30 END BOND\n",
            "M  V30 BEGIN SGROUP\n",
            "M  V30 4 SUP 0 ATOMS=(2 20 30) XBONDS=(1 1) LABEL=Et\n",
            "M  V30 END SGROUP\n",
            "M  V30 BEGIN COLLECTION\n",
            "M  V30 MDLV30/HILITE ATOMS=(1 30) BONDS=(2 1 2)\n",
            "M  V30 END COLLECTION\n",
            "M  V30 END CTAB\n",
            "M  END\n",
        );

        let commands = sdf_to_commands(sdf, RenderMode::Full).unwrap();
        let fragments = fragment_data(&commands);
        assert_eq!(fragments.len(), 2);
        assert!(fragments
            .iter()
            .any(|(element, name, _)| { element == "Et" && name == "a1" }));
        assert_eq!(bond_data(&commands).len(), 1);

        match commands.last().unwrap() {
            Command::Ctfile {
                sgroups,
                highlights,
                ..
            } => {
                assert!(sgroups.is_empty());
                assert_eq!(highlights.len(), 1);
                assert_eq!(highlights[0].atom_indexes, vec![1]);
                assert_eq!(highlights[0].bonds.len(), 1);
                assert_eq!(highlights[0].bonds[0].bond_index, 0);
                assert_eq!(highlights[0].bonds[0].atom1_index, 0);
                assert_eq!(highlights[0].bonds[0].atom2_index, 1);
            }
            command => panic!("expected CTfile overlay, got {command:?}"),
        }
    }

    #[test]
    fn multi_attachment_superatom_uses_boundary_atom_centroid() {
        let sdf = concat!(
            "two attachments\n",
            "  molchemist\n",
            "\n",
            "  0  0  0     0  0            999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 4 3 1 0 0\n",
            "M  V30 BEGIN ATOM\n",
            "M  V30 10 C 0.0000 0.0000 0.0000 0\n",
            "M  V30 20 C 1.5000 0.0000 0.0000 0\n",
            "M  V30 30 C 3.0000 0.0000 0.0000 0\n",
            "M  V30 40 C 4.5000 0.0000 0.0000 0\n",
            "M  V30 END ATOM\n",
            "M  V30 BEGIN BOND\n",
            "M  V30 1 1 10 20\n",
            "M  V30 2 1 20 30\n",
            "M  V30 3 1 30 40\n",
            "M  V30 END BOND\n",
            "M  V30 BEGIN SGROUP\n",
            "M  V30 4 SUP 0 ATOMS=(2 20 30) XBONDS=(2 1 3) LABEL=CH2\n",
            "M  V30 END SGROUP\n",
            "M  V30 END CTAB\n",
            "M  END\n",
        );
        let parsed = parse_sdf_record(sdf, 1).unwrap();
        let mol = render_mol_from_record(&parsed.chemical, vec![None; 4]);
        let projection = project_record_render(&parsed.chemical, mol);

        assert_eq!(projection.mol.atoms.len(), 3);
        assert_eq!(projection.mol.atoms[1].element, "CH₂");
        assert!((projection.mol.atoms[1].x - 2.25).abs() <= 1e-8);
        assert!(projection.mol.atoms[1].y.abs() <= 1e-8);
        assert_eq!(projection.mol.bonds.len(), 2);
        assert_eq!(projection.bond_indexes, vec![Some(0), None, Some(1)]);
    }

    #[test]
    fn explicitly_expanded_multi_atom_superatom_keeps_its_source_graph() {
        let sdf = concat!(
            "expanded superatom\n",
            "  molchemist\n",
            "\n",
            "  0  0  0     0  0            999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 3 2 1 0 0\n",
            "M  V30 BEGIN ATOM\n",
            "M  V30 10 C 0.0000 0.0000 0.0000 0\n",
            "M  V30 20 C 1.5000 0.0000 0.0000 0\n",
            "M  V30 30 O 3.0000 0.0000 0.0000 0\n",
            "M  V30 END ATOM\n",
            "M  V30 BEGIN BOND\n",
            "M  V30 1 1 10 20\n",
            "M  V30 2 1 20 30\n",
            "M  V30 END BOND\n",
            "M  V30 BEGIN SGROUP\n",
            "M  V30 4 SUP 0 ATOMS=(2 20 30) XBONDS=(1 1) LABEL=Et ESTATE=E\n",
            "M  V30 END SGROUP\n",
            "M  V30 END CTAB\n",
            "M  END\n",
        );

        let record = inspect_sdf_record(sdf, 1).unwrap();
        assert_eq!(record.sgroups[0].expanded, Some(true));
        assert!(record.diagnostics.is_empty());
        let commands = sdf_to_commands(sdf, RenderMode::Full).unwrap();
        assert_eq!(fragment_data(&commands).len(), 3);
        assert_eq!(bond_data(&commands).len(), 2);
        match commands.last().unwrap() {
            Command::Ctfile { sgroups, .. } => assert_eq!(sgroups.len(), 1),
            command => panic!("expected expanded SGroup overlay, got {command:?}"),
        }
    }

    #[test]
    fn superatom_labels_use_subscript_digits_without_changing_leading_numbers() {
        assert_eq!(typographic_superatom_label("NO2"), "NO₂");
        assert_eq!(typographic_superatom_label("C12H25"), "C₁₂H₂₅");
        assert_eq!(typographic_superatom_label("1,2-Et2"), "1,2-Et₂");
    }

    #[test]
    fn v2000_atom_lists_rgroups_and_sgroups_are_normalized() {
        let sdf = concat!(
            "v2000 fidelity\n",
            "  molchemist\n",
            "\n",
            "  2  1  0  0  0  0  0  0  0  0  1 V2000\n",
            "    0.0000    0.0000    0.0000 L   0  0  0  1  0  0  0  0  0  0  0  0\n",
            "    1.5000    0.0000    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n",
            "  1  2  1  0  0  0  0\n",
            "M  ALS   1  2 F C   N   \n",
            "M  RGP  1   2   5\n",
            "M  STY  1   3 SRU\n",
            "M  SST  1   3 ALT\n",
            "M  SCN  1   3 HT\n",
            "M  SBT  1   3   1\n",
            "M  SAL   3  2   1   2\n",
            "M  SBL   3  1   1\n",
            "M  SMT   3 n\n",
            "M  SDI   3  4   -0.4000    0.5000   -0.4000   -0.5000\n",
            "M  SDI   3  4    1.9000    0.5000    1.9000   -0.5000\n",
            "M  END\n",
        );

        let record = inspect_sdf_record(sdf, 1).unwrap();
        let query = record.atoms[0].query.as_ref().unwrap();
        assert_eq!(query.elements, Some(vec!["C".to_string(), "N".to_string()]));
        assert_eq!(query.hydrogen_count, Some(0));
        assert_eq!(record.atoms[1].rgroup_labels, vec![5]);
        assert_eq!(record.sgroups[0].kind, "structure-repeat-unit");
        assert_eq!(record.sgroups[0].atom_source_ids, vec![1, 2]);
        assert_eq!(record.sgroups[0].crossing_bond_source_ids, vec![1]);
        assert_eq!(record.sgroups[0].label.as_deref(), Some("n"));
        assert_eq!(record.sgroups[0].subtype.as_deref(), Some("ALT"));
        assert_eq!(record.sgroups[0].connectivity.as_deref(), Some("HT"));
        assert_eq!(record.sgroups[0].bracket_type.as_deref(), Some("1"));
        assert_eq!(record.sgroups[0].brackets.len(), 2);
        assert!(record.diagnostics.is_empty());

        let commands = sdf_to_commands(sdf, RenderMode::Full).unwrap();
        match commands.last().unwrap() {
            Command::Ctfile { atom_queries, .. } => {
                assert_eq!(atom_queries[0].label, "no implicit H");
                assert_eq!(atom_queries[0].compact_label, "(H0)");
            }
            command => panic!("expected CTfile overlay, got {command:?}"),
        }
    }

    #[test]
    fn invalid_v3000_query_count_is_reported() {
        let sdf = concat!(
            "invalid v3000 query\n",
            "  molchemist\n",
            "\n",
            "  0  0  0     0  0            999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 1 0 0 0 0\n",
            "M  V30 BEGIN ATOM\n",
            "M  V30 1 C 0.0000 0.0000 0.0000 0 HCOUNT=6 SUBST=-3 UNSAT=2 RBCNT=-3\n",
            "M  V30 END ATOM\n",
            "M  V30 END CTAB\n",
            "M  END\n",
        );

        let record = inspect_sdf_record(sdf, 1).unwrap();
        assert_eq!(record.diagnostics.len(), 1);
        assert_eq!(record.diagnostics[0].code, "invalid-v3000-query-count");
        assert_eq!(record.diagnostics[0].atom_source_ids, vec![1]);
    }

    #[test]
    fn hydrogen_bonds_do_not_shrink_covalent_bond_lengths() {
        let sdf = concat!(
            "hydrogen bond scale\n",
            "  molchemist\n",
            "\n",
            "  0  0  0     0  0            999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 3 2 0 0 0\n",
            "M  V30 BEGIN ATOM\n",
            "M  V30 1 C 0.0000 0.0000 0.0000 0\n",
            "M  V30 2 N 1.5000 0.0000 0.0000 0\n",
            "M  V30 3 O 46.5000 0.0000 0.0000 0\n",
            "M  V30 END ATOM\n",
            "M  V30 BEGIN BOND\n",
            "M  V30 1 1 1 2\n",
            "M  V30 2 10 2 3\n",
            "M  V30 END BOND\n",
            "M  V30 END CTAB\n",
            "M  END\n",
        );

        let commands = sdf_to_commands(sdf, RenderMode::Full).unwrap();
        let bonds = bond_data(&commands);

        assert_eq!(bonds[0].0, "single");
        assert!((bonds[0].1 - 1.0).abs() < 1e-9);
        assert_eq!(bonds[1].0, "hydrogen");
        assert!((bonds[1].1 - 30.0).abs() < 1e-9);
        assert!(sdf_record_to_layout_input(sdf.as_bytes(), 1)
            .unwrap()
            .is_empty());
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
    fn sdf_requests_layout_only_for_unusable_source_coordinates() {
        let valid = v2000_bond(1, 0);
        let collapsed = degenerate_v2000(0.0, 0);
        let z_only = degenerate_v2000(1.5, 0);

        assert!(sdf_record_to_layout_input(valid.as_bytes(), 1)
            .unwrap()
            .is_empty());
        assert!(!sdf_record_to_layout_input(UNEVEN_V2000.as_bytes(), 1)
            .unwrap()
            .is_empty());
        for sdf in [&collapsed, &z_only] {
            let payload = sdf_record_to_layout_input(sdf.as_bytes(), 1).unwrap();
            let (atoms, bonds, atom_stereo, double_bond_stereo) = decode_layout_input(&payload);

            assert_eq!(atoms, vec![6, 7]);
            assert_eq!(bonds, vec![(0, 1, 1)]);
            assert!(atom_stereo.is_empty());
            assert!(double_bond_stereo.is_empty());
        }
    }

    #[test]
    fn sdf_fallback_coordinates_keep_source_stereo_and_nonzero_scale() {
        let sdf = degenerate_v2000(0.0, 1);
        let coords = coordinate_payload(&[(0.0, 0.0), (1.5, 0.0)], 1);
        let commands =
            sdf_record_to_commands_with_coords(&sdf, &coords, RenderMode::Full, 1).unwrap();
        let bonds = bond_data(&commands);

        assert_eq!(bonds.len(), 1);
        assert_eq!(bonds[0].0, "cram-filled-left");
        assert!((bonds[0].1 - 1.0).abs() < 1e-9);
    }

    #[test]
    fn coordinate_payload_rejects_non_finite_layout_results() {
        let sdf = degenerate_v2000(0.0, 0);
        let coords = coordinate_payload(&[(f32::NAN, 0.0), (1.5, 0.0)], 1);

        assert_eq!(
            sdf_record_to_commands_with_coords(&sdf, &coords, RenderMode::Full, 1).unwrap_err(),
            "Coordinate payload atom 1 has non-finite coordinates"
        );
    }

    #[test]
    fn sdf_preserves_disconnected_components_with_a_boundary() {
        let sdf = concat!(
            "sodium chloride\n",
            "  molchemist\n",
            "\n",
            "  2  0  0  0  0  0  0  0  0  0999 V2000\n",
            "    0.0000    0.0000    0.0000 Na  0  0  0  0  0  0  0  0  0  0  0  0\n",
            "    3.0000    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n",
            "M  CHG  2   1   1   2  -1\n",
            "M  END\n",
        );

        let commands = sdf_to_commands(sdf, RenderMode::Abbreviate).unwrap();
        assert!(matches!(
            commands.as_slice(),
            [
                Command::Fragment { element: sodium, name: sodium_name, .. },
                Command::ComponentBreak,
                Command::Fragment { element: chloride, name: chloride_name, .. },
            ] if sodium == "Na^+"
                && sodium_name == "a0"
                && chloride == "Cl^-"
                && chloride_name == "a1"
        ));
    }

    #[test]
    fn sdf_component_order_uses_the_first_source_atom_even_when_it_is_folded() {
        let sdf = concat!(
            "interleaved components\n",
            "  molchemist\n",
            "\n",
            "  3  1  0  0  0  0  0  0  0  0999 V2000\n",
            "    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n",
            "    4.0000    0.0000    0.0000 Na  0  0  0  0  0  0  0  0  0  0  0  0\n",
            "    1.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n",
            "  1  3  1  0  0  0  0\n",
            "M  END\n",
        );

        let commands = sdf_to_commands(sdf, RenderMode::Abbreviate).unwrap();
        assert!(matches!(
            commands.as_slice(),
            [
                Command::Fragment { element: water, name: oxygen_name, .. },
                Command::ComponentBreak,
                Command::Fragment { element: sodium, name: sodium_name, .. },
            ] if water == "OH"
                && oxygen_name == "a2"
                && sodium == "Na"
                && sodium_name == "a1"
        ));
    }

    #[test]
    fn hydrogen_bridging_heavy_atoms_is_not_folded_or_split() {
        let sdf = concat!(
            "bridging hydrogen\n",
            "  molchemist\n",
            "\n",
            "  3  2  0  0  0  0  0  0  0  0999 V2000\n",
            "   -1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
            "    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n",
            "    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
            "  1  2  1  0  0  0  0\n",
            "  2  3  1  0  0  0  0\n",
            "M  END\n",
        );

        let commands = sdf_to_commands(sdf, RenderMode::Abbreviate).unwrap();
        assert!(contains_fragment(&commands, "H"));
        assert!(!commands
            .iter()
            .any(|command| matches!(command, Command::ComponentBreak)));
    }

    #[test]
    fn disconnected_smiles_keeps_isolated_hydrogen_and_carbon_visible() {
        let coords = coordinate_payload(&[(0.0, 0.0), (2.0, 0.0), (4.0, 0.0)], 0);
        let commands =
            smiles_to_commands_with_coords("[H+].C.[Cl-]", &coords, RenderMode::Skeletal).unwrap();

        assert!(matches!(
            commands.as_slice(),
            [
                Command::Fragment {
                    element: hydrogen,
                    name: hydrogen_name,
                    atom: Some(hydrogen_label),
                    ..
                },
                Command::ComponentBreak,
                Command::Fragment {
                    element: carbon,
                    name: carbon_name,
                    atom: Some(carbon_label),
                    ..
                },
                Command::ComponentBreak,
                Command::Fragment {
                    element: chlorine,
                    name: chlorine_name,
                    atom: Some(chlorine_label),
                    ..
                },
            ] if hydrogen == "H^+"
                && hydrogen_name == "a0"
                && hydrogen_label.symbol == "H"
                && hydrogen_label.charge == 1
                && carbon == "CH_4"
                && carbon_name == "a1"
                && carbon_label.symbol == "C"
                && carbon_label.hydrogen_count == 4
                && chlorine == "Cl^-"
                && chlorine_name == "a2"
                && chlorine_label.symbol == "Cl"
                && chlorine_label.charge == -1
        ));
    }

    #[test]
    fn full_smiles_mode_keeps_component_boundaries_after_expanding_hydrogens() {
        let coords = coordinate_payload(
            &[
                (0.0, 0.0),
                (4.0, 0.0),
                (-1.0, 0.0),
                (0.0, 1.0),
                (1.0, 0.0),
                (0.0, -1.0),
            ],
            4,
        );
        let commands =
            smiles_to_commands_with_coords("C.[Cl-]", &coords, RenderMode::Full).unwrap();

        let breaks = commands
            .iter()
            .filter(|command| matches!(command, Command::ComponentBreak))
            .count();
        let break_index = commands
            .iter()
            .position(|command| matches!(command, Command::ComponentBreak))
            .unwrap();
        assert_eq!(breaks, 1);
        assert!(matches!(
            commands.get(break_index + 1),
            Some(Command::Fragment { element, name, .. })
                if element == "Cl^-" && name == "a1"
        ));
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
        let atom = atom.unwrap();
        assert_eq!(atom.symbol, "C");
        assert_eq!(atom.hydrogen_count, 2);
        assert_eq!(atom.charge, -1);
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
    fn smiles_quadruple_bond_is_preserved_for_layout_and_rendering() {
        for mode in [
            RenderMode::Full,
            RenderMode::Abbreviate,
            RenderMode::Skeletal,
        ] {
            let payload = smiles_layout_input("[Cr]$[Cr]", mode).unwrap();
            let (atoms, bonds, atom_stereo, double_bond_stereo) = decode_layout_input(&payload);
            assert_eq!(atoms, vec![24, 24]);
            assert_eq!(bonds, vec![(0, 1, 4)]);
            assert!(atom_stereo.is_empty());
            assert!(double_bond_stereo.is_empty());

            let coords = coordinate_payload(&[(0.0, 0.0), (1.0, 0.0)], 1);
            let commands = smiles_to_commands_with_coords("[Cr]$[Cr]", &coords, mode).unwrap();
            let mut rendered_bonds = Vec::new();
            collect_bonds(&commands, &mut rendered_bonds);
            assert_eq!(rendered_bonds, vec![("quadruple".to_string(), 1.0, None)]);
        }
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
    fn aromatic_heteroatoms_do_not_receive_kekule_double_bonds() {
        for smiles in ["o1cccc1", "c1occc1", "[nH]1cccc1", "c1cc[nH]c1"] {
            let graph = parse_smiles_graph(smiles).unwrap();
            let hetero = graph
                .atoms
                .iter()
                .position(|atom| atom.element == "O" || atom.element == "N")
                .unwrap();
            assert!(
                graph.bonds.iter().all(|bond| {
                    (bond.atom1 != hetero && bond.atom2 != hetero)
                        || bond.kind == SmilesBondKind::Single
                }),
                "{smiles}"
            );
            assert_eq!(
                graph
                    .bonds
                    .iter()
                    .filter(|bond| bond.kind == SmilesBondKind::Double)
                    .count(),
                2,
                "{smiles}"
            );
        }
    }

    #[test]
    fn mixed_aromatic_and_explicit_kekule_bonds_are_supported() {
        let graph = parse_smiles_graph("OCCc1c(C)[n+](=cs1)Cc2cnc(C)nc(N)2").unwrap();
        assert_eq!(graph.atoms.len(), 18);
        assert!(graph
            .bonds
            .iter()
            .any(|bond| bond.atom1 == 6 && bond.atom2 == 7 && bond.kind == SmilesBondKind::Double));
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

        let ast_bytes = smiles_to_ast(b"[N+]C", &coords, b"abbreviate").unwrap();
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
            rgroup_label: None,
            query_symbol: false,
            query_negated: false,
            stereo_annotation: None,
            force_label: false,
            hide_label: false,
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
    fn tetrahedral_chirality_keeps_incoming_atom_before_bracket_hydrogen() {
        let payload = smiles_to_layout_input(b"N[C@@H](C)C(=O)O").unwrap();
        let (_, _, atom_stereo, _) = decode_layout_input(&payload);

        assert_eq!(atom_stereo.len(), 1);
        let (atom, looking_from, atom1, atom2, direction) = atom_stereo[0];
        assert_eq!(atom, 1);
        assert_eq!(looking_from, 0);
        assert_eq!(atom1, u32::MAX);
        assert_eq!(atom2, 2);
        assert_eq!(direction, 1);
    }

    #[test]
    fn implicit_and_explicit_bracket_hydrogen_use_equivalent_local_order() {
        let implicit = smiles_to_layout_input(b"N[C@@H](C)C(=O)O").unwrap();
        let explicit = smiles_to_layout_input(b"N[C@@]([H])(C)C(=O)O").unwrap();
        let (_, _, implicit_stereo, _) = decode_layout_input(&implicit);
        let (_, _, explicit_stereo, _) = decode_layout_input(&explicit);

        assert_eq!(implicit_stereo, vec![(1, 0, u32::MAX, 2, 1)]);
        assert_eq!(explicit_stereo, vec![(1, 0, 2, 3, 1)]);
    }

    #[test]
    fn bracket_hydrogen_position_follows_component_and_branch_incoming_bonds() {
        let cases = [
            ("[C@H](F)(Cl)Br", 0, 0),
            ("N[C@@H](C)C(=O)O", 1, 1),
            ("C([C@H](F)Cl)O", 1, 1),
            ("[Na+].[C@H](F)(Cl)Br", 1, 0),
            ("N[C@H]1CCCCC1", 1, 1),
        ];

        for (smiles, atom, expected_position) in cases {
            let (_, positions) =
                smiles_neighbor_order(smiles, parse_smiles(smiles).unwrap().nodes().len()).unwrap();
            assert_eq!(positions[atom], expected_position, "{smiles}");
        }
    }

    #[test]
    fn leading_tetrahedral_center_places_implicit_hydrogen_first() {
        let payload = smiles_to_layout_input(b"[C@H](F)(Cl)Br").unwrap();
        let (_, _, atom_stereo, _) = decode_layout_input(&payload);

        assert_eq!(atom_stereo, vec![(0, u32::MAX, 1, 2, 2)]);
    }

    #[test]
    fn tetrahedral_ring_closure_uses_the_ring_tokens_lexical_position() {
        let graph = parse_smiles_graph("[C@]1(Br)(Cl)CCCC(F)C1").unwrap();
        let stereo = graph.atoms[0].stereo.unwrap();

        assert_eq!(stereo.looking_from, Some(8));
        assert_eq!(stereo.atom1, Some(1));
        assert_eq!(stereo.atom2, Some(2));
    }

    #[test]
    fn full_layout_replaces_chiral_implicit_hydrogen_with_its_explicit_atom() {
        let payload = smiles_to_full_layout_input(b"[C@H](F)(Cl)Br").unwrap();
        let (atoms, bonds, atom_stereo, _) = decode_layout_input(&payload);

        assert_eq!(atoms.len(), 5);
        assert_eq!(bonds.last(), Some(&(0, 4, 1)));
        assert_eq!(atom_stereo, vec![(0, 4, 1, 2, 2)]);
    }

    #[test]
    fn full_layout_keeps_nonleading_bracket_hydrogen_after_incoming_atom() {
        let payload = smiles_to_full_layout_input(b"N[C@@H](C)C(=O)O").unwrap();
        let (atoms, _, atom_stereo, _) = decode_layout_input(&payload);
        let [(center, looking_from, hydrogen, atom2, direction)] = atom_stereo.as_slice() else {
            panic!("expected one tetrahedral center");
        };

        assert_eq!(*center, 1);
        assert_eq!(*looking_from, 0);
        assert_eq!(atoms[*hydrogen as usize], 1);
        assert_eq!(*atom2, 2);
        assert_eq!(*direction, 1);
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
    fn square_planar_chirality_is_rendered_without_a_text_annotation() {
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
            .and_then(|(_, value)| value.as_text());

        assert_eq!(first_label, "Pt");
        assert_eq!(first_annotation, None);
    }

    #[test]
    fn square_planar_permutations_produce_u_four_and_z_paths() {
        let cases = [
            ("[Pt@SP1](F)(Cl)(Br)I", [90.0, 180.0, 270.0, 0.0]),
            ("[Pt@SP2](F)(Cl)(Br)I", [90.0, 270.0, 0.0, 180.0]),
            ("[Pt@SP3](F)(Cl)(Br)I", [90.0, 180.0, 0.0, 270.0]),
        ];

        for (smiles, expected) in cases {
            let graph = parse_smiles_graph(smiles).unwrap();
            let mut coords = vec![(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (-1.0, 0.0), (0.0, -1.0)];
            let mut stereo_map = HashMap::new();
            let depicted = apply_extended_stereo_depictions(&graph, &mut coords, &mut stereo_map);

            assert_eq!(depicted, HashSet::from([0]));
            assert!(stereo_map.is_empty());
            for (ligand, expected_angle) in (1..=4).zip(expected) {
                assert_angle(atom_angle(&coords, 0, ligand), expected_angle);
            }
        }
    }

    #[test]
    fn trigonal_bipyramidal_permutation_sets_axis_and_winding() {
        let cases = [
            ("[As@TB5](F)(Cl)(Br)(N)S", [90.0, 210.0, 330.0]),
            ("[As@TB6](F)(Cl)(Br)(N)S", [90.0, 330.0, 210.0]),
        ];

        for (smiles, equatorial_angles) in cases {
            let graph = parse_smiles_graph(smiles).unwrap();
            let mut coords = vec![
                (0.0, 0.0),
                (1.0, 0.0),
                (0.0, 1.0),
                (-1.0, 0.0),
                (0.0, -1.0),
                (1.0, 1.0),
            ];
            let mut stereo_map = HashMap::new();
            let depicted = apply_extended_stereo_depictions(&graph, &mut coords, &mut stereo_map);

            assert_eq!(depicted, HashSet::from([0]));
            assert_eq!(stereo_map.get(&(0, 1)), Some(&(1, true)));
            assert_eq!(stereo_map.get(&(0, 3)), Some(&(6, true)));
            for (ligand, expected_angle) in [2, 4, 5].into_iter().zip(equatorial_angles) {
                assert_angle(atom_angle(&coords, 0, ligand), expected_angle);
            }
        }
    }

    #[test]
    fn octahedral_permutation_sets_axis_shape_and_winding() {
        let graph = parse_smiles_graph("[Co@OH5](F)(Cl)(Br)(I)(N)S").unwrap();
        let mut coords = vec![
            (0.0, 0.0),
            (1.0, 0.0),
            (0.0, 1.0),
            (-1.0, 0.0),
            (0.0, -1.0),
            (1.0, 1.0),
            (-1.0, -1.0),
        ];
        let mut stereo_map = HashMap::new();
        let depicted = apply_extended_stereo_depictions(&graph, &mut coords, &mut stereo_map);

        assert_eq!(depicted, HashSet::from([0]));
        assert_eq!(stereo_map.get(&(0, 1)), Some(&(1, true)));
        assert_eq!(stereo_map.get(&(0, 5)), Some(&(6, true)));
        for (ligand, expected_angle) in [2, 3, 4, 6].into_iter().zip([90.0, 180.0, 0.0, 270.0]) {
            assert_angle(atom_angle(&coords, 0, ligand), expected_angle);
        }
    }

    #[test]
    fn every_tb_and_oh_permutation_has_a_native_depiction() {
        for number in 1..=20 {
            let graph = parse_smiles_graph(&format!("[As@TB{number}](F)(Cl)(Br)(N)S")).unwrap();
            let mut coords = vec![
                (0.0, 0.0),
                (1.0, 0.0),
                (0.0, 1.0),
                (-1.0, 0.0),
                (0.0, -1.0),
                (1.0, 1.0),
            ];
            let mut stereo_map = HashMap::new();

            assert_eq!(
                apply_extended_stereo_depictions(&graph, &mut coords, &mut stereo_map),
                HashSet::from([0]),
                "@TB{number}"
            );
            assert_eq!(stereo_map.len(), 4, "@TB{number}");
        }

        for number in 1..=30 {
            let graph = parse_smiles_graph(&format!("[Co@OH{number}](F)(Cl)(Br)(I)(N)S")).unwrap();
            let mut coords = vec![
                (0.0, 0.0),
                (1.0, 0.0),
                (0.0, 1.0),
                (-1.0, 0.0),
                (0.0, -1.0),
                (1.0, 1.0),
                (-1.0, -1.0),
            ];
            let mut stereo_map = HashMap::new();

            assert_eq!(
                apply_extended_stereo_depictions(&graph, &mut coords, &mut stereo_map),
                HashSet::from([0]),
                "@OH{number}"
            );
            assert_eq!(stereo_map.len(), 4, "@OH{number}");
        }
    }

    #[test]
    fn allene_variants_reverse_terminal_wedges() {
        for (chirality, expected) in [("AL1", (6, 1)), ("AL2", (1, 6))] {
            let graph = parse_smiles_graph(&format!("NC(Br)=[C@{chirality}]=C(O)C")).unwrap();
            let mut coords = vec![(0.0, 0.0); graph.atoms.len()];
            let mut stereo_map = HashMap::new();
            let depicted = apply_extended_stereo_depictions(&graph, &mut coords, &mut stereo_map);

            assert_eq!(depicted, HashSet::from([3]));
            assert_eq!(stereo_map.get(&(4, 5)), Some(&(expected.0, true)));
            assert_eq!(stereo_map.get(&(4, 6)), Some(&(expected.1, true)));
        }
    }

    #[test]
    fn allene_implicit_hydrogen_keeps_the_same_parity_in_full_and_skeletal_modes() {
        let graph = parse_smiles_graph("NC(Br)=[C@AL1]=CO").unwrap();
        let mut skeletal_coords = vec![(0.0, 0.0); graph.atoms.len()];
        let mut skeletal_stereo = HashMap::new();
        let skeletal_depicted =
            apply_extended_stereo_depictions(&graph, &mut skeletal_coords, &mut skeletal_stereo);

        assert_eq!(skeletal_depicted, HashSet::from([3]));
        assert_eq!(graph.atoms[4].hydrogens, 1);
        assert_eq!(skeletal_stereo.get(&(4, 5)), Some(&(1, true)));

        let full_graph = expand_smiles_graph_hydrogens(&graph);
        let mut full_coords = vec![(0.0, 0.0); full_graph.atoms.len()];
        let mut full_stereo = HashMap::new();
        let full_depicted =
            apply_extended_stereo_depictions(&full_graph, &mut full_coords, &mut full_stereo);
        let terminal_hydrogen = full_graph.neighbor_order[4]
            .iter()
            .copied()
            .find(|&neighbor| full_graph.atoms[neighbor].element == "H")
            .unwrap();

        assert_eq!(full_depicted, HashSet::from([3]));
        assert_eq!(full_graph.atoms[terminal_hydrogen].element, "H");
        assert_eq!(full_stereo.get(&(4, terminal_hydrogen)), Some(&(6, true)));
        assert_eq!(full_stereo.get(&(4, 5)), Some(&(1, true)));
    }

    #[test]
    fn unsupported_extended_topology_keeps_its_annotation_fallback() {
        let graph = parse_smiles_graph("[C@AL1](F)(Cl)C").unwrap();
        let mut coords = vec![(0.0, 0.0); graph.atoms.len()];
        let mut stereo_map = HashMap::new();

        assert!(apply_extended_stereo_depictions(&graph, &mut coords, &mut stereo_map).is_empty());
        assert_eq!(graph.atoms[0].stereo_annotation.as_deref(), Some("@AL1"));

        let cyclic = parse_smiles_graph("[Pt@SP1]1(F)(Cl)CCC1").unwrap();
        let mut cyclic_coords = (0..cyclic.atoms.len())
            .map(|atom| {
                let angle = atom as f32;
                (angle.cos(), angle.sin())
            })
            .collect::<Vec<_>>();
        let mut cyclic_stereo_map = HashMap::new();

        assert!(apply_extended_stereo_depictions(
            &cyclic,
            &mut cyclic_coords,
            &mut cyclic_stereo_map,
        )
        .is_empty());
        assert_eq!(cyclic.atoms[0].stereo_annotation.as_deref(), Some("@SP1"));
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
        assert!(matches!(
            al.atoms[0].extended_stereo,
            Some(ExtendedAtomStereoSpec::Allene { .. })
        ));
        assert!(matches!(
            sp.atoms[0].extended_stereo,
            Some(ExtendedAtomStereoSpec::SquarePlanar { .. })
        ));
        assert!(matches!(
            tb.atoms[0].extended_stereo,
            Some(ExtendedAtomStereoSpec::TrigonalBipyramidal { .. })
        ));
        assert!(matches!(
            oh.atoms[0].extended_stereo,
            Some(ExtendedAtomStereoSpec::Octahedral { .. })
        ));
    }
}
