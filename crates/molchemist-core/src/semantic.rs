use sdfrust::{BondOrder, BondStereo, Molecule, SGroupType, SdfFormat, StereoGroupType};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};

pub const CHEMICAL_RECORD_SCHEMA_VERSION: u16 = 1;

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum ChemicalFormat {
    MolfileV2000,
    MolfileV3000,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct ChemicalRecord {
    pub schema_version: u16,
    pub record: usize,
    pub format: ChemicalFormat,
    pub name: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub program_line: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub comment: Option<String>,
    pub atoms: Vec<ChemicalAtom>,
    pub bonds: Vec<ChemicalBond>,
    pub stereo_groups: Vec<ChemicalStereoGroup>,
    pub sgroups: Vec<ChemicalSGroup>,
    #[serde(default)]
    pub link_nodes: Vec<ChemicalLinkNode>,
    pub collections: Vec<ChemicalCollection>,
    pub properties: Vec<SdfProperty>,
    pub diagnostics: Vec<FidelityDiagnostic>,
    /// Exact selected record, excluding the `$$$$` separator.
    pub raw_record: String,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct ChemicalAtom {
    /// Zero-based index used by the current depiction anchors (`a0`, `a1`, ...).
    pub index: usize,
    /// Original V3000 atom ID, or the one-based source position for V2000.
    pub source_id: u32,
    pub element: String,
    pub coordinates: [f64; 3],
    pub formal_charge: i8,
    pub mass_difference: i8,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub isotope: Option<u16>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub stereo_parity: Option<u8>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub hydrogen_count: Option<u8>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub valence: Option<u8>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub atom_map: Option<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub rgroup_label: Option<u8>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub radical: Option<u8>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub query: Option<ChemicalAtomQuery>,
    pub rgroup_labels: Vec<u8>,
    pub attachment_points: Vec<i8>,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct ChemicalAtomQuery {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub elements: Option<Vec<String>>,
    pub is_not_list: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub hydrogen_count: Option<i8>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub valence: Option<i8>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub substitution_count: Option<i8>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub unsaturated: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ring_bond_count: Option<i8>,
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum ChemicalBondOrder {
    Single,
    Double,
    Triple,
    Aromatic,
    SingleOrDouble,
    SingleOrAromatic,
    DoubleOrAromatic,
    Any,
    Coordination,
    Hydrogen,
}

impl From<BondOrder> for ChemicalBondOrder {
    fn from(order: BondOrder) -> Self {
        match order {
            BondOrder::Single => Self::Single,
            BondOrder::Double => Self::Double,
            BondOrder::Triple => Self::Triple,
            BondOrder::Aromatic => Self::Aromatic,
            BondOrder::SingleOrDouble => Self::SingleOrDouble,
            BondOrder::SingleOrAromatic => Self::SingleOrAromatic,
            BondOrder::DoubleOrAromatic => Self::DoubleOrAromatic,
            BondOrder::Any => Self::Any,
            BondOrder::Coordination => Self::Coordination,
            BondOrder::Hydrogen => Self::Hydrogen,
        }
    }
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum ChemicalBondStereo {
    None,
    Up,
    Either,
    Down,
}

impl From<BondStereo> for ChemicalBondStereo {
    fn from(stereo: BondStereo) -> Self {
        match stereo {
            BondStereo::None => Self::None,
            BondStereo::Up => Self::Up,
            BondStereo::Either => Self::Either,
            BondStereo::Down => Self::Down,
        }
    }
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct ChemicalBond {
    /// Zero-based index used by the current depiction anchors (`b0`, `b1`, ...).
    pub index: usize,
    /// Original V3000 bond ID, or the one-based source position for V2000.
    pub source_id: u32,
    pub atom1_index: usize,
    pub atom2_index: usize,
    pub atom1_source_id: u32,
    pub atom2_source_id: u32,
    pub order: ChemicalBondOrder,
    pub stereo: ChemicalBondStereo,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub topology: Option<u8>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub reacting_center: Option<u8>,
    /// V3000 variable-attachment endpoints (`ENDPTS`).
    #[serde(default)]
    pub endpoint_source_ids: Vec<u32>,
    /// V3000 variable-attachment mode (`ATTACH`), normally `ANY` or `ALL`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub attachment_mode: Option<String>,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct ChemicalStereoGroup {
    pub kind: String,
    pub group_number: u32,
    pub atom_source_ids: Vec<u32>,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct ChemicalSGroup {
    pub id: u32,
    /// Original three-letter CTfile SGroup type code.
    pub type_code: String,
    pub kind: String,
    pub atom_source_ids: Vec<u32>,
    pub bond_source_ids: Vec<u32>,
    pub crossing_bond_source_ids: Vec<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub subscript: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub superscript: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub parent_id: Option<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub connectivity: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub subtype: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub bracket_type: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub class: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub component_number: Option<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub expanded: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub field_name: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub field_value: Option<String>,
    pub field_values: Vec<String>,
    pub brackets: Vec<[f64; 4]>,
    #[serde(default)]
    pub attachment_points: Vec<ChemicalSGroupAttachmentPoint>,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct ChemicalSGroupAttachmentPoint {
    pub atom_source_id: u32,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub leaving_atom_source_id: Option<u32>,
    pub id: String,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct ChemicalLinkNodeConnection {
    pub atom_source_id: u32,
    pub neighbor_source_id: u32,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct ChemicalLinkNode {
    pub min_repeat: u32,
    pub max_repeat: u32,
    pub connections: Vec<ChemicalLinkNodeConnection>,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct ChemicalCollection {
    pub id: u32,
    /// Complete CTfile COLLECTION name, including its namespace.
    pub name: String,
    pub kind: String,
    pub atom_source_ids: Vec<u32>,
    pub bond_source_ids: Vec<u32>,
    #[serde(default)]
    pub sgroup_source_ids: Vec<u32>,
    #[serde(default)]
    pub object3d_source_ids: Vec<u32>,
    #[serde(default)]
    pub rgroup_ids: Vec<String>,
    /// Uninterpreted `MEMBERS=(...)` values, in source order.
    #[serde(default)]
    pub members: Vec<String>,
    /// Exact logical V3000 COLLECTION entry, without the `M  V30 ` prefix.
    pub raw_entry: String,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct SdfProperty {
    pub name: String,
    pub value: String,
    /// Original data-header line, including optional registry information.
    pub header: String,
    /// One-based source line containing `header`.
    pub line: usize,
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum DiagnosticSeverity {
    Warning,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct FidelityDiagnostic {
    pub code: String,
    pub severity: DiagnosticSeverity,
    pub feature: String,
    pub message: String,
    pub atom_source_ids: Vec<u32>,
    pub bond_source_ids: Vec<u32>,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum FidelityPolicy {
    Ignore,
    #[default]
    Warn,
    Strict,
}

impl FidelityPolicy {
    pub fn parse(value: &str) -> Result<Self, String> {
        match value {
            "ignore" => Ok(Self::Ignore),
            "warn" => Ok(Self::Warn),
            "strict" => Ok(Self::Strict),
            _ => Err(format!(
                "Invalid fidelity policy {value:?}; expected ignore, warn, or strict"
            )),
        }
    }

    pub fn enforce(self, diagnostics: &[FidelityDiagnostic]) -> Result<(), String> {
        if self != Self::Strict || diagnostics.is_empty() {
            return Ok(());
        }
        Err(format!(
            "faithful depiction is not available: {}",
            diagnostics
                .iter()
                .map(|diagnostic| diagnostic.message.as_str())
                .collect::<Vec<_>>()
                .join("; ")
        ))
    }
}

#[derive(Clone, Debug, Default)]
pub(crate) struct SdfAtomMetadata {
    pub isotope: Option<u16>,
    pub radical: Option<u8>,
    pub atom_map: Option<u32>,
}

impl ChemicalRecord {
    pub(crate) fn from_sdf(
        molecule: &Molecule,
        raw_record: &str,
        record: usize,
        atom_metadata: &[SdfAtomMetadata],
    ) -> Self {
        let atom_source_ids = molecule
            .atoms
            .iter()
            .enumerate()
            .map(|(index, atom)| atom.v3000_id.unwrap_or((index + 1) as u32))
            .collect::<Vec<_>>();
        let bond_source_ids = molecule
            .bonds
            .iter()
            .enumerate()
            .map(|(index, bond)| bond.v3000_id.unwrap_or((index + 1) as u32))
            .collect::<Vec<_>>();
        let ctfile_atom_metadata =
            raw_ctfile_atom_metadata(raw_record, molecule.format_version, &atom_source_ids);
        let ctfile_bond_metadata = raw_v3000_bond_metadata(raw_record, &bond_source_ids);

        let atoms = molecule
            .atoms
            .iter()
            .enumerate()
            .map(|(index, atom)| {
                let metadata = atom_metadata.get(index).cloned().unwrap_or_default();
                let ctfile = ctfile_atom_metadata
                    .get(&atom_source_ids[index])
                    .cloned()
                    .unwrap_or_default();
                let hydrogen_count = match molecule.format_version {
                    SdfFormat::V2000 => ctfile.hydrogen_count.or_else(|| {
                        atom.hydrogen_count
                            .map(|value| value.saturating_sub(1) as i8)
                    }),
                    SdfFormat::V3000 => ctfile
                        .hydrogen_count
                        .or_else(|| atom.hydrogen_count.map(|value| value as i8)),
                };
                let valence = ctfile
                    .valence
                    .or_else(|| atom.valence.map(|value| value as i8));
                let query = (ctfile.elements.is_some()
                    || hydrogen_count.is_some()
                    || valence.is_some()
                    || ctfile.substitution_count.is_some()
                    || ctfile.unsaturated.is_some()
                    || ctfile.ring_bond_count.is_some())
                .then_some(ChemicalAtomQuery {
                    elements: ctfile.elements,
                    is_not_list: ctfile.is_not_list,
                    hydrogen_count,
                    valence,
                    substitution_count: ctfile.substitution_count,
                    unsaturated: ctfile.unsaturated,
                    ring_bond_count: ctfile.ring_bond_count,
                });
                let rgroup_labels = if ctfile.rgroup_labels.is_empty() {
                    atom.rgroup_label.into_iter().collect()
                } else {
                    ctfile.rgroup_labels
                };
                ChemicalAtom {
                    index,
                    source_id: atom_source_ids[index],
                    element: atom.element.clone(),
                    coordinates: [atom.x, atom.y, atom.z],
                    formal_charge: atom.formal_charge,
                    mass_difference: atom.mass_difference,
                    isotope: metadata.isotope,
                    stereo_parity: atom.stereo_parity,
                    hydrogen_count: hydrogen_count.and_then(|value| u8::try_from(value).ok()),
                    valence: atom.valence,
                    atom_map: atom.atom_atom_mapping.or(metadata.atom_map),
                    rgroup_label: rgroup_labels.first().copied().or(atom.rgroup_label),
                    radical: atom.radical.or(metadata.radical),
                    query,
                    rgroup_labels,
                    attachment_points: ctfile.attachment_points,
                }
            })
            .collect::<Vec<_>>();

        let bonds = molecule
            .bonds
            .iter()
            .enumerate()
            .map(|(index, bond)| {
                let metadata = ctfile_bond_metadata
                    .get(&bond_source_ids[index])
                    .cloned()
                    .unwrap_or_default();
                ChemicalBond {
                    index,
                    source_id: bond_source_ids[index],
                    atom1_index: bond.atom1,
                    atom2_index: bond.atom2,
                    atom1_source_id: atom_source_ids[bond.atom1],
                    atom2_source_id: atom_source_ids[bond.atom2],
                    order: bond.order.into(),
                    stereo: bond.stereo.into(),
                    topology: bond.topology,
                    reacting_center: bond.reacting_center,
                    endpoint_source_ids: metadata.endpoint_source_ids,
                    attachment_mode: metadata.attachment_mode,
                }
            })
            .collect::<Vec<_>>();

        let mut stereo_groups = molecule
            .stereogroups
            .iter()
            .map(|group| ChemicalStereoGroup {
                kind: stereo_group_kind(group.group_type).to_string(),
                group_number: group.group_number,
                atom_source_ids: source_ids(&group.atoms, &atom_source_ids),
            })
            .collect::<Vec<_>>();
        for group in raw_v3000_stereo_groups(raw_record, &atom_source_ids) {
            if !stereo_groups.contains(&group) {
                stereo_groups.push(group);
            }
        }

        let mut sgroups = raw_ctfile_sgroups(
            raw_record,
            molecule.format_version,
            &atom_source_ids,
            &bond_source_ids,
        );
        if sgroups.is_empty() {
            sgroups = molecule
                .sgroups
                .iter()
                .map(|group| ChemicalSGroup {
                    id: group.id,
                    type_code: sgroup_code(group.sgroup_type).to_string(),
                    kind: sgroup_kind(group.sgroup_type).to_string(),
                    atom_source_ids: source_ids(&group.atoms, &atom_source_ids),
                    bond_source_ids: source_ids(&group.bonds, &bond_source_ids),
                    crossing_bond_source_ids: source_ids(&group.crossing_bonds, &bond_source_ids),
                    label: group.label.clone(),
                    subscript: group.subscript.clone(),
                    superscript: group.superscript.clone(),
                    parent_id: group.parent_id,
                    connectivity: group.connectivity.clone(),
                    subtype: None,
                    bracket_type: None,
                    class: None,
                    component_number: None,
                    expanded: None,
                    field_name: group.field_name.clone(),
                    field_value: group.field_value.clone(),
                    field_values: group.field_value.clone().into_iter().collect(),
                    brackets: group
                        .brackets
                        .iter()
                        .map(|&(x1, y1, x2, y2)| [x1, y1, x2, y2])
                        .collect(),
                    attachment_points: Vec::new(),
                })
                .collect();
        }

        // Parse COLLECTION entries from the source rather than sdfrust's
        // compatibility aliases. Atom lists, R-group labels, and attachment
        // points are atom/SGROUP attributes in standard CTfile syntax and must
        // not be reconstructed from private MDLV30 collection names.
        let collections = raw_v3000_collections(raw_record);
        let link_nodes = raw_ctfile_link_nodes(raw_record, molecule.format_version);

        let properties = ordered_sdf_properties(raw_record);
        let diagnostics = fidelity_diagnostics(
            molecule.format_version,
            raw_record,
            &atoms,
            &bonds,
            &sgroups,
            &link_nodes,
            &collections,
        );

        Self {
            schema_version: CHEMICAL_RECORD_SCHEMA_VERSION,
            record,
            format: match molecule.format_version {
                SdfFormat::V2000 => ChemicalFormat::MolfileV2000,
                SdfFormat::V3000 => ChemicalFormat::MolfileV3000,
            },
            name: molecule.name.clone(),
            program_line: molecule.program_line.clone(),
            comment: molecule.comment.clone(),
            atoms,
            bonds,
            stereo_groups,
            sgroups,
            link_nodes,
            collections,
            properties,
            diagnostics,
            raw_record: raw_record.to_string(),
        }
    }
}

fn source_ids(indices: &[usize], ids: &[u32]) -> Vec<u32> {
    indices
        .iter()
        .filter_map(|&index| ids.get(index).copied())
        .collect()
}

fn stereo_group_kind(kind: StereoGroupType) -> &'static str {
    match kind {
        StereoGroupType::Absolute => "absolute",
        StereoGroupType::Or => "or",
        StereoGroupType::And => "and",
    }
}

fn sgroup_kind(kind: SGroupType) -> &'static str {
    match kind {
        SGroupType::Superatom => "superatom",
        SGroupType::Multiple => "multiple",
        SGroupType::StructureRepeatUnit => "structure-repeat-unit",
        SGroupType::Data => "data",
        SGroupType::Generic => "generic",
        SGroupType::Monomer => "monomer",
        SGroupType::Mer => "mer",
        SGroupType::Copolymer => "copolymer",
        SGroupType::Component => "component",
        SGroupType::Mixture => "mixture",
        SGroupType::Formulation => "formulation",
    }
}

fn sgroup_code(kind: SGroupType) -> &'static str {
    match kind {
        SGroupType::Superatom => "SUP",
        SGroupType::Multiple => "MUL",
        SGroupType::StructureRepeatUnit => "SRU",
        SGroupType::Data => "DAT",
        SGroupType::Generic => "GEN",
        SGroupType::Monomer => "MON",
        SGroupType::Mer => "MER",
        SGroupType::Copolymer => "COP",
        SGroupType::Component => "COM",
        SGroupType::Mixture => "MIX",
        SGroupType::Formulation => "FOR",
    }
}

#[derive(Clone, Debug, Default)]
struct CtfileAtomMetadata {
    elements: Option<Vec<String>>,
    is_not_list: bool,
    hydrogen_count: Option<i8>,
    valence: Option<i8>,
    substitution_count: Option<i8>,
    unsaturated: Option<bool>,
    ring_bond_count: Option<i8>,
    rgroup_labels: Vec<u8>,
    attachment_points: Vec<i8>,
}

#[derive(Clone, Debug, Default)]
struct CtfileBondMetadata {
    endpoint_source_ids: Vec<u32>,
    attachment_mode: Option<String>,
}

fn raw_v3000_bond_metadata(
    raw_record: &str,
    bond_source_ids: &[u32],
) -> HashMap<u32, CtfileBondMetadata> {
    let valid_bonds = bond_source_ids.iter().copied().collect::<HashSet<_>>();
    let mut metadata = HashMap::new();
    let mut in_bond_block = false;
    for line in v3000_logical_lines(raw_record) {
        if line == "BEGIN BOND" {
            in_bond_block = true;
            continue;
        }
        if line == "END BOND" {
            break;
        }
        if !in_bond_block {
            continue;
        }
        let Some(source_id) = line
            .split_whitespace()
            .next()
            .and_then(|value| value.parse::<u32>().ok())
            .filter(|id| valid_bonds.contains(id))
        else {
            continue;
        };
        let endpoint_source_ids = v3000_id_list(&line, "ENDPTS");
        let attachment_mode = v3000_value(&line, "ATTACH")
            .map(|value| value.to_ascii_uppercase())
            .filter(|value| !value.is_empty());
        if !endpoint_source_ids.is_empty() || attachment_mode.is_some() {
            metadata.insert(
                source_id,
                CtfileBondMetadata {
                    endpoint_source_ids,
                    attachment_mode,
                },
            );
        }
    }
    metadata
}

fn raw_ctfile_atom_metadata(
    raw_record: &str,
    format: SdfFormat,
    atom_source_ids: &[u32],
) -> HashMap<u32, CtfileAtomMetadata> {
    match format {
        SdfFormat::V2000 => raw_v2000_atom_metadata(raw_record, atom_source_ids),
        SdfFormat::V3000 => raw_v3000_atom_metadata(raw_record, atom_source_ids),
    }
}

fn raw_v3000_atom_metadata(
    raw_record: &str,
    atom_source_ids: &[u32],
) -> HashMap<u32, CtfileAtomMetadata> {
    let valid_ids = atom_source_ids.iter().copied().collect::<HashSet<_>>();
    let mut metadata = HashMap::new();
    let mut in_atom_block = false;

    for line in v3000_logical_lines(raw_record) {
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

        let Some(source_id) = line
            .split_whitespace()
            .next()
            .and_then(|value| value.parse::<u32>().ok())
        else {
            continue;
        };
        if !valid_ids.contains(&source_id) {
            continue;
        }

        let mut atom = CtfileAtomMetadata::default();
        if let Some((elements, is_not_list)) = v3000_atom_list(&line) {
            atom.elements = Some(elements);
            atom.is_not_list = is_not_list;
        }
        atom.hydrogen_count = v3000_hydrogen_count(&line);
        atom.valence = v3000_integer(&line, "VAL").filter(|value| *value != 0);
        atom.substitution_count =
            v3000_integer(&line, "SUBST").and_then(normalize_substitution_count);
        atom.unsaturated = v3000_integer(&line, "UNSAT")
            .filter(|value| *value != 0)
            .map(|_| true);
        atom.ring_bond_count = v3000_integer(&line, "RBCNT").and_then(normalize_ring_bond_count);
        atom.attachment_points = v3000_integer(&line, "ATTCHPT").into_iter().collect();
        atom.rgroup_labels = v3000_id_list(&line, "RGROUPS")
            .into_iter()
            .filter_map(|value| u8::try_from(value).ok())
            .collect();
        metadata.insert(source_id, atom);
    }
    metadata
}

fn v3000_hydrogen_count(line: &str) -> Option<i8> {
    match v3000_integer(line, "HCOUNT")? {
        0 => None,
        // V3000 stores the depicted query count directly: -1 is H0, while
        // positive values 1..=4 are H1..H4.
        -1 => Some(0),
        value @ 1..=4 => Some(value),
        // Preserve an invalid value so strict diagnostics can reject it.
        value => Some(value),
    }
}

fn normalize_substitution_count(value: i8) -> Option<i8> {
    match value {
        0 => None,
        value if value >= 6 => Some(6),
        value => Some(value),
    }
}

fn normalize_ring_bond_count(value: i8) -> Option<i8> {
    match value {
        0 => None,
        value if value >= 4 => Some(4),
        value => Some(value),
    }
}

fn v3000_atom_list(line: &str) -> Option<(Vec<String>, bool)> {
    let after_id = line
        .trim_start()
        .split_once(char::is_whitespace)?
        .1
        .trim_start();
    let after_id = after_id.strip_prefix('"').unwrap_or(after_id);
    let (list, is_not_list) = if let Some(list) = after_id.strip_prefix("NOT ") {
        (list, true)
    } else {
        (after_id, false)
    };
    let list = list.strip_prefix('[')?;
    let end = list.find(']')?;
    let elements = list[..end]
        .split(',')
        .map(str::trim)
        .filter(|element| !element.is_empty())
        .map(ToString::to_string)
        .collect::<Vec<_>>();
    (!elements.is_empty()).then_some((elements, is_not_list))
}

fn v3000_integer(line: &str, key: &str) -> Option<i8> {
    v3000_value(line, key)?.parse().ok()
}

fn raw_v2000_atom_metadata(
    raw_record: &str,
    atom_source_ids: &[u32],
) -> HashMap<u32, CtfileAtomMetadata> {
    let valid_ids = atom_source_ids.iter().copied().collect::<HashSet<_>>();
    let mut metadata = HashMap::<u32, CtfileAtomMetadata>::new();

    for line in raw_record.lines() {
        if line.starts_with("M  ALS") {
            let values = line.split_whitespace().collect::<Vec<_>>();
            let Some(source_id) = values.get(2).and_then(|value| value.parse::<u32>().ok()) else {
                continue;
            };
            if !valid_ids.contains(&source_id) {
                continue;
            }
            let count = values
                .get(3)
                .and_then(|value| value.parse::<usize>().ok())
                .unwrap_or(0);
            let is_not_list = values.get(4).is_some_and(|value| *value == "T");
            let elements = values
                .iter()
                .skip(5)
                .take(count)
                .map(|value| (*value).to_string())
                .collect::<Vec<_>>();
            let atom = metadata.entry(source_id).or_default();
            atom.elements = Some(elements);
            atom.is_not_list = is_not_list;
            continue;
        }

        let property = if line.starts_with("M  SUB") {
            Some("SUB")
        } else if line.starts_with("M  UNS") {
            Some("UNS")
        } else if line.starts_with("M  RBC") {
            Some("RBC")
        } else if line.starts_with("M  RGP") {
            Some("RGP")
        } else if line.starts_with("M  APO") {
            Some("APO")
        } else {
            None
        };
        let Some(property) = property else {
            continue;
        };
        let values = line.split_whitespace().skip(3).collect::<Vec<_>>();
        for pair in values.chunks_exact(2) {
            let Some(source_id) = pair[0].parse::<u32>().ok() else {
                continue;
            };
            let Some(value) = pair[1].parse::<i8>().ok() else {
                continue;
            };
            if !valid_ids.contains(&source_id) {
                continue;
            }
            let atom = metadata.entry(source_id).or_default();
            match property {
                "SUB" => atom.substitution_count = normalize_substitution_count(value),
                "UNS" if value != 0 => atom.unsaturated = Some(true),
                "RBC" => atom.ring_bond_count = normalize_ring_bond_count(value),
                "RGP" => {
                    if let Ok(label) = u8::try_from(value) {
                        atom.rgroup_labels.push(label);
                    }
                }
                "APO" if value != 0 => atom.attachment_points.push(value),
                "UNS" | "APO" => {}
                _ => unreachable!(),
            }
        }
    }
    metadata
}

fn raw_ctfile_sgroups(
    raw_record: &str,
    format: SdfFormat,
    atom_source_ids: &[u32],
    bond_source_ids: &[u32],
) -> Vec<ChemicalSGroup> {
    match format {
        SdfFormat::V2000 => raw_v2000_sgroups(raw_record, atom_source_ids, bond_source_ids),
        SdfFormat::V3000 => raw_v3000_sgroups(raw_record, atom_source_ids, bond_source_ids),
    }
}

fn empty_sgroup(id: u32, type_code: String) -> ChemicalSGroup {
    let kind = sgroup_kind_from_code(&type_code)
        .unwrap_or("unknown")
        .to_string();
    ChemicalSGroup {
        id,
        type_code,
        kind,
        atom_source_ids: Vec::new(),
        bond_source_ids: Vec::new(),
        crossing_bond_source_ids: Vec::new(),
        label: None,
        subscript: None,
        superscript: None,
        parent_id: None,
        connectivity: None,
        subtype: None,
        bracket_type: None,
        class: None,
        component_number: None,
        expanded: None,
        field_name: None,
        field_value: None,
        field_values: Vec::new(),
        brackets: Vec::new(),
        attachment_points: Vec::new(),
    }
}

fn raw_v3000_sgroups(
    raw_record: &str,
    atom_source_ids: &[u32],
    bond_source_ids: &[u32],
) -> Vec<ChemicalSGroup> {
    let valid_atoms = atom_source_ids.iter().copied().collect::<HashSet<_>>();
    let valid_bonds = bond_source_ids.iter().copied().collect::<HashSet<_>>();
    let mut groups = Vec::new();
    let mut in_sgroup_block = false;

    for line in v3000_logical_lines(raw_record) {
        if line == "BEGIN SGROUP" {
            in_sgroup_block = true;
            continue;
        }
        if line == "END SGROUP" {
            break;
        }
        if !in_sgroup_block || line.starts_with("DEFAULT") {
            continue;
        }
        let values = line.split_whitespace().collect::<Vec<_>>();
        let (Some(id), Some(type_code)) = (
            values.first().and_then(|value| value.parse::<u32>().ok()),
            values.get(1),
        ) else {
            continue;
        };
        let mut group = empty_sgroup(id, type_code.to_ascii_uppercase());
        group.atom_source_ids = v3000_id_list(&line, "ATOMS")
            .into_iter()
            .filter(|id| valid_atoms.contains(id))
            .collect();
        group.bond_source_ids = v3000_id_list(&line, "CBONDS")
            .into_iter()
            .filter(|id| valid_bonds.contains(id))
            .collect();
        group.crossing_bond_source_ids = v3000_id_list(&line, "XBONDS")
            .into_iter()
            .filter(|id| valid_bonds.contains(id))
            .collect();
        group.label = v3000_value(&line, "LABEL");
        group.subscript = group.label.clone();
        group.parent_id = v3000_value(&line, "PARENT").and_then(|value| value.parse().ok());
        group.connectivity = v3000_value(&line, "CONNECT");
        group.subtype = v3000_value(&line, "SUBTYPE");
        group.bracket_type = v3000_value(&line, "BRKTYP");
        group.class = v3000_value(&line, "CLASS");
        group.component_number = v3000_value(&line, "COMPNO").and_then(|value| value.parse().ok());
        group.expanded = v3000_value(&line, "ESTATE").map(|value| {
            matches!(
                value.to_ascii_uppercase().as_str(),
                "E" | "EXP" | "EXPANDED"
            )
        });
        group.field_name = v3000_value(&line, "FIELDNAME");
        group.field_values = v3000_values(&line, "FIELDDATA");
        group.field_value = group.field_values.first().cloned();
        group.brackets = v3000_brackets(&line);
        group.attachment_points = v3000_sgroup_attachment_points(&line, &valid_atoms);
        groups.push(group);
    }
    groups
}

fn raw_v2000_sgroups(
    raw_record: &str,
    atom_source_ids: &[u32],
    bond_source_ids: &[u32],
) -> Vec<ChemicalSGroup> {
    let valid_atoms = atom_source_ids.iter().copied().collect::<HashSet<_>>();
    let valid_bonds = bond_source_ids.iter().copied().collect::<HashSet<_>>();
    let mut groups = HashMap::<u32, ChemicalSGroup>::new();
    let mut order = Vec::new();

    for line in raw_record.lines() {
        if line.starts_with("M  STY") {
            let values = line.split_whitespace().skip(3).collect::<Vec<_>>();
            for pair in values.chunks_exact(2) {
                let Some(id) = pair[0].parse::<u32>().ok() else {
                    continue;
                };
                if !groups.contains_key(&id) {
                    order.push(id);
                }
                groups.insert(id, empty_sgroup(id, pair[1].to_ascii_uppercase()));
            }
            continue;
        }
        let values = line.split_whitespace().collect::<Vec<_>>();
        if line.starts_with("M  SCN")
            || line.starts_with("M  SST")
            || line.starts_with("M  SPL")
            || line.starts_with("M  SNC")
            || line.starts_with("M  SBT")
        {
            for pair in values.get(3..).unwrap_or_default().chunks_exact(2) {
                let Some(group) = pair[0]
                    .parse::<u32>()
                    .ok()
                    .and_then(|id| groups.get_mut(&id))
                else {
                    continue;
                };
                if line.starts_with("M  SCN") {
                    group.connectivity = Some(pair[1].to_string());
                } else if line.starts_with("M  SST") {
                    group.subtype = Some(pair[1].to_string());
                } else if line.starts_with("M  SPL") {
                    group.parent_id = pair[1].parse().ok();
                } else if line.starts_with("M  SNC") {
                    group.component_number = pair[1].parse().ok();
                } else {
                    group.bracket_type = Some(pair[1].to_string());
                }
            }
            continue;
        }
        if line.starts_with("M  SDS") && values.get(2) == Some(&"EXP") {
            for id in values.iter().skip(4).filter_map(|value| value.parse().ok()) {
                if let Some(group) = groups.get_mut(&id) {
                    group.expanded = Some(true);
                }
            }
            continue;
        }
        let group_id = values.get(2).and_then(|value| value.parse::<u32>().ok());
        let Some(group_id) = group_id else {
            continue;
        };
        let Some(group) = groups.get_mut(&group_id) else {
            continue;
        };

        if line.starts_with("M  SAL") {
            group.atom_source_ids.extend(
                values
                    .iter()
                    .skip(4)
                    .filter_map(|value| value.parse::<u32>().ok())
                    .filter(|id| valid_atoms.contains(id)),
            );
        } else if line.starts_with("M  SBL") {
            let bonds = values
                .iter()
                .skip(4)
                .filter_map(|value| value.parse::<u32>().ok())
                .filter(|id| valid_bonds.contains(id));
            if group.kind == "data" {
                group.bond_source_ids.extend(bonds);
            } else {
                group.crossing_bond_source_ids.extend(bonds);
            }
        } else if line.starts_with("M  SMT") {
            let label = values.iter().skip(3).copied().collect::<Vec<_>>().join(" ");
            group.label = (!label.is_empty()).then_some(label);
            group.subscript = group.label.clone();
        } else if line.starts_with("M  SCL") {
            let class = values.iter().skip(3).copied().collect::<Vec<_>>().join(" ");
            group.class = (!class.is_empty()).then_some(class);
        } else if line.starts_with("M  SAP") {
            group
                .attachment_points
                .extend(v2000_sgroup_attachment_points(line, &valid_atoms));
        } else if line.starts_with("M  SDI") {
            let coordinates = values
                .iter()
                .skip(4)
                .filter_map(|value| value.parse::<f64>().ok())
                .collect::<Vec<_>>();
            if coordinates.len() >= 4 {
                group.brackets.push([
                    coordinates[0],
                    coordinates[1],
                    coordinates[2],
                    coordinates[3],
                ]);
            }
        } else if line.starts_with("M  SDT") {
            let field = line.get(10..40).unwrap_or_default().trim();
            if !field.is_empty() {
                group.field_name = Some(field.to_string());
            }
        } else if line.starts_with("M  SCD") || line.starts_with("M  SED") {
            let value = line.get(10..).unwrap_or_default().trim_end();
            if let Some(previous) = group.field_values.last_mut() {
                if line.starts_with("M  SCD") {
                    previous.push_str(value);
                } else {
                    previous.push_str(value);
                    group.field_value = Some(previous.clone());
                }
            } else {
                group.field_values.push(value.to_string());
                if line.starts_with("M  SED") {
                    group.field_value = Some(value.to_string());
                }
            }
        }
    }
    order
        .into_iter()
        .filter_map(|id| groups.remove(&id))
        .collect()
}

fn v2000_sgroup_attachment_points(
    line: &str,
    valid_atoms: &HashSet<u32>,
) -> Vec<ChemicalSGroupAttachmentPoint> {
    let values = line.split_whitespace().collect::<Vec<_>>();
    let count = values
        .get(3)
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(0);
    values
        .get(4..)
        .unwrap_or_default()
        .chunks(3)
        .take(count)
        .filter_map(|entry| {
            let atom_source_id = entry.first()?.parse::<u32>().ok()?;
            if !valid_atoms.contains(&atom_source_id) {
                return None;
            }
            let leaving_atom_source_id = entry
                .get(1)
                .and_then(|value| value.parse::<u32>().ok())
                .filter(|id| *id != 0 && valid_atoms.contains(id));
            Some(ChemicalSGroupAttachmentPoint {
                atom_source_id,
                leaving_atom_source_id,
                id: entry.get(2).copied().unwrap_or("").trim().to_string(),
            })
        })
        .collect()
}

fn v3000_sgroup_attachment_points(
    line: &str,
    valid_atoms: &HashSet<u32>,
) -> Vec<ChemicalSGroupAttachmentPoint> {
    v3000_parenthesized_values(line, "SAP")
        .into_iter()
        .filter_map(|values| {
            if values.len() < 4 || values[0] != "3" {
                return None;
            }
            let atom_source_id = values[1].parse::<u32>().ok()?;
            if !valid_atoms.contains(&atom_source_id) {
                return None;
            }
            let leaving_atom_source_id = if values[2].eq_ignore_ascii_case("AIDX") {
                Some(atom_source_id)
            } else {
                values[2]
                    .parse::<u32>()
                    .ok()
                    .filter(|id| *id != 0 && valid_atoms.contains(id))
            };
            Some(ChemicalSGroupAttachmentPoint {
                atom_source_id,
                leaving_atom_source_id,
                id: values[3].trim_matches('"').to_string(),
            })
        })
        .collect()
}

fn sgroup_kind_from_code(value: &str) -> Option<&'static str> {
    match value.to_ascii_uppercase().as_str() {
        "SUP" => Some("superatom"),
        "MUL" => Some("multiple"),
        "SRU" => Some("structure-repeat-unit"),
        "DAT" => Some("data"),
        "GEN" => Some("generic"),
        "MON" => Some("monomer"),
        "MER" => Some("mer"),
        "COP" => Some("copolymer"),
        "COM" => Some("component"),
        "MIX" => Some("mixture"),
        "FOR" => Some("formulation"),
        "CRO" => Some("crosslink"),
        "MOD" => Some("modification"),
        "GRA" => Some("graft"),
        "ANY" => Some("any"),
        _ => None,
    }
}

fn raw_v3000_stereo_groups(raw_record: &str, atom_source_ids: &[u32]) -> Vec<ChemicalStereoGroup> {
    let valid_ids = atom_source_ids.iter().copied().collect::<HashSet<_>>();
    let mut groups = Vec::new();
    for line in v3000_logical_lines(raw_record) {
        let Some(name) = line.split_whitespace().next() else {
            continue;
        };
        let Some(tag) = internal_collection_subname(name) else {
            continue;
        };
        let tag = tag.to_ascii_uppercase();
        let (kind, group_number) = if tag == "STEABS" {
            ("absolute", 0)
        } else if let Some(number) = tag.strip_prefix("STEREL") {
            ("or", number.parse().unwrap_or(1))
        } else if let Some(number) = tag.strip_prefix("STERAC") {
            ("and", number.parse().unwrap_or(1))
        } else {
            continue;
        };
        let atom_source_ids = v3000_id_list(&line, "ATOMS")
            .into_iter()
            .filter(|id| valid_ids.contains(id))
            .collect();
        groups.push(ChemicalStereoGroup {
            kind: kind.to_string(),
            group_number,
            atom_source_ids,
        });
    }
    groups
}

fn raw_v3000_collections(raw_record: &str) -> Vec<ChemicalCollection> {
    let mut collections = Vec::new();
    let mut in_collection = false;
    for line in v3000_logical_lines(raw_record) {
        if line == "BEGIN COLLECTION" {
            in_collection = true;
            continue;
        }
        if line == "END COLLECTION" {
            break;
        }
        if !in_collection {
            continue;
        }
        let entry = line.strip_prefix("DEFAULT ").unwrap_or(&line);
        let Some(name) = v3000_collection_name(entry) else {
            continue;
        };
        let Some(internal_name) = internal_collection_subname(&name) else {
            collections.push(ChemicalCollection {
                id: (collections.len() + 1) as u32,
                name,
                kind: "user-defined".to_string(),
                atom_source_ids: v3000_id_list(&line, "ATOMS"),
                bond_source_ids: v3000_id_list(&line, "BONDS"),
                sgroup_source_ids: v3000_id_list(&line, "SGROUPS"),
                object3d_source_ids: v3000_id_list(&line, "OBJ3DS"),
                rgroup_ids: v3000_list_values(&line, "RGROUPS"),
                members: v3000_list_values(&line, "MEMBERS"),
                raw_entry: line,
            });
            continue;
        };
        let upper = internal_name.to_ascii_uppercase();
        if upper == "STEABS" || upper.starts_with("STEREL") || upper.starts_with("STERAC") {
            // Enhanced-stereo collections are normalized into `stereo_groups`.
            continue;
        }
        collections.push(ChemicalCollection {
            id: (collections.len() + 1) as u32,
            name,
            kind: if upper == "HILITE" {
                "highlight"
            } else {
                "internal-unknown"
            }
            .to_string(),
            atom_source_ids: v3000_id_list(&line, "ATOMS"),
            bond_source_ids: v3000_id_list(&line, "BONDS"),
            sgroup_source_ids: v3000_id_list(&line, "SGROUPS"),
            object3d_source_ids: v3000_id_list(&line, "OBJ3DS"),
            rgroup_ids: v3000_list_values(&line, "RGROUPS"),
            members: v3000_list_values(&line, "MEMBERS"),
            raw_entry: line,
        });
    }
    collections
}

fn internal_collection_subname(name: &str) -> Option<&str> {
    let prefix = name.get(.."MDLV30/".len())?;
    prefix
        .eq_ignore_ascii_case("MDLV30/")
        .then(|| &name["MDLV30/".len()..])
}

fn v3000_collection_name(line: &str) -> Option<String> {
    let line = line.trim_start();
    if let Some(quoted) = line.strip_prefix('"') {
        let end = quoted.find('"')?;
        Some(quoted[..end].to_string())
    } else {
        line.split_whitespace().next().map(ToString::to_string)
    }
}

fn raw_ctfile_link_nodes(raw_record: &str, format: SdfFormat) -> Vec<ChemicalLinkNode> {
    match format {
        SdfFormat::V3000 => v3000_logical_lines(raw_record)
            .into_iter()
            .filter_map(|line| line.strip_prefix("LINKNODE ").map(parse_v3000_link_node))
            .flatten()
            .collect(),
        SdfFormat::V2000 => raw_record
            .lines()
            .filter(|line| line.starts_with("M  LIN"))
            .flat_map(parse_v2000_link_nodes)
            .collect(),
    }
}

fn parse_v3000_link_node(line: &str) -> Option<ChemicalLinkNode> {
    let values = line
        .split_whitespace()
        .filter_map(|value| value.parse::<u32>().ok())
        .collect::<Vec<_>>();
    let (&min_repeat, &max_repeat, &count) = (values.first()?, values.get(1)?, values.get(2)?);
    let connections = values
        .get(3..)?
        .chunks_exact(2)
        .take(count as usize)
        .map(|pair| ChemicalLinkNodeConnection {
            atom_source_id: pair[0],
            neighbor_source_id: pair[1],
        })
        .collect::<Vec<_>>();
    (connections.len() == count as usize).then_some(ChemicalLinkNode {
        min_repeat,
        max_repeat,
        connections,
    })
}

fn parse_v2000_link_nodes(line: &str) -> Vec<ChemicalLinkNode> {
    let count = line
        .get(6..9)
        .and_then(|value| value.trim().parse::<usize>().ok())
        .unwrap_or(0);
    let mut nodes = Vec::new();
    let mut offset = 9;
    for _ in 0..count {
        let field = |start: usize| {
            line.get(start..start + 4)
                .and_then(|value| value.trim().parse::<u32>().ok())
        };
        let (Some(atom), Some(max_repeat), Some(neighbor1)) =
            (field(offset), field(offset + 4), field(offset + 8))
        else {
            break;
        };
        let neighbor2 = field(offset + 12).filter(|value| *value != 0);
        let mut connections = vec![ChemicalLinkNodeConnection {
            atom_source_id: atom,
            neighbor_source_id: neighbor1,
        }];
        if let Some(neighbor_source_id) = neighbor2 {
            connections.push(ChemicalLinkNodeConnection {
                atom_source_id: atom,
                neighbor_source_id,
            });
        }
        nodes.push(ChemicalLinkNode {
            min_repeat: 1,
            max_repeat,
            connections,
        });
        offset += 16;
    }
    nodes
}

fn link_node_atoms_are_valid(node: &ChemicalLinkNode, valid: &HashSet<u32>) -> bool {
    !node.connections.is_empty()
        && node.connections.iter().all(|connection| {
            valid.contains(&connection.atom_source_id)
                && valid.contains(&connection.neighbor_source_id)
        })
}

fn v3000_logical_lines(raw_record: &str) -> Vec<String> {
    let mut lines = Vec::new();
    let mut current = String::new();
    for line in raw_record.lines() {
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

fn v3000_id_list(line: &str, key: &str) -> Vec<u32> {
    let marker = format!("{key}=(");
    let Some(start) = line.find(&marker) else {
        return Vec::new();
    };
    let values = &line[start + marker.len()..];
    let Some(end) = values.find(')') else {
        return Vec::new();
    };
    values[..end]
        .split_whitespace()
        .skip(1)
        .filter_map(|value| value.parse().ok())
        .collect()
}

fn v3000_list_values(line: &str, key: &str) -> Vec<String> {
    let Some(values) = v3000_parenthesized_values(line, key).into_iter().next() else {
        return Vec::new();
    };
    let Some(count) = values.first().and_then(|value| value.parse::<usize>().ok()) else {
        return Vec::new();
    };
    values
        .into_iter()
        .skip(1)
        .take(count)
        .map(|value| value.trim_matches('"').to_string())
        .filter(|value| !value.is_empty())
        .collect()
}

fn v3000_parenthesized_values(line: &str, key: &str) -> Vec<Vec<String>> {
    let marker = format!("{key}=(");
    let mut result = Vec::new();
    let mut offset = 0;
    while let Some(relative) = line[offset..].find(&marker) {
        let start = offset + relative + marker.len();
        let Some(end) = line[start..].find(')') else {
            break;
        };
        result.push(
            line[start..start + end]
                .split_whitespace()
                .map(ToString::to_string)
                .collect(),
        );
        offset = start + end + 1;
    }
    result
}

fn v3000_value(line: &str, key: &str) -> Option<String> {
    let marker = format!("{key}=");
    let start = line.find(&marker)? + marker.len();
    let value = &line[start..];
    if let Some(value) = value.strip_prefix('"') {
        let mut result = String::new();
        let mut characters = value.chars().peekable();
        while let Some(character) = characters.next() {
            if character == '"' {
                if characters.peek() == Some(&'"') {
                    characters.next();
                    result.push('"');
                } else {
                    return Some(result);
                }
            } else {
                result.push(character);
            }
        }
        None
    } else {
        let end = value.find(char::is_whitespace).unwrap_or(value.len());
        Some(value[..end].to_string())
    }
}

fn v3000_values(line: &str, key: &str) -> Vec<String> {
    let marker = format!("{key}=");
    let mut values = Vec::new();
    let mut offset = 0;
    while let Some(relative) = line[offset..].find(&marker) {
        let start = offset + relative;
        if let Some(value) = v3000_value(&line[start..], key) {
            values.push(value);
        }
        offset = start + marker.len();
    }
    values
}

fn v3000_brackets(line: &str) -> Vec<[f64; 4]> {
    let marker = "BRKXYZ=(";
    let mut brackets = Vec::new();
    let mut offset = 0;
    while let Some(relative) = line[offset..].find(marker) {
        let start = offset + relative + marker.len();
        let Some(end) = line[start..].find(')') else {
            break;
        };
        let values = line[start..start + end]
            .split_whitespace()
            .skip(1)
            .filter_map(|value| value.parse::<f64>().ok())
            .collect::<Vec<_>>();
        if values.len() >= 6 {
            brackets.push([values[0], values[1], values[3], values[4]]);
        }
        offset = start + end + 1;
    }
    brackets
}

fn ordered_sdf_properties(raw_record: &str) -> Vec<SdfProperty> {
    let lines = raw_record.lines().collect::<Vec<_>>();
    let Some(data_start) = lines.iter().position(|line| line.trim_end() == "M  END") else {
        return Vec::new();
    };
    let mut properties = Vec::new();
    let mut index = data_start + 1;

    while index < lines.len() {
        let header = lines[index];
        if header.trim() == "$$$$" {
            break;
        }
        let Some(name) = property_name(header) else {
            index += 1;
            continue;
        };
        let line = index + 1;
        index += 1;
        let mut value_lines = Vec::new();
        while index < lines.len() {
            let value_line = lines[index];
            if value_line.trim() == "$$$$" || value_line.is_empty() {
                break;
            }
            if property_name(value_line).is_some() {
                break;
            }
            value_lines.push(value_line);
            index += 1;
        }
        properties.push(SdfProperty {
            name,
            value: value_lines.join("\n"),
            header: header.to_string(),
            line,
        });
        if index < lines.len() && lines[index].is_empty() {
            index += 1;
        }
    }

    properties
}

fn property_name(header: &str) -> Option<String> {
    let header = header.trim_start();
    if !header.starts_with('>') {
        return None;
    }
    let start = header.find('<')? + 1;
    let end = header[start..].find('>')? + start;
    let name = header[start..end].trim();
    (!name.is_empty()).then(|| name.to_string())
}

fn fidelity_diagnostics(
    format: SdfFormat,
    raw_record: &str,
    atoms: &[ChemicalAtom],
    bonds: &[ChemicalBond],
    sgroups: &[ChemicalSGroup],
    link_nodes: &[ChemicalLinkNode],
    collections: &[ChemicalCollection],
) -> Vec<FidelityDiagnostic> {
    let mut diagnostics = Vec::new();
    let mut seen = HashSet::new();

    let contractible_superatoms = contractible_superatom_indexes(sgroups)
        .into_iter()
        .collect::<HashSet<_>>();
    let uncontracted_superatoms = sgroups
        .iter()
        .enumerate()
        .filter(|group| {
            let (index, group) = group;
            group.kind == "superatom"
                && group.atom_source_ids.len() > 1
                && group.expanded != Some(true)
                && !contractible_superatoms.contains(index)
        })
        .map(|(_, group)| group)
        .collect::<Vec<_>>();
    if !uncontracted_superatoms.is_empty() {
        push_diagnostic(
            &mut diagnostics,
            &mut seen,
            "superatom-contraction-not-depicted",
            "superatom",
            "Multi-atom superatom SGroups without a usable label or with overlapping atom membership cannot be contracted",
            uncontracted_superatoms
                .iter()
                .flat_map(|group| group.atom_source_ids.iter().copied())
                .collect(),
            uncontracted_superatoms
                .iter()
                .flat_map(|group| group.crossing_bond_source_ids.iter().copied())
                .collect(),
        );
    }

    let misplaced_attachment_atoms = atoms
        .iter()
        .filter(|atom| !atom.attachment_points.is_empty())
        .map(|atom| atom.source_id)
        .collect::<Vec<_>>();
    if !misplaced_attachment_atoms.is_empty() {
        push_diagnostic(
            &mut diagnostics,
            &mut seen,
            "rgroup-attachment-context-not-supported",
            "attachment-point",
            "M APO and ATTCHPT values belong to R-group member CTABs; values found on the main CTAB are preserved but are not a standard main-structure depiction",
            misplaced_attachment_atoms,
            Vec::new(),
        );
    }

    let atom_ids = atoms
        .iter()
        .map(|atom| atom.source_id)
        .collect::<HashSet<_>>();
    let invalid_multicenter_bonds = bonds
        .iter()
        .filter(|bond| {
            let has_endpoints = !bond.endpoint_source_ids.is_empty();
            let valid_mode = matches!(bond.attachment_mode.as_deref(), Some("ALL" | "ANY"));
            has_endpoints != bond.attachment_mode.is_some()
                || (has_endpoints
                    && (bond.endpoint_source_ids.len() < 2
                        || !valid_mode
                        || bond
                            .endpoint_source_ids
                            .iter()
                            .any(|id| !atom_ids.contains(id))))
        })
        .collect::<Vec<_>>();
    if !invalid_multicenter_bonds.is_empty() {
        push_diagnostic(
            &mut diagnostics,
            &mut seen,
            "invalid-multicenter-bond-semantics",
            "multicenter-bond",
            "ENDPTS requires at least two valid atom IDs and an ATTACH mode of ALL or ANY",
            invalid_multicenter_bonds
                .iter()
                .flat_map(|bond| bond.endpoint_source_ids.iter().copied())
                .collect(),
            invalid_multicenter_bonds
                .iter()
                .map(|bond| bond.source_id)
                .collect(),
        );
    }

    let invalid_link_nodes = link_nodes
        .iter()
        .filter(|node| {
            node.min_repeat != 1
                || node.max_repeat == 0
                || !link_node_atoms_are_valid(node, &atom_ids)
        })
        .collect::<Vec<_>>();
    if !invalid_link_nodes.is_empty() {
        push_diagnostic(
            &mut diagnostics,
            &mut seen,
            "invalid-link-node-semantics",
            "link-node",
            "LINKNODE requires minRepeat=1, a positive maxRepeat, and valid atom connections",
            invalid_link_nodes
                .iter()
                .flat_map(|node| {
                    node.connections.iter().flat_map(|connection| {
                        [connection.atom_source_id, connection.neighbor_source_id]
                    })
                })
                .collect(),
            Vec::new(),
        );
    }

    for group in sgroups.iter().filter(|group| group.kind == "unknown") {
        push_diagnostic(
            &mut diagnostics,
            &mut seen,
            "unknown-sgroup-semantics",
            "sgroup",
            &format!(
                "SGroup type {} is preserved by inspection but is not depicted because its semantics are unknown",
                group.type_code
            ),
            group.atom_source_ids.clone(),
            group
                .bond_source_ids
                .iter()
                .chain(&group.crossing_bond_source_ids)
                .copied()
                .collect(),
        );
    }

    for collection in collections
        .iter()
        .filter(|collection| collection.kind != "highlight")
    {
        let kind = collection.kind.as_str();
        push_diagnostic(
            &mut diagnostics,
            &mut seen,
            &format!("{kind}-collection-not-depicted"),
            kind,
            &format!(
                "COLLECTION {} is preserved by inspection but has no standard default depiction",
                collection.name
            ),
            collection.atom_source_ids.clone(),
            collection.bond_source_ids.clone(),
        );
    }

    let valid_bond_ids = bonds
        .iter()
        .map(|bond| bond.source_id)
        .collect::<HashSet<_>>();
    let invalid_collection_members = collections
        .iter()
        .filter(|collection| {
            collection
                .atom_source_ids
                .iter()
                .any(|id| !atom_ids.contains(id))
                || collection
                    .bond_source_ids
                    .iter()
                    .any(|id| !valid_bond_ids.contains(id))
        })
        .collect::<Vec<_>>();
    if !invalid_collection_members.is_empty() {
        push_diagnostic(
            &mut diagnostics,
            &mut seen,
            "invalid-collection-members",
            "collection",
            "A COLLECTION references atom or bond IDs outside the current connection-table context",
            invalid_collection_members
                .iter()
                .flat_map(|collection| collection.atom_source_ids.iter().copied())
                .collect(),
            invalid_collection_members
                .iter()
                .flat_map(|collection| collection.bond_source_ids.iter().copied())
                .collect(),
        );
    }

    let partially_depicted_highlights = collections
        .iter()
        .filter(|collection| {
            collection.kind == "highlight"
                && (!collection.sgroup_source_ids.is_empty()
                    || !collection.object3d_source_ids.is_empty()
                    || !collection.rgroup_ids.is_empty()
                    || !collection.members.is_empty())
        })
        .collect::<Vec<_>>();
    if !partially_depicted_highlights.is_empty() {
        push_diagnostic(
            &mut diagnostics,
            &mut seen,
            "highlight-members-not-depicted",
            "highlight",
            "Non-atom/bond HILITE members are preserved by inspection but are not depicted",
            partially_depicted_highlights
                .iter()
                .flat_map(|collection| collection.atom_source_ids.iter().copied())
                .collect(),
            partially_depicted_highlights
                .iter()
                .flat_map(|collection| collection.bond_source_ids.iter().copied())
                .collect(),
        );
    }

    let reacting_center_bonds = bonds
        .iter()
        .filter_map(|bond| bond.reacting_center.map(|_| bond.source_id))
        .collect::<Vec<_>>();
    if !reacting_center_bonds.is_empty() || raw_record.contains(" RXCTR=") {
        push_diagnostic(
            &mut diagnostics,
            &mut seen,
            "reaction-center-not-depicted",
            "reaction-center",
            "Reaction-center metadata is preserved by inspection but is not depicted",
            Vec::new(),
            reacting_center_bonds,
        );
    }

    if format == SdfFormat::V3000 {
        let mut invalid_query_atoms = atoms
            .iter()
            .filter(|atom| {
                atom.query.as_ref().is_some_and(|query| {
                    query
                        .hydrogen_count
                        .is_some_and(|value| !(0..=4).contains(&value))
                        || query
                            .valence
                            .is_some_and(|value| !(-1..=14).contains(&value))
                        || query
                            .substitution_count
                            .is_some_and(|value| !(-2..=6).contains(&value))
                        || query
                            .ring_bond_count
                            .is_some_and(|value| !(-2..=4).contains(&value))
                })
            })
            .map(|atom| atom.source_id)
            .chain(raw_v3000_invalid_unsaturated_atoms(raw_record))
            .collect::<HashSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        invalid_query_atoms.sort_unstable();
        if !invalid_query_atoms.is_empty() {
            push_diagnostic(
                &mut diagnostics,
                &mut seen,
                "invalid-v3000-query-count",
                "query-constraint",
                "One or more V3000 query values are outside the CTfile value range",
                invalid_query_atoms,
                Vec::new(),
            );
        }
    }

    diagnostics
}

fn raw_v3000_invalid_unsaturated_atoms(raw_record: &str) -> Vec<u32> {
    let mut invalid = Vec::new();
    let mut in_atom_block = false;
    for line in v3000_logical_lines(raw_record) {
        if line == "BEGIN ATOM" {
            in_atom_block = true;
            continue;
        }
        if line == "END ATOM" {
            break;
        }
        if !in_atom_block
            || !v3000_integer(&line, "UNSAT").is_some_and(|value| !matches!(value, 0 | 1))
        {
            continue;
        }
        if let Some(source_id) = line
            .split_whitespace()
            .next()
            .and_then(|value| value.parse().ok())
        {
            invalid.push(source_id);
        }
    }
    invalid
}

pub(crate) fn contractible_superatom_indexes(sgroups: &[ChemicalSGroup]) -> Vec<usize> {
    let mut claimed_atoms = HashSet::new();
    let mut indexes = Vec::new();

    for (index, group) in sgroups.iter().enumerate() {
        let unique_atoms = group
            .atom_source_ids
            .iter()
            .copied()
            .collect::<HashSet<_>>();
        let has_label = group
            .label
            .as_ref()
            .or(group.subscript.as_ref())
            .is_some_and(|label| !label.trim().is_empty());
        if group.kind != "superatom"
            || unique_atoms.len() <= 1
            || group.expanded == Some(true)
            || !has_label
            || unique_atoms
                .iter()
                .any(|source_id| claimed_atoms.contains(source_id))
        {
            continue;
        }

        claimed_atoms.extend(unique_atoms);
        indexes.push(index);
    }

    indexes
}

#[allow(clippy::too_many_arguments)]
fn push_diagnostic(
    diagnostics: &mut Vec<FidelityDiagnostic>,
    seen: &mut HashSet<String>,
    code: &str,
    feature: &str,
    message: &str,
    atom_source_ids: Vec<u32>,
    bond_source_ids: Vec<u32>,
) {
    if !seen.insert(code.to_string()) {
        return;
    }
    diagnostics.push(FidelityDiagnostic {
        code: code.to_string(),
        severity: DiagnosticSeverity::Warning,
        feature: feature.to_string(),
        message: message.to_string(),
        atom_source_ids,
        bond_source_ids,
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn properties_preserve_order_duplicates_multiline_values_and_headers() {
        let sdf = concat!(
            "record\n  molchemist\ncomment\n",
            "  0  0  0  0  0  0  0  0  0  0  0 V2000\n",
            "M  END\n",
            ">  <NAME> (registry)\nfirst\nline two\n\n",
            "> <NAME>\nsecond\n\n",
            "> <EMPTY>\n\n",
            "$$$$\n",
        );

        assert_eq!(
            ordered_sdf_properties(sdf),
            vec![
                SdfProperty {
                    name: "NAME".to_string(),
                    value: "first\nline two".to_string(),
                    header: ">  <NAME> (registry)".to_string(),
                    line: 6,
                },
                SdfProperty {
                    name: "NAME".to_string(),
                    value: "second".to_string(),
                    header: "> <NAME>".to_string(),
                    line: 10,
                },
                SdfProperty {
                    name: "EMPTY".to_string(),
                    value: String::new(),
                    header: "> <EMPTY>".to_string(),
                    line: 13,
                },
            ]
        );
    }

    #[test]
    fn strict_policy_rejects_diagnostics_and_other_policies_do_not() {
        let diagnostic = FidelityDiagnostic {
            code: "sgroup-not-depicted".to_string(),
            severity: DiagnosticSeverity::Warning,
            feature: "sgroup".to_string(),
            message: "SGroup is not depicted".to_string(),
            atom_source_ids: vec![1],
            bond_source_ids: Vec::new(),
        };
        assert!(FidelityPolicy::Ignore
            .enforce(std::slice::from_ref(&diagnostic))
            .is_ok());
        assert!(FidelityPolicy::Warn
            .enforce(std::slice::from_ref(&diagnostic))
            .is_ok());
        assert_eq!(
            FidelityPolicy::Strict.enforce(&[diagnostic]).unwrap_err(),
            "faithful depiction is not available: SGroup is not depicted"
        );
    }

    #[test]
    fn v3000_data_sgroup_preserves_ordered_field_values_and_display_metadata() {
        let sdf = concat!(
            "M  V30 BEGIN SGROUP\n",
            "M  V30 8 DAT 0 ATOMS=(1 10) FIELDNAME=Temperature FIELDDATA=\"25 C\" FIELDDATA=ambient BRKTYP=BRACKET CLASS=measurement COMPNO=2\n",
            "M  V30 END SGROUP\n",
        );
        let groups = raw_v3000_sgroups(sdf, &[10], &[]);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].kind, "data");
        assert_eq!(groups[0].field_name.as_deref(), Some("Temperature"));
        assert_eq!(groups[0].field_value.as_deref(), Some("25 C"));
        assert_eq!(groups[0].field_values, vec!["25 C", "ambient"]);
        assert_eq!(groups[0].bracket_type.as_deref(), Some("BRACKET"));
        assert_eq!(groups[0].class.as_deref(), Some("measurement"));
        assert_eq!(groups[0].component_number, Some(2));
    }

    #[test]
    fn v3000_sap_entries_preserve_order_ids_and_aidx_leaving_atom() {
        let sdf = concat!(
            "M  V30 BEGIN SGROUP\n",
            "M  V30 3 SUP 0 ATOMS=(2 10 20) LABEL=A SAP=(3 10 AIDX Al) SAP=(3 20 30 Br)\n",
            "M  V30 END SGROUP\n",
        );
        let groups = raw_v3000_sgroups(sdf, &[10, 20, 30], &[]);
        assert_eq!(
            groups[0].attachment_points,
            vec![
                ChemicalSGroupAttachmentPoint {
                    atom_source_id: 10,
                    leaving_atom_source_id: Some(10),
                    id: "Al".to_string(),
                },
                ChemicalSGroupAttachmentPoint {
                    atom_source_id: 20,
                    leaving_atom_source_id: Some(30),
                    id: "Br".to_string(),
                },
            ]
        );
    }

    #[test]
    fn v3000_hcount_uses_ctfile_query_encoding() {
        assert_eq!(v3000_hydrogen_count("1 C 0 0 0 0 HCOUNT=-1"), Some(0));
        assert_eq!(v3000_hydrogen_count("1 C 0 0 0 0 HCOUNT=1"), Some(1));
        assert_eq!(v3000_hydrogen_count("1 C 0 0 0 0 HCOUNT=2"), Some(2));
        assert_eq!(v3000_hydrogen_count("1 C 0 0 0 0 HCOUNT=3"), Some(3));
        assert_eq!(v3000_hydrogen_count("1 C 0 0 0 0 HCOUNT=4"), Some(4));
        assert_eq!(v3000_hydrogen_count("1 C 0 0 0 0 HCOUNT=5"), Some(5));
        assert_eq!(v3000_hydrogen_count("1 C 0 0 0 0 HCOUNT=0"), None);
    }

    #[test]
    fn query_count_sentinels_and_upper_buckets_are_normalized() {
        assert_eq!(normalize_substitution_count(-2), Some(-2));
        assert_eq!(normalize_substitution_count(-1), Some(-1));
        assert_eq!(normalize_substitution_count(0), None);
        assert_eq!(normalize_substitution_count(5), Some(5));
        assert_eq!(normalize_substitution_count(6), Some(6));
        assert_eq!(normalize_substitution_count(12), Some(6));

        assert_eq!(normalize_ring_bond_count(-2), Some(-2));
        assert_eq!(normalize_ring_bond_count(-1), Some(-1));
        assert_eq!(normalize_ring_bond_count(0), None);
        assert_eq!(normalize_ring_bond_count(3), Some(3));
        assert_eq!(normalize_ring_bond_count(4), Some(4));
        assert_eq!(normalize_ring_bond_count(12), Some(4));
    }

    #[test]
    fn v3000_collections_preserve_standard_and_user_defined_members() {
        let sdf = concat!(
            "M  V30 BEGIN COLLECTION\n",
            "M  V30 MDLV30/HILITE ATOMS=(1 10) BONDS=(1 7) SGROUPS=(1 3) OBJ3DS=(1 5) RGROUPS=(1 R2) MEMBERS=(1 ROOT)\n",
            "M  V30 \"AUDIT/SET 1\" ATOMS=(2 10 20) MEMBERS=(2 ROOT R2M3)\n",
            "M  V30 MDLV30/FUTURE BONDS=(1 7)\n",
            "M  V30 mdlv30/hilite ATOMS=(1 20)\n",
            "M  V30 END COLLECTION\n",
        );
        let collections = raw_v3000_collections(sdf);
        assert_eq!(collections.len(), 4);
        assert_eq!(collections[0].name, "MDLV30/HILITE");
        assert_eq!(collections[0].kind, "highlight");
        assert_eq!(collections[0].atom_source_ids, vec![10]);
        assert_eq!(collections[0].bond_source_ids, vec![7]);
        assert_eq!(collections[0].sgroup_source_ids, vec![3]);
        assert_eq!(collections[0].object3d_source_ids, vec![5]);
        assert_eq!(collections[0].rgroup_ids, vec!["R2"]);
        assert_eq!(collections[0].members, vec!["ROOT"]);
        assert_eq!(collections[1].name, "AUDIT/SET 1");
        assert_eq!(collections[1].kind, "user-defined");
        assert_eq!(collections[1].members, vec!["ROOT", "R2M3"]);
        assert_eq!(collections[2].kind, "internal-unknown");
        assert_eq!(collections[2].raw_entry, "MDLV30/FUTURE BONDS=(1 7)");
        assert_eq!(collections[3].name, "mdlv30/hilite");
        assert_eq!(collections[3].kind, "highlight");
    }

    #[test]
    fn standard_and_unknown_sgroup_type_codes_are_not_dropped() {
        let sdf = concat!(
            "M  V30 BEGIN SGROUP\n",
            "M  V30 1 CRO 0 ATOMS=(1 10)\n",
            "M  V30 2 MOD 0 ATOMS=(1 10)\n",
            "M  V30 3 GRA 0 ATOMS=(1 10)\n",
            "M  V30 4 ANY 0 ATOMS=(1 10)\n",
            "M  V30 5 ZZZ 0 ATOMS=(1 10)\n",
            "M  V30 END SGROUP\n",
        );
        let groups = raw_v3000_sgroups(sdf, &[10], &[]);
        assert_eq!(groups.len(), 5);
        assert_eq!(groups[0].kind, "crosslink");
        assert_eq!(groups[1].kind, "modification");
        assert_eq!(groups[2].kind, "graft");
        assert_eq!(groups[3].kind, "any");
        assert_eq!(groups[4].kind, "unknown");
        assert_eq!(groups[4].type_code, "ZZZ");
    }
}
