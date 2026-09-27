//! Input adapters. These preserve source text and reject malformed records
//! before passing normalized CTABs to the CTfile semantic parser.
use std::collections::{BTreeMap, HashSet};
use std::fmt::Write;

use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct MoleculeInput {
    pub format: String,
    pub data: String,
}

pub fn mol2_to_sdf(input: &str, record: usize) -> Result<String, String> {
    if record == 0 {
        return Err("MOL2 record numbers are one-based".into());
    }
    let records = input.split("@<TRIPOS>MOLECULE").skip(1).collect::<Vec<_>>();
    let data = records
        .get(record - 1)
        .ok_or("MOL2 record does not exist")?;
    let mut sections = BTreeMap::<String, Vec<&str>>::new();
    let mut section = "MOLECULE".to_string();
    for line in data.trim_start_matches(['\r', '\n']).lines() {
        if let Some(name) = line.trim().strip_prefix("@<TRIPOS>") {
            section = name.to_string();
        } else {
            sections.entry(section.clone()).or_default().push(line);
        }
    }
    let header = sections.get("MOLECULE").ok_or("Missing MOL2 header")?;
    let counts = header
        .get(1)
        .ok_or("Missing MOL2 counts")?
        .split_whitespace()
        .collect::<Vec<_>>();
    let count = |i: usize| -> Result<usize, String> {
        counts
            .get(i)
            .ok_or("Missing MOL2 count")?
            .parse()
            .map_err(|_| "Invalid MOL2 count".into())
    };
    let rows = |name: &str| {
        sections
            .get(name)
            .into_iter()
            .flatten()
            .copied()
            .filter(|s| !s.trim().is_empty() && !s.trim_start().starts_with('#'))
            .collect::<Vec<_>>()
    };
    let atom_rows = rows("ATOM");
    let bond_rows = rows("BOND");
    if count(0)? != atom_rows.len() || count(1)? != bond_rows.len() {
        return Err("MOL2 atom/bond counts do not match the records".into());
    }
    if atom_rows.is_empty() {
        return Err("MOL2 structure is empty".into());
    }
    let mut atoms = Vec::new();
    let mut ids = HashSet::new();
    for row in atom_rows {
        let p = row.split_whitespace().collect::<Vec<_>>();
        if p.len() < 6 {
            return Err("Incomplete MOL2 atom".into());
        }
        let id = p[0].parse::<u32>().map_err(|_| "Invalid MOL2 atom ID")?;
        if id == 0 || !ids.insert(id) {
            return Err("Duplicate or zero MOL2 atom ID".into());
        }
        for value in &p[2..5] {
            if !value.parse::<f64>().is_ok_and(f64::is_finite) {
                return Err("Non-finite or invalid MOL2 coordinate".into());
            }
        }
        if p.get(8)
            .is_some_and(|v| !v.parse::<f64>().is_ok_and(f64::is_finite))
        {
            return Err("Invalid MOL2 partial charge".into());
        }
        let element = p[5].split('.').next().unwrap();
        let element = match element {
            "Du" | "Any" | "LP" | "Hal" | "Het" | "Hev" => "*",
            other => other,
        };
        if element.parse::<crate::AtomSymbol>().is_err() {
            return Err(format!("Unknown MOL2 atom type {}", p[5]));
        }
        // MOL2 charge columns are partial charges, not formal charges.
        let charge = if p[5] == "N.4" { " CHG=1" } else { "" };
        atoms.push(format!(
            "M  V30 {id} {element} {} {} {} 0{charge}\n",
            p[2], p[3], p[4]
        ));
    }
    let mut bond_ids = HashSet::new();
    let mut edges = HashSet::new();
    let mut bonds = Vec::new();
    for row in bond_rows {
        let p = row.split_whitespace().collect::<Vec<_>>();
        if p.len() < 4 {
            return Err("Incomplete MOL2 bond".into());
        }
        let parse = |i: usize| {
            p[i].parse::<u32>()
                .map_err(|_| "Invalid MOL2 bond ID or endpoint")
        };
        let (id, a, b) = (parse(0)?, parse(1)?, parse(2)?);
        if id == 0
            || !bond_ids.insert(id)
            || a == b
            || !ids.contains(&a)
            || !ids.contains(&b)
            || !edges.insert((a.min(b), a.max(b)))
        {
            return Err("Duplicate, self, or out-of-range MOL2 bond".into());
        }
        let order = match p[3] {
            "1" | "am" => 1,
            "2" => 2,
            "3" => 3,
            "ar" => 4,
            "du" | "un" => 8,
            "nc" => continue,
            other => return Err(format!("Unknown MOL2 bond type {other}")),
        };
        bonds.push(format!("M  V30 {id} {order} {a} {b}\n"));
    }
    let mut out = format!("{}\n  molchemist\nMOL2 conversion; partial charges retained in MOL2_SOURCE\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS {} {} 0 0 0\nM  V30 BEGIN ATOM\n", header[0], atoms.len(), bonds.len());
    for atom in atoms {
        out.push_str(&atom);
    }
    out.push_str("M  V30 END ATOM\nM  V30 BEGIN BOND\n");
    for bond in bonds {
        out.push_str(&bond);
    }
    out.push_str("M  V30 END BOND\nM  V30 END CTAB\nM  END\n> <MOL2_SOURCE>\n");
    // SDF properties end at blank lines. A single JSON string also safely
    // preserves embedded blank lines, record separators and CRLF in MOL2.
    out.push_str(
        &serde_json::to_string(&format!("@<TRIPOS>MOLECULE{data}")).map_err(|e| e.to_string())?,
    );
    out.push_str("\n\n$$$$\n");
    crate::inspect_sdf_record(&out, 1)?;
    Ok(out)
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ReactionInput {
    pub raw_source: String,
    pub reactants: Vec<MoleculeInput>,
    pub agents: Vec<MoleculeInput>,
    pub products: Vec<MoleculeInput>,
}

pub fn parse_reaction(input: &str, format: &str) -> Result<ReactionInput, String> {
    if format == "reaction-smiles" {
        let parts = input.trim().split('>').collect::<Vec<_>>();
        if parts.len() != 3 || parts[0].is_empty() || parts[2].is_empty() {
            return Err("Reaction SMILES requires reactants>agents>products".into());
        }
        let side = |text: &str| -> Result<Vec<MoleculeInput>, String> {
            if text.is_empty() {
                return Ok(Vec::new());
            }
            text.split('.')
                .map(|data| {
                    crate::parse(data).map_err(|e| e.to_string())?;
                    Ok(MoleculeInput {
                        format: "smiles".into(),
                        data: data.into(),
                    })
                })
                .collect()
        };
        return Ok(ReactionInput {
            raw_source: input.to_string(),
            reactants: side(parts[0])?,
            agents: side(parts[1])?,
            products: side(parts[2])?,
        });
    }
    if format != "rxn" || !input.trim_start().starts_with("$RXN") {
        return Err("Expected RXN or reaction-smiles input".into());
    }
    let raw_source = input.to_string();
    let input = input.replace("\r\n", "\n");
    if input.lines().next().is_some_and(|l| l.contains("V3000")) {
        let mut result = ReactionInput {
            raw_source,
            reactants: Vec::new(),
            agents: Vec::new(),
            products: Vec::new(),
        };
        let mut side = "";
        let mut ctab = String::new();
        let mut depth = 0;
        let mut counts = None;
        for line in input.lines() {
            let text = line.strip_prefix("M  V30 ").unwrap_or(line);
            if depth == 0 && text.starts_with("COUNTS ") {
                let c = text
                    .split_whitespace()
                    .skip(1)
                    .map(str::parse::<usize>)
                    .collect::<Result<Vec<_>, _>>()
                    .map_err(|_| "Invalid RXN counts")?;
                if c.len() < 2 {
                    return Err("Missing RXN counts".into());
                }
                counts = Some((c[0], c[1], *c.get(2).unwrap_or(&0)));
            }
            match text {
                "BEGIN REACTANT" => side = "reactants",
                "BEGIN PRODUCT" => side = "products",
                "BEGIN AGENT" => side = "agents",
                "BEGIN CTAB" => {
                    depth += 1;
                    if depth == 1 {
                        ctab.clear();
                    }
                }
                _ => {}
            }
            if depth > 0 {
                writeln!(ctab, "{line}").unwrap();
            }
            if text == "END CTAB" {
                if depth == 0 {
                    return Err("Unmatched RXN CTAB end".into());
                }
                depth -= 1;
                if depth == 0 {
                    let data = format!("RXN component\n  molchemist\n\n  0  0  0  0  0  0            999 V3000\n{ctab}M  END\n");
                    crate::inspect_sdf_record(&data, 1)?;
                    let molecule = MoleculeInput {
                        format: "mol".into(),
                        data,
                    };
                    match side {
                        "reactants" => result.reactants.push(molecule),
                        "products" => result.products.push(molecule),
                        "agents" => result.agents.push(molecule),
                        _ => return Err("RXN CTAB outside a reaction side".into()),
                    }
                }
            }
        }
        if depth != 0
            || counts
                != Some((
                    result.reactants.len(),
                    result.products.len(),
                    result.agents.len(),
                ))
        {
            return Err("RXN counts or CTAB boundaries do not match".into());
        }
        return Ok(result);
    }
    let mut chunks = input.split("$MOL");
    let header = chunks.next().unwrap();
    let line = header.lines().nth(4).ok_or("Missing RXN counts")?;
    let c = line
        .split_whitespace()
        .map(str::parse::<usize>)
        .collect::<Result<Vec<_>, _>>()
        .map_err(|_| "Invalid RXN counts")?;
    if c.len() < 2 {
        return Err("Missing RXN counts".into());
    }
    let mut molecules = chunks
        .map(|s| {
            let data = s.strip_prefix('\n').unwrap_or(s).to_string();
            crate::inspect_sdf_record(&data, 1)?;
            Ok(MoleculeInput {
                format: "mol".into(),
                data,
            })
        })
        .collect::<Result<Vec<_>, String>>()?;
    let agents = *c.get(2).unwrap_or(&0);
    if molecules.len() != c[0] + c[1] + agents {
        return Err("RXN molecule count mismatch".into());
    }
    let rest = molecules.split_off(c[0]);
    let mut products = rest;
    let agents = products.split_off(c[1]);
    Ok(ReactionInput {
        raw_source,
        reactants: molecules,
        products,
        agents,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    const MOL2: &str = "@<TRIPOS>MOLECULE\nwater\n3 2 0 0 0\nSMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n10 O 0 0 0 O.3 1 HOH -0.8\n20 H 1 0 0 H 1 HOH 0.4\n30 H 0 1 0 H 1 HOH 0.4\n@<TRIPOS>BOND\n1 10 20 1\n2 10 30 1\n";
    #[test]
    fn mol2_preserves_ids_partial_charges_and_validates() {
        let sdf = mol2_to_sdf(MOL2, 1).unwrap();
        let record = crate::inspect_sdf_record(&sdf, 1).unwrap();
        assert_eq!(record.atoms[0].source_id, 10);
        assert_eq!(record.atoms[0].formal_charge, 0);
        assert!(sdf.contains("-0.8"));
        let source = MOL2.replace("USER_CHARGES\n", "USER_CHARGES\r\n\r\n");
        let normalized = mol2_to_sdf(&source, 1).unwrap();
        let inspected = crate::inspect_sdf_record(&normalized, 1).unwrap();
        let property = inspected
            .properties
            .iter()
            .find(|p| p.name == "MOL2_SOURCE")
            .unwrap();
        assert_eq!(
            serde_json::from_str::<String>(&property.value).unwrap(),
            source
        );
        assert!(mol2_to_sdf(&MOL2.replace("2 10 30 1", "2 10 40 1"), 1).is_err());
        assert!(mol2_to_sdf(&MOL2.replace("O 0 0 0", "O NaN 0 0"), 1).is_err());
    }
    #[test]
    fn reaction_smiles_checks_each_component() {
        let r = parse_reaction("[CH3:1][OH:2]>O>[CH2:1]=[O:2]", "reaction-smiles").unwrap();
        assert_eq!(r.agents.len(), 1);
        assert!(parse_reaction("C..O>>CO", "reaction-smiles").is_err());
    }
}
