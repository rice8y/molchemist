//! RGfile panels retain the root, alternatives, attachment sites and logic.
//! Each member remains associated with its source group and attachment sites.
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct RGroupMember {
    pub group: u32,
    pub data: String,
    pub attachments: Vec<(usize, Vec<i8>)>,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct RGroupCondition {
    pub group: u32,
    pub then_group: u32,
    pub rest_h: bool,
    pub occurrence: String,
}

fn parse_condition(group: u32, fields: &str) -> Result<RGroupCondition, String> {
    let mut values = fields.split_whitespace();
    let then_group = values
        .next()
        .ok_or("Missing R-group dependency")?
        .parse()
        .map_err(|_| "Invalid R-group dependency")?;
    let rest_h = match values.next() {
        Some("0") => false,
        Some("1") => true,
        _ => return Err("Invalid R-group RestH condition".into()),
    };
    let occurrence = values
        .collect::<Vec<_>>()
        .join("")
        .trim_matches('"')
        .to_string();
    let occurrence = if occurrence.is_empty() {
        ">0".into()
    } else {
        occurrence
    };
    Ok(RGroupCondition {
        group,
        then_group,
        rest_h,
        occurrence,
    })
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct RGroupDocument {
    pub raw_source: String,
    pub root: String,
    pub members: Vec<RGroupMember>,
    pub logic: Vec<String>,
    pub conditions: Vec<RGroupCondition>,
}

pub fn parse_rgroups(input: &str) -> Result<RGroupDocument, String> {
    let text = input.replace("\r\n", "\n");
    let mut root = None;
    let mut members = Vec::new();
    let mut logic = Vec::new();
    let mut conditions = Vec::new();
    let mut group = 0u32;
    let mut in_ctab = false;
    let mut ctab = Vec::new();
    let mut pending_group = false;
    let v3000 = !text.trim_start().starts_with("$MDL");
    for line in text.lines() {
        let value = line.strip_prefix("M  V30 ").unwrap_or(line).trim();
        if value.starts_with("BEGIN RGROUP ") {
            group = value
                .trim_start_matches("BEGIN RGROUP ")
                .parse()
                .map_err(|_| "Invalid RGROUP number")?;
            if group == 0 {
                return Err("RGROUP number must be positive".into());
            }
        } else if value == "$RGP" {
            pending_group = true;
        } else if pending_group {
            group = value.parse().map_err(|_| "Invalid RGP number")?;
            pending_group = false;
        } else if value == "END RGROUP" || value == "$END RGP" {
            group = 0;
        }
        if value.starts_with("RLOGIC ") || line.starts_with("M  LOG") {
            logic.push(format!("R{group}: {value}"));
            if let Some(fields) = value.strip_prefix("RLOGIC ") {
                conditions.push(parse_condition(group, fields)?);
            } else {
                let fields = line
                    .get(6..)
                    .unwrap_or("")
                    .split_whitespace()
                    .collect::<Vec<_>>();
                if fields.len() < 4 || fields[0] != "1" {
                    return Err("Invalid V2000 R-group logic line".into());
                }
                let id = fields[1]
                    .parse()
                    .map_err(|_| "Invalid R-group logic number")?;
                conditions.push(parse_condition(id, &fields[2..].join(" "))?);
            }
        }
        if value == "BEGIN CTAB" || value == "$CTAB" {
            if in_ctab {
                return Err("Unexpected nested CTAB".into());
            }
            in_ctab = true;
            ctab.clear();
            if v3000 {
                ctab.push(line.to_string());
            }
        } else if value == "END CTAB" || value == "$END CTAB" {
            if !in_ctab {
                return Err("Unmatched CTAB end".into());
            }
            if v3000 {
                ctab.push(line.to_string());
            }
            let body = ctab.join("\n");
            let mut data = if v3000 {
                format!("R-group member\n  molchemist\n\n  0  0  0  0  0  0            999 V3000\n{body}\nM  END\n")
            } else {
                format!("R-group member\n  molchemist\n\n{body}\n")
            };
            if !data.contains("M  END") {
                data.push_str("M  END\n");
            }
            let inspection = crate::inspect_sdf_record(&data, 1)?;
            if root.is_none() && group == 0 {
                root = Some(data);
            } else {
                if group == 0 {
                    return Err("R-group member outside a numbered group".into());
                }
                let attachments = inspection
                    .atoms
                    .iter()
                    .filter(|a| !a.attachment_points.is_empty())
                    .map(|a| (a.index, a.attachment_points.clone()))
                    .collect();
                members.push(RGroupMember {
                    group,
                    data,
                    attachments,
                });
            }
            in_ctab = false;
        } else if in_ctab {
            ctab.push(line.to_string());
        }
    }
    if in_ctab || pending_group {
        return Err("Truncated RGfile".into());
    }
    let root = root.ok_or("RGfile has no root CTAB")?;
    for member in &members {
        if !conditions.iter().any(|c| c.group == member.group) {
            conditions.push(parse_condition(member.group, "0 0 >0")?);
        }
    }
    Ok(RGroupDocument {
        raw_source: input.to_string(),
        root,
        members,
        logic,
        conditions,
    })
}

/// Only the depiction copy changes; inspection retains the exact source.
pub fn expand_superatoms(input: &str, record: usize) -> Result<String, String> {
    let inspected = crate::inspect_sdf_record(input, record)?;
    let mut out = String::new();
    let ids = inspected
        .sgroups
        .iter()
        .filter(|g| g.kind == "superatom")
        .map(|g| g.id)
        .collect::<Vec<_>>();
    let mut in_groups = false;
    let mut lines = inspected.raw_record.lines();
    while let Some(physical) = lines.next() {
        let mut line = physical.to_string();
        while line.starts_with("M  V30 ") && line.ends_with('-') {
            line.pop();
            line.push_str(
                lines
                    .next()
                    .and_then(|next| next.strip_prefix("M  V30 "))
                    .ok_or("Truncated V3000 continuation")?,
            );
        }
        if line == "M  V30 BEGIN SGROUP" {
            in_groups = true;
        }
        if line == "M  V30 END SGROUP" {
            in_groups = false;
        }
        if in_groups && line.split_whitespace().nth(3) == Some("SUP") {
            // Locate a top-level attribute, leaving quoted labels unchanged.
            let mut quoted = false;
            let mut state = None;
            for (index, ch) in line.char_indices() {
                if ch == '"' {
                    quoted = !quoted;
                }
                if !quoted && line[index..].starts_with(" ESTATE=") {
                    let start = index + 1;
                    let end = line[start..]
                        .find(char::is_whitespace)
                        .map_or(line.len(), |n| start + n);
                    state = Some(start..end);
                    break;
                }
            }
            if let Some(range) = state {
                line.replace_range(range, "ESTATE=E");
            } else {
                line.push_str(" ESTATE=E");
            }
            out.push_str(&line);
            out.push('\n');
        } else {
            if line == "M  END" && inspected.format == crate::ChemicalFormat::MolfileV2000 {
                for chunk in ids.chunks(8) {
                    out.push_str(&format!("M  SDS EXP{:>3}", chunk.len()));
                    for id in chunk {
                        out.push_str(&format!("{id:>4}"));
                    }
                    out.push('\n');
                }
            }
            out.push_str(&line);
            out.push('\n');
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn member_attachment_points_have_context() {
        let ctab = |atom: &str| {
            format!("M  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 {atom}\nM  V30 END ATOM\nM  V30 END CTAB\n")
        };
        let input=format!("rgroup\n\n\n  0  0  0  0  0  0            999 V3000\n{}M  V30 BEGIN RGROUP 1\nM  V30 RLOGIC 0 1 >0\n{}M  V30 END RGROUP\nM  END\n",ctab("1 R# 0 0 0 0 RGROUPS=(1 1)"),ctab("1 C 0 0 0 0 ATTCHPT=1"));
        let doc = parse_rgroups(&input).unwrap();
        assert_eq!(doc.members[0].attachments, vec![(0, vec![1])]);
        assert_eq!(doc.logic.len(), 1);
        assert!(doc.conditions[0].rest_h);
        assert_eq!(doc.conditions[0].occurrence, ">0");
        let condition = parse_condition(2, "3 0 \"> 0\"").unwrap();
        assert_eq!(condition.then_group, 3);
        assert_eq!(condition.occurrence, ">0");
    }

    #[test]
    fn expansion_preserves_continued_superatom_labels() {
        let mol = "SUP\n\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 1 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 2 O 1 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 END BOND\nM  V30 BEGIN SGROUP\nM  V30 1 SUP 0 ATOMS=(2 1 2) -\nM  V30 LABEL=\"test  label\" ESTATE=C\nM  V30 END SGROUP\nM  V30 END CTAB\nM  END\n";
        let expanded = expand_superatoms(mol, 1).unwrap();
        assert!(expanded.contains("LABEL=\"test  label\" ESTATE=E"));
        assert!(crate::inspect_sdf_record(&expanded, 1)
            .unwrap()
            .diagnostics
            .is_empty());
        assert!(expand_superatoms(&expanded, 1)
            .unwrap()
            .contains("ESTATE=E"));
    }
}
