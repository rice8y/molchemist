use std::fmt::Write;

use crate::{AtomLabel, Command, LinkData};

pub const DEFAULT_ALCHEMIST_IMPORT: &str = "@preview/alchemist:0.2.0";

const EXTENDED_BOND_DEFINITIONS: &str = r#"#import "@preview/cetz:0.5.2"

#let _molchemist-dashed-stroke(stroke, dash) = {
  if stroke == none or stroke == auto {
    stroke
  } else if type(stroke) == dictionary {
    stroke + (dash: dash)
  } else if type(stroke) == color {
    (paint: stroke, dash: dash)
  } else if type(stroke) == length {
    (thickness: stroke, dash: dash)
  } else {
    (paint: stroke.paint, thickness: stroke.thickness, dash: dash)
  }
}

#let _molchemist-partial-double(dash) = build-link((length, ctx, cetz-ctx, args) => {
  let args = args
  let offset = args.at("offset", default: ctx.config.double.offset)
  let key = if offset == "left" { "stroke-left" } else { "stroke-right" }
  let stroke = args.at(key, default: args.at("stroke", default: ctx.config.double.stroke))
  args.insert(key, _molchemist-dashed-stroke(stroke, dash))
  (double(..args).first().draw)(length, ctx, cetz-ctx)
})

#let _molchemist-dashed-single(dash) = build-link((length, ctx, cetz-ctx, args) => {
  let args = args
  let stroke = args.at("stroke", default: ctx.config.single.stroke)
  args.insert("stroke", _molchemist-dashed-stroke(stroke, dash))
  (single(..args).first().draw)(length, ctx, cetz-ctx)
})

#let _molchemist-dashed-double = build-link((length, ctx, cetz-ctx, args) => {
  let args = args
  let stroke = args.at("stroke", default: ctx.config.double.stroke)
  args.insert("stroke-left", _molchemist-dashed-stroke(args.at("stroke-left", default: stroke), "dashed"))
  args.insert("stroke-right", _molchemist-dashed-stroke(args.at("stroke-right", default: stroke), "dashed"))
  (double(..args).first().draw)(length, ctx, cetz-ctx)
})

#let _molchemist-wavy = build-link((length, ctx, _, args) => {
  import cetz.draw: *
  cetz.decorations.wave(
    line((0, 0), (length, 0), stroke: none),
    segments: 8,
    amplitude: .12,
    stroke: args.at("stroke", default: ctx.config.single.stroke),
  )
})

#let _molchemist-coordination-right = build-link((length, ctx, _, args) => {
  import cetz.draw: *
  line(
    (0, 0),
    (length, 0),
    stroke: args.at("stroke", default: ctx.config.single.stroke),
    mark: (end: ">"),
  )
})

#let _molchemist-coordination-left = build-link((length, ctx, _, args) => {
  import cetz.draw: *
  line(
    (0, 0),
    (length, 0),
    stroke: args.at("stroke", default: ctx.config.single.stroke),
    mark: (start: ">"),
  )
})

#let _molchemist-aromatic = _molchemist-partial-double("dashed")
#let _molchemist-single-or-double = _molchemist-partial-double("dotted")
#let _molchemist-single-or-aromatic = _molchemist-dashed-single("dashed")
#let _molchemist-double-or-aromatic = _molchemist-dashed-double
#let _molchemist-hydrogen = _molchemist-dashed-single("dotted")

"#;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct StandaloneOptions {
    pub alchemist_import: String,
    pub page_margin: String,
}

impl Default for StandaloneOptions {
    fn default() -> Self {
        Self {
            alchemist_import: DEFAULT_ALCHEMIST_IMPORT.to_string(),
            page_margin: "3mm".to_string(),
        }
    }
}

pub fn format_alchemist(commands: &[Command], base_sep: &str, indent_width: usize) -> String {
    let mut output = format!("#let base-sep = {base_sep}\n");
    if has_extended_bonds(commands) {
        output.push_str(EXTENDED_BOND_DEFINITIONS);
    }
    output.push_str("#skeletize({\n");
    format_commands(&mut output, commands, 1, indent_width);
    output.push_str("})");
    output
}

pub fn format_standalone(
    commands: &[Command],
    base_sep: &str,
    indent_width: usize,
    options: &StandaloneOptions,
) -> String {
    format_standalone_code(&format_alchemist(commands, base_sep, indent_width), options)
}

pub fn format_standalone_code(code: &str, options: &StandaloneOptions) -> String {
    format!(
        "#import \"{}\": *\n\n#set page(width: auto, height: auto, margin: {})\n\n{}",
        escape_string(&options.alchemist_import),
        options.page_margin,
        code,
    )
}

fn format_commands(output: &mut String, commands: &[Command], depth: usize, indent_width: usize) {
    let indent = " ".repeat(depth * indent_width);
    for command in commands {
        match command {
            Command::Fragment {
                element,
                name,
                links,
                atom,
                ..
            } => {
                let links_text = format_links(links, depth, indent_width);
                if !element.is_empty() {
                    let mut arguments = Vec::new();
                    if !name.is_empty() {
                        arguments.push(format!("name: \"{}\"", escape_string(name)));
                    }
                    if !links_text.is_empty() {
                        arguments.push(links_text);
                    }

                    let body = atom.as_ref().map_or_else(
                        || format!("\"{}\"", escape_string(element)),
                        format_atom_label,
                    );
                    write!(output, "{indent}fragment({body}").unwrap();
                    if !arguments.is_empty() {
                        write!(output, ", {}", arguments.join(", ")).unwrap();
                    }
                    output.push_str(")\n");
                } else {
                    if !name.is_empty() {
                        writeln!(output, "{indent}hook(\"{}\")", escape_string(name)).unwrap();
                    }
                    if !links_text.is_empty() {
                        writeln!(output, "{indent}branch({{").unwrap();
                        let inner = " ".repeat((depth + 1) * indent_width);
                        writeln!(
                            output,
                            "{inner}single(absolute: 0deg, atom-sep: 0pt, stroke: none, name: \"{}-links\", {links_text})",
                            escape_string(name),
                        )
                        .unwrap();
                        writeln!(output, "{indent}}})").unwrap();
                    }
                }
            }
            Command::Bond {
                name,
                bond_type,
                angle,
                offset,
                length_scale,
            } => {
                let angle = typst_number(*angle);
                let length_scale = typst_number(*length_scale);
                write!(
                    output,
                    "{indent}{}(absolute: {angle}deg, atom-sep: base-sep * {length_scale}",
                    bond_function_name(bond_type),
                )
                .unwrap();
                if let Some(offset) = offset {
                    write!(output, ", offset: \"{}\"", escape_string(offset)).unwrap();
                }
                if !name.is_empty() {
                    write!(output, ", name: \"{}\"", escape_string(name)).unwrap();
                }
                output.push_str(")\n");
            }
            Command::Branch { body } => {
                writeln!(output, "{indent}branch({{").unwrap();
                format_commands(output, body, depth + 1, indent_width);
                writeln!(output, "{indent}}})").unwrap();
            }
            Command::ComponentBreak => {
                writeln!(output, "{indent}operator(none, margin: base-sep * 0.5)").unwrap();
            }
        }
    }
}

fn format_atom_label(atom: &AtomLabel) -> String {
    let symbol = escape_content(&atom.symbol);
    let base = match atom.hydrogen_count {
        0 => format!("[{symbol}]"),
        1 => format!("[{symbol}H]"),
        count => format!("[{symbol}#math.attach([H], b: [{count}])]"),
    };
    let mut attachments = Vec::new();

    if let Some(isotope) = atom.isotope {
        attachments.push(format!("tl: [{isotope}]"));
    }

    let top_right = format!("{}{}", charge_text(atom.charge), radical_text(atom.radical));
    if !top_right.is_empty() {
        attachments.push(format!("tr: [{}]", escape_content(&top_right)));
    }

    if let Some(atom_map) = atom.atom_map {
        attachments.push(format!("br: [:{atom_map}]"));
    }

    if attachments.is_empty() {
        format!("math.equation(math.attach({base}))")
    } else {
        format!(
            "math.equation(math.attach({base}, {}))",
            attachments.join(", ")
        )
    }
}

fn charge_text(charge: i8) -> String {
    match charge {
        0 => String::new(),
        1 => "+".to_string(),
        -1 => "−".to_string(),
        charge if charge > 1 => format!("{charge}+"),
        charge => format!("{}−", charge.abs()),
    }
}

fn radical_text(radical: Option<u8>) -> String {
    match radical {
        None | Some(0) => String::new(),
        Some(1) => ":".to_string(),
        Some(2) => "•".to_string(),
        Some(3) => "••".to_string(),
        Some(radical) => format!("rad{radical}"),
    }
}

fn escape_content(value: &str) -> String {
    value
        .replace('\\', "\\\\")
        .replace('#', "\\#")
        .replace('*', "\\*")
        .replace(']', "\\]")
}

fn format_links(links: &[LinkData], depth: usize, indent_width: usize) -> String {
    if links.is_empty() {
        return String::new();
    }

    let indent = " ".repeat(depth * indent_width);
    let item_indent = " ".repeat((depth + 1) * indent_width);
    let mut output = String::from("links: (\n");
    for link in links {
        let angle = typst_number(link.angle);
        let length_scale = typst_number(link.length_scale);
        write!(
            output,
            "{item_indent}\"{}\": {}(absolute: {}deg, atom-sep: base-sep * {}",
            escape_string(&link.target),
            bond_function_name(&link.bond_type),
            angle,
            length_scale,
        )
        .unwrap();
        if let Some(offset) = &link.offset {
            write!(output, ", offset: \"{}\"", escape_string(offset)).unwrap();
        }
        if !link.name.is_empty() {
            write!(output, ", name: \"{}\"", escape_string(&link.name)).unwrap();
        }
        output.push_str("),\n");
    }
    write!(output, "{indent})").unwrap();
    output
}

fn has_extended_bonds(commands: &[Command]) -> bool {
    commands.iter().any(|command| match command {
        Command::Fragment { links, .. } => {
            links.iter().any(|link| is_extended_bond(&link.bond_type))
        }
        Command::Bond { bond_type, .. } => is_extended_bond(bond_type),
        Command::Branch { body } => has_extended_bonds(body),
        Command::ComponentBreak => false,
    })
}

fn is_extended_bond(bond_type: &str) -> bool {
    matches!(
        bond_type,
        "aromatic"
            | "single-or-double"
            | "single-or-aromatic"
            | "double-or-aromatic"
            | "any"
            | "either"
            | "coordination-right"
            | "coordination-left"
            | "hydrogen"
    )
}

fn bond_function_name(bond_type: &str) -> &str {
    match bond_type {
        "aromatic" => "_molchemist-aromatic",
        "single-or-double" => "_molchemist-single-or-double",
        "single-or-aromatic" => "_molchemist-single-or-aromatic",
        "double-or-aromatic" => "_molchemist-double-or-aromatic",
        "any" | "either" => "_molchemist-wavy",
        "coordination-right" => "_molchemist-coordination-right",
        "coordination-left" => "_molchemist-coordination-left",
        "hydrogen" => "_molchemist-hydrogen",
        bond_type => bond_type,
    }
}

fn typst_number(value: f64) -> String {
    let value = value.to_string();
    if let Some(value) = value.strip_prefix('-') {
        format!("−{value}")
    } else {
        value
    }
}

fn escape_string(value: &str) -> String {
    value
        .replace('\\', "\\\\")
        .replace('"', "\\\"")
        .replace('\n', "\\n")
        .replace('\r', "\\r")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn commands() -> Vec<Command> {
        vec![
            Command::Fragment {
                element: "O".to_string(),
                name: "a0".to_string(),
                links: Vec::new(),
                atom: None,
                annotation: None,
            },
            Command::Bond {
                name: "b0".to_string(),
                bond_type: "double".to_string(),
                angle: 90.0,
                offset: Some("right".to_string()),
                length_scale: 1.25,
            },
            Command::Branch {
                body: vec![Command::Fragment {
                    element: "C".to_string(),
                    name: "a1".to_string(),
                    links: vec![LinkData {
                        target: "a0".to_string(),
                        name: "b1".to_string(),
                        bond_type: "single".to_string(),
                        angle: 180.0,
                        offset: None,
                        length_scale: 1.0,
                    }],
                    atom: None,
                    annotation: None,
                }],
            },
        ]
    }

    #[test]
    fn raw_format_matches_typst_dump_shape() {
        assert_eq!(
            format_alchemist(&commands(), "3em", 2),
            concat!(
                "#let base-sep = 3em\n",
                "#skeletize({\n",
                "  fragment(\"O\", name: \"a0\")\n",
                "  double(absolute: 90deg, atom-sep: base-sep * 1.25, offset: \"right\", name: \"b0\")\n",
                "  branch({\n",
                "    fragment(\"C\", name: \"a1\", links: (\n",
                "      \"a0\": single(absolute: 180deg, atom-sep: base-sep * 1, name: \"b1\"),\n",
                "    ))\n",
                "  })\n",
                "})",
            )
        );
    }

    #[test]
    fn standalone_format_adds_only_document_wrapper() {
        let output = format_standalone(&commands(), "3em", 2, &StandaloneOptions::default());
        assert!(output.starts_with(
            "#import \"@preview/alchemist:0.2.0\": *\n\n#set page(width: auto, height: auto, margin: 3mm)\n\n"
        ));
        assert!(output.ends_with("})"));
    }

    #[test]
    fn standalone_code_wrapper_preserves_generated_code() {
        let code = "#skeletize({\n  fragment(\"O\")\n})";
        let output = format_standalone_code(code, &StandaloneOptions::default());
        assert!(output.ends_with(code));
    }

    #[test]
    fn structured_atom_metadata_uses_math_attachments() {
        let commands = vec![Command::Fragment {
            element: "CH_3^+".to_string(),
            name: "a0".to_string(),
            links: Vec::new(),
            atom: Some(AtomLabel {
                symbol: "C".to_string(),
                hydrogen_count: 3,
                charge: 1,
                isotope: Some(13),
                radical: Some(2),
                atom_map: Some(7),
            }),
            annotation: None,
        }];

        let output = format_alchemist(&commands, "3em", 2);

        assert!(output.contains(concat!(
            "fragment(math.equation(math.attach(",
            "[C#math.attach([H], b: [3])], ",
            "tl: [13], tr: [+•], br: [:7]",
            ")), name: \"a0\")",
        )));
    }

    #[test]
    fn component_breaks_reset_placement_without_a_visible_operator() {
        let commands = vec![
            Command::Fragment {
                element: "Na^+".to_string(),
                name: "a0".to_string(),
                links: Vec::new(),
                atom: None,
                annotation: None,
            },
            Command::ComponentBreak,
            Command::Fragment {
                element: "Cl^-".to_string(),
                name: "a1".to_string(),
                links: Vec::new(),
                atom: None,
                annotation: None,
            },
        ];

        assert_eq!(
            format_alchemist(&commands, "3em", 2),
            concat!(
                "#let base-sep = 3em\n",
                "#skeletize({\n",
                "  fragment(\"Na^+\", name: \"a0\")\n",
                "  operator(none, margin: base-sep * 0.5)\n",
                "  fragment(\"Cl^-\", name: \"a1\")\n",
                "})",
            )
        );
    }

    #[test]
    fn extended_bonds_emit_self_contained_typst_helpers() {
        let commands = vec![
            Command::Fragment {
                element: "C".to_string(),
                name: "a0".to_string(),
                links: vec![LinkData {
                    target: "a2".to_string(),
                    name: "b1".to_string(),
                    bond_type: "any".to_string(),
                    angle: 180.0,
                    offset: None,
                    length_scale: 1.0,
                }],
                atom: None,
                annotation: None,
            },
            Command::Bond {
                name: "b0".to_string(),
                bond_type: "aromatic".to_string(),
                angle: 0.0,
                offset: Some("left".to_string()),
                length_scale: 1.0,
            },
            Command::Branch {
                body: vec![Command::Bond {
                    name: "b2".to_string(),
                    bond_type: "coordination-left".to_string(),
                    angle: 90.0,
                    offset: None,
                    length_scale: 1.0,
                }],
            },
        ];

        let output = format_alchemist(&commands, "3em", 2);

        assert_eq!(output.matches("#import \"@preview/cetz:0.5.2\"").count(), 1);
        assert!(output.contains("#let _molchemist-wavy = build-link"));
        assert!(output.contains("\"a2\": _molchemist-wavy("));
        assert!(output.contains("_molchemist-aromatic(absolute: 0deg"));
        assert!(output.contains("_molchemist-coordination-left(absolute: 90deg"));
    }

    #[test]
    fn negative_numbers_match_typst_stringification() {
        assert_eq!(typst_number(-90.0), "−90");
        assert_eq!(typst_number(-0.25), "−0.25");
    }
}
