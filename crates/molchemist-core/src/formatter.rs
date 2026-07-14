use std::fmt::Write;

use crate::{Command, LinkData};

pub const DEFAULT_ALCHEMIST_IMPORT: &str = "@preview/alchemist:0.2.0";

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
    let mut output = format!("#let base-sep = {base_sep}\n#skeletize({{\n");
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

                    write!(output, "{indent}fragment(\"{}\"", escape_string(element)).unwrap();
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
                    "{indent}{bond_type}(absolute: {angle}deg, atom-sep: base-sep * {length_scale}"
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
        }
    }
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
            link.bond_type,
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
    fn negative_numbers_match_typst_stringification() {
        assert_eq!(typst_number(-90.0), "−90");
        assert_eq!(typst_number(-0.25), "−0.25");
    }
}
