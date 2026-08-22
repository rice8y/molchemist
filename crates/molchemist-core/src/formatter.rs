use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write;

use crate::{AtomLabel, Command, CtfilePoint, LinkData};

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

#let _molchemist-crossed-double = build-link((length, ctx, cetz-ctx, args) => {
  import cetz.draw: *
  let gap = utils.convert-length(
    cetz-ctx,
    args.at("gap", default: ctx.config.double.gap),
  ) / 2
  let stroke = args.at("stroke", default: ctx.config.double.stroke)
  line(
    (0, -gap),
    (length, gap),
    stroke: args.at("stroke-right", default: stroke),
  )
  line(
    (0, gap),
    (length, -gap),
    stroke: args.at("stroke-left", default: stroke),
  )
})

#let _molchemist-quadruple = build-link((length, ctx, cetz-ctx, args) => {
  import cetz.draw: *
  let gap = utils.convert-length(
    cetz-ctx,
    args.at("gap", default: ctx.config.double.gap),
  )
  let stroke = args.at("stroke", default: ctx.config.double.stroke)
  for offset in (-1.5, -0.5, 0.5, 1.5) {
    line(
      (0, offset * gap),
      (length, offset * gap),
      stroke: stroke,
    )
  }
})

#let _molchemist-coordination-right = build-link((length, ctx, _, args) => {
  import cetz.draw: *
  line(
    (0, 0),
    (length, 0),
    stroke: args.at("stroke", default: ctx.config.single.stroke),
    mark: (end: ">", fill: black),
  )
})

#let _molchemist-coordination-left = build-link((length, ctx, _, args) => {
  import cetz.draw: *
  line(
    (0, 0),
    (length, 0),
    stroke: args.at("stroke", default: ctx.config.single.stroke),
    mark: (start: ">", fill: black),
  )
})

#let _molchemist-aromatic = _molchemist-partial-double("dashed")
#let _molchemist-single-or-double = _molchemist-partial-double("dotted")
#let _molchemist-single-or-aromatic = _molchemist-dashed-single("dashed")
#let _molchemist-double-or-aromatic = _molchemist-dashed-double
#let _molchemist-hydrogen = _molchemist-dashed-single("dotted")

"#;

const CTFILE_DEFINITIONS: &str = r#"
#let _molchemist-atom-center(centers, index) = centers.at(index)

#let _molchemist-highlight-path(bonds, atoms, centers, labels, paint) = {
  cetz.draw.get-ctx(cetz-ctx => {
    let atom-radius = utils.convert-length(cetz-ctx, base-sep * 0.28)
    let bond-radius = utils.convert-length(cetz-ctx, base-sep * 0.08)
    let label-padding = bond-radius
    cetz.draw.compound-path(fill: paint, stroke: none, fill-rule: "non-zero", {
      for (start-index, end-index) in bonds {
        let (_, start) = cetz.coordinate.resolve(cetz-ctx, _molchemist-atom-center(centers, start-index))
        let (_, end) = cetz.coordinate.resolve(cetz-ctx, _molchemist-atom-center(centers, end-index))
        let dx = end.at(0) - start.at(0)
        let dy = end.at(1) - start.at(1)
        let length = calc.sqrt(dx * dx + dy * dy)
        if length > 0 {
          let nx = -dy / length * bond-radius
          let ny = dx / length * bond-radius
          // Match the counter-clockwise CeTZ circle/rect winding so the
          // non-zero compound fill produces a union at atom-bond joins.
          cetz.draw.line(
            (start.at(0) + nx, start.at(1) + ny),
            (start.at(0) - nx, start.at(1) - ny),
            (end.at(0) - nx, end.at(1) - ny),
            (end.at(0) + nx, end.at(1) + ny),
            close: true,
          )
          cetz.draw.circle(start, radius: bond-radius)
          cetz.draw.circle(end, radius: bond-radius)
        }
      }
      for atom-index in atoms {
        let (_, center) = cetz.coordinate.resolve(cetz-ctx, _molchemist-atom-center(centers, atom-index))
        let label = labels.at(atom-index)
        if label == none {
          cetz.draw.circle(center, radius: atom-radius)
        } else {
          let prefix = "a" + str(atom-index) + ".0."
          let (_, west) = cetz.coordinate.resolve(cetz-ctx, (name: "molchemist-structure", anchor: prefix + "west"))
          let (_, east) = cetz.coordinate.resolve(cetz-ctx, (name: "molchemist-structure", anchor: prefix + "east"))
          let (_, north) = cetz.coordinate.resolve(cetz-ctx, (name: "molchemist-structure", anchor: prefix + "north"))
          let (_, south) = cetz.coordinate.resolve(cetz-ctx, (name: "molchemist-structure", anchor: prefix + "south"))
          let width = east.at(0) - west.at(0)
          let height = north.at(1) - south.at(1)
          if width <= atom-radius and height <= atom-radius {
            cetz.draw.circle(center, radius: atom-radius)
          } else {
            cetz.draw.rect(
              (west.at(0) - label-padding, south.at(1) - label-padding),
              (east.at(0) + label-padding, north.at(1) + label-padding),
              radius: label-padding,
            )
          }
        }
      }
    })
  })
}

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
    let extended_bonds = has_extended_bonds(commands);
    let has_ctfile = commands
        .iter()
        .any(|command| matches!(command, Command::Ctfile { .. }));
    if extended_bonds {
        output.push_str(EXTENDED_BOND_DEFINITIONS);
    } else if has_ctfile {
        output.push_str("#import \"@preview/cetz:0.5.2\"\n\n");
    }
    if has_ctfile {
        output.push_str(CTFILE_DEFINITIONS);
    }
    let annotations = collect_stereo_annotations(commands);
    if !annotations.is_empty() {
        output.push_str("#let _molchemist-structure = ");
    } else {
        output.push('#');
    }
    if has_ctfile {
        output.push_str("cetz.canvas({\n");
        format_ctfile_atom_metadata(&mut output, commands, 1, indent_width);
        output.push_str("  draw-skeleton(name: \"molchemist-structure\", {\n");
        format_commands(&mut output, commands, 2, indent_width);
        output.push_str("  })\n");
        format_ctfile_overlays(&mut output, commands, 1, indent_width);
        output.push_str("})");
    } else {
        output.push_str("skeletize({\n");
        format_commands(&mut output, commands, 1, indent_width);
        output.push_str("})");
    }
    if !annotations.is_empty() {
        output.push_str("\n#stack(\n");
        output.push_str("  dir: ttb,\n");
        output.push_str("  spacing: 0.45em,\n");
        output.push_str("  _molchemist-structure,\n");
        writeln!(
            output,
            "  text(size: 0.8em, fill: luma(80%))[Stereo annotations: {}],",
            escape_content(&annotations.join(", ")),
        )
        .unwrap();
        output.push(')');
    }
    output
}

fn collect_stereo_annotations(commands: &[Command]) -> Vec<String> {
    let mut annotations = Vec::new();
    for command in commands {
        match command {
            Command::Fragment {
                element,
                name,
                annotation: Some(annotation),
                ..
            } => {
                let label = if element.is_empty() { name } else { element };
                annotations.push(format!("{label} {annotation} ({name})"));
            }
            Command::Branch { body } => annotations.extend(collect_stereo_annotations(body)),
            _ => {}
        }
    }
    annotations
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

fn atom_command_index(name: &str) -> Option<usize> {
    name.strip_prefix('a')?.parse().ok()
}

fn fragment_body(element: &str, atom: Option<&AtomLabel>) -> String {
    atom.map_or_else(
        || format!("\"{}\"", escape_string(element)),
        format_atom_label,
    )
}

fn collect_ctfile_atom_metadata(
    commands: &[Command],
    current_atom: Option<usize>,
    atoms: &mut BTreeSet<usize>,
    centers: &mut BTreeMap<usize, String>,
    labels: &mut BTreeMap<usize, String>,
) {
    let mut current_atom = current_atom;
    let mut pending_bond: Option<&str> = None;

    for command in commands {
        match command {
            Command::Fragment {
                element,
                name,
                atom,
                ..
            } => {
                let Some(atom_index) = atom_command_index(name) else {
                    pending_bond = None;
                    continue;
                };
                atoms.insert(atom_index);
                if !element.is_empty() {
                    labels.insert(atom_index, fragment_body(element, atom.as_ref()));
                }
                if let (Some(previous_atom), Some(bond_name)) = (current_atom, pending_bond) {
                    centers
                        .entry(previous_atom)
                        .or_insert_with(|| format!("{bond_name}-start-anchor"));
                    centers
                        .entry(atom_index)
                        .or_insert_with(|| format!("{bond_name}-end-anchor"));
                }
                current_atom = Some(atom_index);
                pending_bond = None;
            }
            Command::Bond { name, .. } => pending_bond = Some(name),
            Command::Branch { body } => {
                collect_ctfile_atom_metadata(body, current_atom, atoms, centers, labels);
            }
            Command::ComponentBreak => {
                current_atom = None;
                pending_bond = None;
            }
            Command::Ctfile { .. } => {}
        }
    }
}

fn format_ctfile_atom_metadata(
    output: &mut String,
    commands: &[Command],
    depth: usize,
    indent_width: usize,
) {
    let mut atoms = BTreeSet::new();
    let mut centers = BTreeMap::new();
    let mut labels = BTreeMap::new();
    collect_ctfile_atom_metadata(commands, None, &mut atoms, &mut centers, &mut labels);
    let indent = " ".repeat(depth * indent_width);
    let item_indent = " ".repeat((depth + 1) * indent_width);

    writeln!(output, "{indent}let _molchemist-atom-centers = (").unwrap();
    for atom_index in &atoms {
        let anchor = centers
            .get(atom_index)
            .cloned()
            .unwrap_or_else(|| format!("a{atom_index}.0.mid"));
        writeln!(
            output,
            "{item_indent}(name: \"molchemist-structure\", anchor: \"{}\"),",
            escape_string(&anchor)
        )
        .unwrap();
    }
    writeln!(output, "{indent})").unwrap();
    writeln!(output, "{indent}let _molchemist-atom-labels = (").unwrap();
    for atom_index in atoms {
        let label = labels
            .get(&atom_index)
            .cloned()
            .unwrap_or_else(|| "none".to_string());
        writeln!(output, "{item_indent}{label},").unwrap();
    }
    writeln!(output, "{indent})").unwrap();
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

                    let body = fragment_body(element, atom.as_ref());
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
            Command::Ctfile { .. } => {}
        }
    }
}

fn format_ctfile_overlays(
    output: &mut String,
    commands: &[Command],
    depth: usize,
    indent_width: usize,
) {
    let indent = " ".repeat(depth * indent_width);
    for command in commands {
        let Command::Ctfile {
            sgroups,
            highlights,
            atom_queries,
            bond_queries,
            atom_annotations,
            variable_attachments,
        } = command
        else {
            continue;
        };
        writeln!(
            output,
            "{indent}let _molchemist-highlight = rgb(\"#ffd43b\").transparentize(45%)",
        )
        .unwrap();
        let highlight_bonds = highlights
            .iter()
            .flat_map(|highlight| &highlight.bonds)
            .map(|bond| format!("({}, {})", bond.atom1_index, bond.atom2_index))
            .collect::<Vec<_>>();
        let highlight_atoms = highlights
            .iter()
            .flat_map(|highlight| &highlight.atom_indexes)
            .map(ToString::to_string)
            .collect::<Vec<_>>();
        if !highlight_bonds.is_empty() || !highlight_atoms.is_empty() {
            writeln!(
                output,
                "{indent}cetz.draw.on-layer(−1, _molchemist-highlight-path({}, {}, _molchemist-atom-centers, _molchemist-atom-labels, _molchemist-highlight))",
                format_typst_array(&highlight_bonds),
                format_typst_array(&highlight_atoms),
            )
            .unwrap();
        }
        for group in sgroups {
            let name = format!("molchemist-sgroup-{}", group.id);
            let left_top = format_ctfile_point(&group.left_top);
            let left_bottom = format_ctfile_point(&group.left_bottom);
            let right_top = format_ctfile_point(&group.right_top);
            let right_bottom = format_ctfile_point(&group.right_bottom);
            writeln!(
                output,
                "{indent}cetz.draw.line({left_top}, {left_bottom}, name: \"{name}-left\", stroke: black)",
            )
            .unwrap();
            writeln!(
                output,
                "{indent}cetz.draw.line({left_top}, (rel: (base-sep * 0.18, 0pt), to: {left_top}), stroke: black)",
            )
            .unwrap();
            writeln!(
                output,
                "{indent}cetz.draw.line({left_bottom}, (rel: (base-sep * 0.18, 0pt), to: {left_bottom}), stroke: black)",
            )
            .unwrap();
            writeln!(
                output,
                "{indent}cetz.draw.line({right_top}, {right_bottom}, name: \"{name}-right\", stroke: black)",
            )
            .unwrap();
            writeln!(
                output,
                "{indent}cetz.draw.line({right_top}, (rel: (base-sep * −0.18, 0pt), to: {right_top}), stroke: black)",
            )
            .unwrap();
            writeln!(
                output,
                "{indent}cetz.draw.line({right_bottom}, (rel: (base-sep * −0.18, 0pt), to: {right_bottom}), stroke: black)",
            )
            .unwrap();
            if let Some(label) = &group.label {
                writeln!(
                    output,
                    "{indent}cetz.draw.content((rel: (base-sep * 0.01, base-sep * 0.035), to: {right_bottom}), text(size: 0.8em)[{}], anchor: \"north-west\")",
                    escape_content(label),
                )
                .unwrap();
            }
        }
        for query in atom_queries {
            if !query.compact_label.is_empty() {
                writeln!(
                    output,
                    "{indent}cetz.draw.content((rel: (0pt, base-sep * −0.36), to: _molchemist-atom-center(_molchemist-atom-centers, {})), text(size: 0.44em, fill: luma(38%))[{}], anchor: \"north\")",
                    query.atom_index,
                    escape_content(&query.compact_label),
                )
                .unwrap();
            }
        }
        for query in bond_queries {
            writeln!(
                output,
                "{indent}cetz.draw.content((rel: (0pt, base-sep * 0.2), to: (name: \"molchemist-structure\", anchor: \"b{}.50%\")), text(size: 0.56em, fill: luma(32%))[{}], anchor: \"center\")",
                query.bond_index,
                escape_content(&query.label),
            )
            .unwrap();
        }
        for attachment in variable_attachments {
            let dash = if attachment.mode.eq_ignore_ascii_case("ANY") {
                ", dash: \"dotted\""
            } else {
                ""
            };
            for atom_index in &attachment.endpoint_atom_indexes {
                writeln!(
                    output,
                    "{indent}cetz.draw.line((name: \"molchemist-structure\", anchor: \"b{}.50%\"), _molchemist-atom-center(_molchemist-atom-centers, {}), stroke: (thickness: 0.7pt, paint: luma(8%){}))",
                    attachment.bond_index,
                    atom_index,
                    dash,
                )
                .unwrap();
            }
        }
        for annotation in atom_annotations {
            let (x, y, size) = if annotation.kind == "link-node" {
                (0.0, 0.3, 0.82)
            } else {
                (0.22, 0.18, 0.82)
            };
            writeln!(
                output,
                "{indent}cetz.draw.content((rel: (base-sep * {x}, base-sep * {y}), to: _molchemist-atom-center(_molchemist-atom-centers, {})), text(size: {size}em, fill: luma(18%))[{}], anchor: \"center\")",
                annotation.atom_index,
                escape_content(&annotation.label),
            )
            .unwrap();
        }
    }
}

fn format_typst_array(items: &[String]) -> String {
    match items {
        [] => "()".to_string(),
        [item] => format!("({item},)"),
        items => format!("({})", items.join(", ")),
    }
}

fn format_ctfile_point(point: &CtfilePoint) -> String {
    format!(
        "(rel: (base-sep * {}, base-sep * {}), to: _molchemist-atom-center(_molchemist-atom-centers, {}))",
        typst_number(point.offset[0]),
        typst_number(point.offset[1]),
        point.atom_index,
    )
}

fn format_atom_label(atom: &AtomLabel) -> String {
    let symbol = escape_content(&atom.symbol);
    let base = if atom.query_symbol {
        format!(
            "text(size: 0.62em)[{}{symbol}]",
            if atom.query_negated { "!" } else { "" }
        )
    } else if let Some(rgroup_label) = atom.rgroup_label {
        format!("math.attach([R], b: [{rgroup_label}], t: std.hide([{rgroup_label}]))")
    } else if atom.hidden {
        format!("std.hide([{symbol}])")
    } else {
        match atom.hydrogen_count {
            0 => format!("[{symbol}]"),
            1 => format!("[{symbol}H]"),
            count => format!("[{symbol}#math.attach([H], b: [{count}], t: std.hide([{count}]))]"),
        }
    };
    let mut attachments = Vec::new();

    if let Some(isotope) = atom.isotope {
        attachments.push(format!("tl: [{isotope}]"));
        attachments.push(format!("bl: std.hide([{isotope}])"));
    }

    let top_right = format!("{}{}", charge_text(atom.charge), radical_text(atom.radical));
    if !top_right.is_empty() {
        attachments.push(format!("tr: [{}]", escape_content(&top_right)));
        attachments.push(format!("br: std.hide([{}])", escape_content(&top_right)));
    }

    let chemical_label = if attachments.is_empty() {
        format!("math.attach({base})")
    } else {
        format!("math.attach({base}, {})", attachments.join(", "))
    };

    if let Some(atom_map) = atom.atom_map {
        format!("math.equation(math.attach({chemical_label} + [:{atom_map}]))")
    } else {
        format!("math.equation({chemical_label})")
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
        .replace('[', "\\[")
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
        Command::ComponentBreak | Command::Ctfile { .. } => false,
    })
}

fn is_extended_bond(bond_type: &str) -> bool {
    matches!(
        bond_type,
        "aromatic"
            | "quadruple"
            | "single-or-double"
            | "single-or-aromatic"
            | "double-or-aromatic"
            | "any"
            | "either"
            | "crossed-double"
            | "coordination-right"
            | "coordination-left"
            | "hydrogen"
    )
}

fn bond_function_name(bond_type: &str) -> &str {
    match bond_type {
        "aromatic" => "_molchemist-aromatic",
        "quadruple" => "_molchemist-quadruple",
        "single-or-double" => "_molchemist-single-or-double",
        "single-or-aromatic" => "_molchemist-single-or-aromatic",
        "double-or-aromatic" => "_molchemist-double-or-aromatic",
        "any" | "either" => "_molchemist-wavy",
        "crossed-double" => "_molchemist-crossed-double",
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
    fn ctfile_highlight_subpaths_share_counter_clockwise_winding() {
        assert!(CTFILE_DEFINITIONS.contains(concat!(
            "(start.at(0) + nx, start.at(1) + ny),\n",
            "            (start.at(0) - nx, start.at(1) - ny),\n",
            "            (end.at(0) - nx, end.at(1) - ny),\n",
            "            (end.at(0) + nx, end.at(1) + ny),",
        )));
        assert!(CTFILE_DEFINITIONS.contains(concat!(
            "let bond-radius = utils.convert-length(cetz-ctx, base-sep * 0.08)\n",
            "    let label-padding = bond-radius",
        )));
        assert!(CTFILE_DEFINITIONS.contains("anchor: prefix + \"west\""));
        assert!(CTFILE_DEFINITIONS.contains("south.at(1) - label-padding"));
        assert!(!CTFILE_DEFINITIONS.contains("let label-size = measure(label)"));
    }

    #[test]
    fn ctfile_atom_metadata_keeps_hidden_centers_and_visible_glyphs_separate() {
        let commands = vec![
            Command::Fragment {
                element: String::new(),
                name: "a0".to_string(),
                links: Vec::new(),
                atom: None,
                annotation: None,
            },
            Command::Bond {
                bond_type: "single".to_string(),
                angle: 0.0,
                length_scale: 1.0,
                offset: None,
                name: "b0".to_string(),
            },
            Command::Fragment {
                element: "[C,N,O,S,P]".to_string(),
                name: "a1".to_string(),
                links: Vec::new(),
                atom: None,
                annotation: None,
            },
        ];
        let mut atoms = BTreeSet::new();
        let mut centers = BTreeMap::new();
        let mut labels = BTreeMap::new();

        collect_ctfile_atom_metadata(&commands, None, &mut atoms, &mut centers, &mut labels);

        assert_eq!(atoms, BTreeSet::from([0, 1]));
        assert_eq!(centers.get(&0).unwrap(), "b0-start-anchor");
        assert_eq!(centers.get(&1).unwrap(), "b0-end-anchor");
        assert!(!labels.contains_key(&0));
        assert_eq!(labels.get(&1).unwrap(), "\"[C,N,O,S,P]\"");
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
                rgroup_label: None,
                query_symbol: false,
                query_negated: false,
                hidden: false,
            }),
            annotation: None,
        }];

        let output = format_alchemist(&commands, "3em", 2);

        assert!(output.contains(concat!(
            "fragment(math.equation(math.attach(math.attach(",
            "[C#math.attach([H], b: [3], t: std.hide([3]))], ",
            "tl: [13], bl: std.hide([13]), tr: [+•], br: std.hide([+•])",
            ") + [:7])), name: \"a0\")",
        )));
        assert!(!output.contains("br: [:7]"));
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
        assert!(output.contains("mark: (end: \">\", fill: black)"));
        assert!(output.contains("mark: (start: \">\", fill: black)"));
    }

    #[test]
    fn crossed_double_emits_its_self_contained_typst_helper() {
        let commands = vec![Command::Bond {
            name: "b0".to_string(),
            bond_type: "crossed-double".to_string(),
            angle: 0.0,
            offset: None,
            length_scale: 1.0,
        }];

        let output = format_alchemist(&commands, "3em", 2);

        assert!(output.contains("#let _molchemist-crossed-double = build-link"));
        assert!(output.contains("_molchemist-crossed-double(absolute: 0deg"));
    }

    #[test]
    fn quadruple_emits_its_self_contained_typst_helper() {
        let commands = vec![Command::Bond {
            name: "b0".to_string(),
            bond_type: "quadruple".to_string(),
            angle: 0.0,
            offset: None,
            length_scale: 1.0,
        }];

        let output = format_alchemist(&commands, "3em", 2);

        assert!(output.contains("#let _molchemist-quadruple = build-link"));
        assert!(output.contains("_molchemist-quadruple(absolute: 0deg"));
    }

    #[test]
    fn stereo_annotations_survive_formatted_output() {
        let commands = vec![Command::Fragment {
            element: "Pt".to_string(),
            name: "a0".to_string(),
            links: Vec::new(),
            atom: None,
            annotation: Some("@SP1".to_string()),
        }];

        let output = format_alchemist(&commands, "3em", 2);

        assert!(output.contains("#let _molchemist-structure = skeletize({"));
        assert!(output.contains("Stereo annotations: Pt @SP1 (a0)"));
        assert!(output.ends_with(')'));
    }

    #[test]
    fn negative_numbers_match_typst_stringification() {
        assert_eq!(typst_number(-90.0), "−90");
        assert_eq!(typst_number(-0.25), "−0.25");
    }
}
