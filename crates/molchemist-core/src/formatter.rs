//! Standalone Typst source assembled from the package's drawing modules.
//!
//! Local imports are resolved by the ordered module list. The generated source
//! retains one Alchemist import and contains no workspace-relative file paths.

const RENDER_MODULES: &[(&str, &str)] = &[
    (
        "anchors.typ",
        include_str!("../../../package/src/anchors.typ"),
    ),
    ("bonds.typ", include_str!("../../../package/src/bonds.typ")),
    (
        "labels.typ",
        include_str!("../../../package/src/labels.typ"),
    ),
    (
        "highlights.typ",
        include_str!("../../../package/src/highlights.typ"),
    ),
    (
        "annotations.typ",
        include_str!("../../../package/src/annotations.typ"),
    ),
    (
        "annotation-renderer.typ",
        include_str!("../../../package/src/annotation-renderer.typ"),
    ),
    (
        "indices.typ",
        include_str!("../../../package/src/indices.typ"),
    ),
    (
        "ctfile.typ",
        include_str!("../../../package/src/ctfile.typ"),
    ),
    (
        "layout.typ",
        include_str!("../../../package/src/layout.typ"),
    ),
    (
        "render.typ",
        include_str!("../../../package/src/render.typ"),
    ),
];

/// Shared reaction and R-group panel composition for standalone documents.
pub fn composition_code() -> &'static str {
    include_str!("../../../package/src/composition.typ")
}

/// Serialize a scene and configuration with the same renderer as the package.
/// A text configuration is a Typst expression supplied by the Typst caller;
/// other configuration values are serialized from CBOR by this function.
pub fn format_coordinate_code(ast: &[u8], config: &[u8], base_sep: &str) -> Result<String, String> {
    let ast: ciborium::Value = ciborium::from_reader(ast).map_err(|e| e.to_string())?;
    let config: ciborium::Value = ciborium::from_reader(config).map_err(|e| e.to_string())?;
    let config = if let ciborium::Value::Text(source) = config {
        source
    } else {
        typst_value(&config)?
    };
    let mut code = String::from("#import \"@preview/alchemist:0.2.0\": *\n\n");
    for (name, module) in RENDER_MODULES {
        for line in module.lines() {
            if let Some(import) = line.strip_prefix("#import \"") {
                let (path, _) = import
                    .split_once('"')
                    .ok_or_else(|| format!("Invalid import in {name}"))?;
                if path == "@preview/alchemist:0.2.0" {
                    continue;
                }
                if !path.starts_with('@') {
                    if !RENDER_MODULES.iter().any(|(name, _)| *name == path) {
                        return Err(format!(
                            "Renderer module {name} imports unbundled module {path}"
                        ));
                    }
                    continue;
                }
            }
            code.push_str(line);
            code.push('\n');
        }
        code.push('\n');
    }
    code.push_str(&format!("#let _scene = {}\n#let _scene-config = {}\n#_render-with-stereo-annotations(_scene, _render-graphic(_scene, {}, config: _scene-config))", typst_value(&ast)?, config, base_sep));
    Ok(code)
}

fn typst_value(value: &ciborium::Value) -> Result<String, String> {
    use ciborium::Value;
    Ok(match value {
        Value::Null => "none".into(),
        Value::Bool(b) => b.to_string(),
        Value::Integer(n) => i128::from(*n).to_string(),
        Value::Float(n) if n.is_finite() => typst_number(*n),
        Value::Text(s) => format!("\"{}\"", escape_string(s)),
        Value::Array(a) => {
            if a.is_empty() {
                "()".into()
            } else {
                format!(
                    "({},)",
                    a.iter()
                        .map(typst_value)
                        .collect::<Result<Vec<_>, _>>()?
                        .join(", ")
                )
            }
        }
        Value::Map(m) => {
            let mut items = m
                .iter()
                .map(|(k, v)| Ok((typst_value(k)?, typst_value(v)?)))
                .collect::<Result<Vec<_>, String>>()?;
            items.sort_by(|a, b| a.0.cmp(&b.0));
            if items.is_empty() {
                "(:)".into()
            } else {
                format!(
                    "({},)",
                    items
                        .into_iter()
                        .map(|(k, v)| format!("{k}: {v}"))
                        .collect::<Vec<_>>()
                        .join(", ")
                )
            }
        }
        _ => return Err("Cannot export this value to Typst source".into()),
    })
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
    use ciborium::Value;

    #[test]
    fn source_values_preserve_text_and_reject_nonfinite_numbers() {
        assert_eq!(
            typst_value(&Value::Text("C\"\\\n".into())).unwrap(),
            "\"C\\\"\\\\\\n\""
        );
        assert!(typst_value(&Value::Float(f64::INFINITY)).is_err());
        assert!(typst_value(&Value::Float(f64::NAN)).is_err());
        assert_eq!(typst_value(&Value::Float(-0.25)).unwrap(), "−0.25");
    }
}
