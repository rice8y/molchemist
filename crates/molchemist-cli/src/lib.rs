//! WebAssembly-backed generation engine used by the `molchemist` executable.

mod runtime;

pub use runtime::Generator;

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
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct StandaloneOptions {
    pub page_margin: String,
}

impl Default for StandaloneOptions {
    fn default() -> Self {
        Self {
            page_margin: "3mm".to_string(),
        }
    }
}

pub fn format_standalone_code(code: &str, options: &StandaloneOptions) -> String {
    format!(
        "#set page(width: auto, height: auto, margin: {})\n\n{}",
        options.page_margin, code,
    )
}

fn escape_typst_string(value: &str) -> String {
    value
        .replace('\\', "\\\\")
        .replace('"', "\\\"")
        .replace('\n', "\\n")
        .replace('\r', "\\r")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn standalone_wrapper_preserves_generated_source() {
        let source = "#let _scene = ()\n#_render-graphic(_scene, 3em)";
        let document = format_standalone_code(source, &StandaloneOptions::default());
        assert!(document.starts_with("#set page(width: auto, height: auto, margin: 3mm)\n\n"));
        assert!(document.ends_with(source));
    }
}
