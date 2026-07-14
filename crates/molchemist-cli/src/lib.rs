//! WebAssembly-backed generation engine used by the `molchemist` executable.

mod runtime;

pub use runtime::Generator;

pub const DEFAULT_ALCHEMIST_IMPORT: &str = "@preview/alchemist:0.2.0";

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

pub fn format_standalone_code(code: &str, options: &StandaloneOptions) -> String {
    format!(
        "#import \"{}\": *\n\n#set page(width: auto, height: auto, margin: {})\n\n{}",
        escape_typst_string(&options.alchemist_import),
        options.page_margin,
        code,
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
        let source = "#skeletize({\n  fragment(\"O\")\n})";
        let document = format_standalone_code(source, &StandaloneOptions::default());
        assert!(document.starts_with(
            "#import \"@preview/alchemist:0.2.0\": *\n\n#set page(width: auto, height: auto, margin: 3mm)\n\n"
        ));
        assert!(document.ends_with(source));
    }
}
