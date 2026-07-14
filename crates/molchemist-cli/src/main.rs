use std::fs;
use std::io::{self, Read, Write};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::{Args, Parser, Subcommand, ValueEnum};
use molchemist_cli::{
    format_standalone_code, Generator, RenderMode, StandaloneOptions, DEFAULT_ALCHEMIST_IMPORT,
};

#[derive(Parser)]
#[command(name = "molchemist", version, about)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Dump Alchemist Typst code for one molecule.
    Dump(DumpArgs),
}

#[derive(Args)]
struct DumpArgs {
    /// Input file. Use '-' or omit it to read standard input.
    #[arg(value_name = "INPUT", conflicts_with_all = ["text", "smiles"])]
    input: Option<PathBuf>,

    /// Read Molfile/SDF or SMILES content directly from this argument.
    #[arg(long, value_name = "TEXT", conflicts_with_all = ["input", "smiles"])]
    text: Option<String>,

    /// Read a SMILES string directly from this argument.
    #[arg(long, value_name = "SMILES", conflicts_with_all = ["input", "text"])]
    smiles: Option<String>,

    /// Input format. Auto uses the explicit input kind, file extension, then content.
    #[arg(short, long, value_enum, default_value_t = InputFormat::Auto)]
    format: InputFormat,

    /// Rendering mode, matching render-mol and render-smiles.
    #[arg(short, long, value_enum, default_value_t = Mode::Full)]
    mode: Mode,

    /// One-based record number for a multi-record SDF input.
    #[arg(long, default_value = "1")]
    record: NonZeroUsize,

    /// Write generated code to this file instead of standard output.
    #[arg(short, long, value_name = "PATH")]
    output: Option<PathBuf>,

    /// Generate a complete, directly compilable Typst document.
    #[arg(long)]
    standalone: bool,

    /// Alchemist atom separation used by the generated code.
    #[arg(long, default_value = "3em", value_parser = parse_typst_length)]
    atom_sep: String,

    /// Standalone document page margin.
    #[arg(long, default_value = "3mm", requires = "standalone", value_parser = parse_typst_length)]
    page_margin: String,

    /// Alchemist package import used by standalone documents.
    #[arg(long, default_value = DEFAULT_ALCHEMIST_IMPORT, requires = "standalone", value_parser = parse_import)]
    alchemist_import: String,

    /// Number of spaces used for each indentation level.
    #[arg(long, default_value_t = 2, value_parser = parse_indent)]
    indent: usize,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, ValueEnum)]
enum InputFormat {
    #[default]
    Auto,
    Mol,
    Sdf,
    Smiles,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, ValueEnum)]
enum Mode {
    #[default]
    Full,
    Abbreviate,
    Skeletal,
}

impl From<Mode> for RenderMode {
    fn from(mode: Mode) -> Self {
        match mode {
            Mode::Full => Self::Full,
            Mode::Abbreviate => Self::Abbreviate,
            Mode::Skeletal => Self::Skeletal,
        }
    }
}

struct Input {
    content: String,
    path: Option<PathBuf>,
    explicit_smiles: bool,
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error}");
            ExitCode::FAILURE
        }
    }
}

fn run(cli: Cli) -> Result<(), String> {
    match cli.command {
        Command::Dump(args) => dump(args),
    }
}

fn dump(args: DumpArgs) -> Result<(), String> {
    let input = read_input(&args)?;
    let format = detect_format(&input, args.format)?;
    let mode = args.mode.into();
    let mut generator = Generator::new()
        .map_err(|error| format!("could not initialize the conversion engine: {error}"))?;

    let generated = match format {
        InputFormat::Mol => {
            require_first_record(args.record, "Molfile")?;
            generator.sdf_to_code(&input.content, mode, &args.atom_sep, args.indent)
        }
        InputFormat::Sdf => {
            let record = select_sdf_record(&input.content, args.record.get())?;
            generator.sdf_to_code(record, mode, &args.atom_sep, args.indent)
        }
        InputFormat::Smiles => {
            require_first_record(args.record, "SMILES")?;
            let smiles = input.content.trim();
            if smiles.is_empty() {
                return Err("SMILES input is empty".to_string());
            }
            generator.smiles_to_code(smiles, mode, &args.atom_sep, args.indent)
        }
        InputFormat::Auto => unreachable!("auto format is resolved before conversion"),
    }
    .map_err(|error| format!("failed to convert {format}: {error}"))?;

    let generated = if args.standalone {
        format_standalone_code(
            &generated,
            &StandaloneOptions {
                alchemist_import: args.alchemist_import,
                page_margin: args.page_margin,
            },
        )
    } else {
        generated
    };

    write_output(args.output.as_deref(), generated.as_bytes())
}

fn read_input(args: &DumpArgs) -> Result<Input, String> {
    if let Some(smiles) = &args.smiles {
        return Ok(Input {
            content: smiles.clone(),
            path: None,
            explicit_smiles: true,
        });
    }

    if let Some(text) = &args.text {
        return Ok(Input {
            content: text.clone(),
            path: None,
            explicit_smiles: false,
        });
    }

    if let Some(path) = &args.input {
        if path != Path::new("-") {
            let content = fs::read_to_string(path)
                .map_err(|error| format!("could not read {}: {error}", path.display()))?;
            return Ok(Input {
                content,
                path: Some(path.clone()),
                explicit_smiles: false,
            });
        }
    }

    let mut content = String::new();
    io::stdin()
        .read_to_string(&mut content)
        .map_err(|error| format!("could not read standard input: {error}"))?;
    Ok(Input {
        content,
        path: None,
        explicit_smiles: false,
    })
}

fn detect_format(input: &Input, requested: InputFormat) -> Result<InputFormat, String> {
    if input.explicit_smiles {
        return match requested {
            InputFormat::Auto | InputFormat::Smiles => Ok(InputFormat::Smiles),
            _ => Err("--smiles conflicts with a non-SMILES --format".to_string()),
        };
    }

    if requested != InputFormat::Auto {
        return Ok(requested);
    }

    if let Some(extension) = input
        .path
        .as_deref()
        .and_then(Path::extension)
        .and_then(|extension| extension.to_str())
        .map(str::to_ascii_lowercase)
    {
        match extension.as_str() {
            "mol" => return Ok(InputFormat::Mol),
            "sdf" => return Ok(InputFormat::Sdf),
            "smi" | "smiles" => return Ok(InputFormat::Smiles),
            _ => {}
        }
    }

    if input.content.contains("M  END")
        && (input.content.contains("V2000") || input.content.contains("V3000"))
    {
        return Ok(if input.content.contains("$$$$") {
            InputFormat::Sdf
        } else {
            InputFormat::Mol
        });
    }

    let non_empty_lines = input
        .content
        .lines()
        .filter(|line| !line.trim().is_empty())
        .count();
    if non_empty_lines == 1 {
        return Ok(InputFormat::Smiles);
    }

    Err("could not detect input format; pass --format mol, sdf, or smiles".to_string())
}

fn select_sdf_record(content: &str, record: usize) -> Result<&str, String> {
    let records = content
        .split("$$$$")
        .filter(|part| !part.trim().is_empty())
        .map(|part| {
            part.strip_prefix("\r\n")
                .or_else(|| part.strip_prefix('\n'))
                .unwrap_or(part)
        })
        .collect::<Vec<_>>();
    records.get(record - 1).copied().ok_or_else(|| {
        format!(
            "SDF record {record} does not exist; input contains {} record(s)",
            records.len()
        )
    })
}

fn require_first_record(record: NonZeroUsize, format: &str) -> Result<(), String> {
    if record.get() == 1 {
        Ok(())
    } else {
        Err(format!(
            "--record can only exceed 1 for SDF input, not {format}"
        ))
    }
}

fn write_output(path: Option<&Path>, content: &[u8]) -> Result<(), String> {
    if let Some(path) = path {
        fs::write(path, content)
            .map_err(|error| format!("could not write {}: {error}", path.display()))
    } else {
        let mut stdout = io::stdout().lock();
        stdout
            .write_all(content)
            .and_then(|()| stdout.flush())
            .map_err(|error| format!("could not write standard output: {error}"))
    }
}

fn parse_typst_length(value: &str) -> Result<String, String> {
    const UNITS: &[&str] = &["pt", "mm", "cm", "in", "em"];
    let unit = UNITS
        .iter()
        .find(|unit| value.ends_with(**unit))
        .ok_or_else(|| "expected a Typst length using pt, mm, cm, in, or em".to_string())?;
    let number = &value[..value.len() - unit.len()];
    let parsed = number
        .parse::<f64>()
        .map_err(|_| "expected a numeric Typst length such as 3em or 2.5mm".to_string())?;
    if !parsed.is_finite() || parsed < 0.0 {
        return Err("Typst length must be a finite, non-negative value".to_string());
    }
    Ok(value.to_string())
}

fn parse_import(value: &str) -> Result<String, String> {
    if value.is_empty() || value.contains(['"', '\n', '\r']) {
        Err("alchemist import must be a non-empty package or file specifier".to_string())
    } else {
        Ok(value.to_string())
    }
}

fn parse_indent(value: &str) -> Result<usize, String> {
    let width = value
        .parse::<usize>()
        .map_err(|_| "indent must be an integer from 1 through 8".to_string())?;
    if (1..=8).contains(&width) {
        Ok(width)
    } else {
        Err("indent must be an integer from 1 through 8".to_string())
    }
}

impl std::fmt::Display for InputFormat {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter.write_str(match self {
            Self::Auto => "auto input",
            Self::Mol => "Molfile input",
            Self::Sdf => "SDF input",
            Self::Smiles => "SMILES input",
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn content_detection_distinguishes_molfile_and_smiles() {
        let mol = Input {
            content: "example\n  molchemist\n\n  0  0  0  0  0  0  0  0  0  0  0 V2000\nM  END\n"
                .to_string(),
            path: None,
            explicit_smiles: false,
        };
        assert_eq!(
            detect_format(&mol, InputFormat::Auto).unwrap(),
            InputFormat::Mol
        );

        let smiles = Input {
            content: "c1ccccc1\n".to_string(),
            path: None,
            explicit_smiles: false,
        };
        assert_eq!(
            detect_format(&smiles, InputFormat::Auto).unwrap(),
            InputFormat::Smiles
        );
    }

    #[test]
    fn selects_one_based_sdf_record() {
        let input = "first\n$$$$\nsecond\n$$$$\n";
        assert_eq!(select_sdf_record(input, 2).unwrap().trim(), "second");
        assert!(select_sdf_record(input, 3).is_err());
    }

    #[test]
    fn validates_typst_lengths() {
        assert_eq!(parse_typst_length("2.5mm").unwrap(), "2.5mm");
        assert!(parse_typst_length("calc(2mm)").is_err());
        assert!(parse_typst_length("-1em").is_err());
    }
}
