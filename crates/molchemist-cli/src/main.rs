use std::fs;
use std::io::{self, Read, Write};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::{Args, Parser, Subcommand, ValueEnum};
use molchemist_cli::{format_standalone_code, Generator, RenderMode, StandaloneOptions};

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
    /// Inspect one Molfile/SDF record as loss-aware semantic JSON.
    Inspect(InspectArgs),
}

#[derive(Args)]
struct DumpArgs {
    /// Expand all known superatoms instead of contracting their labels.
    #[arg(long)]
    expand_superatoms: bool,
    /// Infer tetrahedral wedge orientation from nondegenerate 3D coordinates.
    #[arg(long)]
    infer_stereo: bool,
    /// Coordinate policy. Avoid enlarges spacing using the actual font metrics.
    #[arg(long, default_value = "avoid", value_parser = ["coordinates", "avoid", "reflow"])]
    layout: String,

    /// Preserve source component positions or pack disconnected components.
    #[arg(long, default_value = "pack", value_parser = ["pack", "preserve"])]
    components: String,

    /// Reject unresolved label collisions at Typst compilation time.
    #[arg(long)]
    strict_collisions: bool,

    /// Infer missing reaction atom correspondences; ambiguity is reported.
    #[arg(long)]
    infer_mapping: bool,

    /// Maximum number of states visited by the atom-correspondence search.
    #[arg(long, default_value = "200000")]
    mapping_search_limit: NonZeroUsize,

    /// Text placed below the reaction arrow.
    #[arg(long, default_value = "")]
    conditions: String,
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

    /// Handling of parsed features that the current renderer cannot depict faithfully.
    #[arg(long, value_enum, default_value_t = FidelityMode::Warn)]
    fidelity: FidelityMode,
}

#[derive(Args)]
struct InspectArgs {
    #[arg(long)]
    infer_mapping: bool,
    #[arg(long, default_value = "200000")]
    mapping_search_limit: NonZeroUsize,
    /// Molfile/SDF input file. Use '-' or omit it to read standard input.
    #[arg(value_name = "INPUT", conflicts_with = "text")]
    input: Option<PathBuf>,

    /// Read Molfile/SDF content directly from this argument.
    #[arg(long, value_name = "TEXT", conflicts_with = "input")]
    text: Option<String>,

    /// Input format. Auto uses the file extension, then content.
    #[arg(short, long, value_enum, default_value_t = InputFormat::Auto)]
    format: InputFormat,

    /// One-based record number for a multi-record SDF input.
    #[arg(long, default_value = "1")]
    record: NonZeroUsize,

    /// Write JSON to this file instead of standard output.
    #[arg(short, long, value_name = "PATH")]
    output: Option<PathBuf>,

    /// Emit compact JSON instead of indented JSON.
    #[arg(long)]
    compact: bool,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, ValueEnum)]
enum InputFormat {
    #[default]
    Auto,
    Mol,
    Sdf,
    Smiles,
    Mol2,
    Rxn,
    ReactionSmiles,
    Rgroup,
}

impl InputFormat {
    fn key(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::Mol => "mol",
            Self::Sdf => "sdf",
            Self::Smiles => "smiles",
            Self::Mol2 => "mol2",
            Self::Rxn => "rxn",
            Self::ReactionSmiles => "reaction-smiles",
            Self::Rgroup => "rgroup",
        }
    }
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, ValueEnum)]
enum Mode {
    #[default]
    Full,
    Abbreviate,
    Skeletal,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, ValueEnum)]
enum FidelityMode {
    Ignore,
    #[default]
    Warn,
    Strict,
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
        Command::Inspect(args) => inspect(args),
    }
}

fn dump(args: DumpArgs) -> Result<(), String> {
    let mut input = read_input(&args)?;
    let mut format = detect_format(&input, args.format)?;
    let mode = args.mode.into();
    let mut generator = Generator::new()
        .map_err(|error| format!("could not initialize the conversion engine: {error}"))?;

    let config = serde_json::json!({"layout": args.layout, "infer-stereo": args.infer_stereo, "components": args.components, "collision-policy": if args.strict_collisions { "error" } else { "report" }});
    if format == InputFormat::Rgroup {
        require_first_record(args.record, "RGfile")?;
        let document = generator.inspect_rgroup(&input.content)?;
        if args.fidelity != FidelityMode::Ignore {
            let root = document["root"].as_str().ok_or("Missing R-group root")?;
            apply_fidelity_policy(
                &generator.inspect_sdf_record_json(root, 1, false)?,
                args.fidelity,
            )?;
            for member in document["members"]
                .as_array()
                .ok_or("Missing R-group members")?
            {
                let data = member["data"].as_str().ok_or("Missing R-group member")?;
                let mut inspection: serde_json::Value =
                    serde_json::from_str(&generator.inspect_sdf_record_json(data, 1, false)?)
                        .map_err(|e| e.to_string())?;
                inspection["diagnostics"]
                    .as_array_mut()
                    .ok_or("Missing diagnostics")?
                    .retain(|d| d["code"] != "rgroup-attachment-context-not-supported");
                apply_fidelity_policy(&inspection.to_string(), args.fidelity)?;
            }
        }
        let style = config.clone();
        let code = generator.rgroup_to_code(&input.content, mode, &args.atom_sep, &style)?;
        let code = if args.standalone {
            format_standalone_code(
                &code,
                &StandaloneOptions {
                    page_margin: args.page_margin,
                },
            )
        } else {
            code
        };
        return write_output(args.output.as_deref(), code.as_bytes());
    }
    if matches!(format, InputFormat::Rxn | InputFormat::ReactionSmiles) {
        require_first_record(args.record, "reaction")?;
        let analysis = generator.inspect_reaction(
            &input.content,
            format.key(),
            args.infer_mapping,
            args.mapping_search_limit.get(),
        )?;
        if args.fidelity != FidelityMode::Ignore {
            for side in ["reactants", "agents", "products"] {
                for molecule in analysis["reaction"][side]
                    .as_array()
                    .ok_or("Missing reaction side")?
                {
                    if molecule["format"] == "mol" {
                        let data = molecule["data"]
                            .as_str()
                            .ok_or("Missing reaction molecule")?;
                        apply_fidelity_policy(
                            &generator.inspect_sdf_record_json(data, 1, false)?,
                            args.fidelity,
                        )?;
                    }
                }
            }
        }
        if analysis["ambiguous"] == true || analysis["searchComplete"] == false {
            let message = "Reaction mapping is ambiguous or the search limit was reached; inspect the mapping before relying on its reaction centers";
            if args.fidelity == FidelityMode::Strict {
                return Err(message.into());
            }
            if args.fidelity == FidelityMode::Warn {
                eprintln!("warning[mapping]: {message}");
            }
        }
        let style = config.clone();
        let code = generator.reaction_to_code(
            &analysis,
            mode,
            &args.atom_sep,
            &style,
            &args.conditions,
        )?;
        let code = if args.standalone {
            format_standalone_code(
                &code,
                &StandaloneOptions {
                    page_margin: args.page_margin,
                },
            )
        } else {
            code
        };
        return write_output(args.output.as_deref(), code.as_bytes());
    }
    let mut record = if format == InputFormat::Mol2 {
        input.content = generator.mol2_to_sdf(&input.content, args.record.get())?;
        format = InputFormat::Sdf;
        1
    } else {
        args.record.get()
    };
    if args.expand_superatoms {
        if !matches!(format, InputFormat::Mol | InputFormat::Sdf) {
            return Err("--expand-superatoms requires CTfile or MOL2 input".into());
        }
        input.content = generator.expand_superatoms(&input.content, record)?;
        record = 1;
    }

    if args.fidelity != FidelityMode::Ignore
        && matches!(format, InputFormat::Mol | InputFormat::Sdf)
    {
        let inspection = generator
            .inspect_sdf_record_json(&input.content, record, false)
            .map_err(|error| format!("failed to inspect {format}: {error}"))?;
        apply_fidelity_policy(&inspection, args.fidelity)?;
    }

    if format != InputFormat::Sdf {
        require_first_record(args.record, &format.to_string())?;
    }
    if format == InputFormat::Smiles {
        input.content = input.content.trim().to_string();
        if input.content.is_empty() {
            return Err("SMILES input is empty".into());
        }
    }
    let generated = generator
        .coordinate_to_code(
            &input.content,
            format.key(),
            mode,
            record,
            &args.atom_sep,
            &config,
        )
        .map_err(|error| format!("failed to convert {format}: {error}"))?;

    let generated = if args.standalone {
        format_standalone_code(
            &generated,
            &StandaloneOptions {
                page_margin: args.page_margin,
            },
        )
    } else {
        generated
    };

    write_output(args.output.as_deref(), generated.as_bytes())
}

fn inspect(args: InspectArgs) -> Result<(), String> {
    let input = read_inspect_input(&args)?;
    let format = detect_format(&input, args.format)?;
    let mut generator = Generator::new().map_err(|e| e.to_string())?;
    if matches!(format, InputFormat::Rxn | InputFormat::ReactionSmiles) {
        require_first_record(args.record, "reaction")?;
        let analysis = generator.inspect_reaction(
            &input.content,
            format.key(),
            args.infer_mapping,
            args.mapping_search_limit.get(),
        )?;
        let json = if args.compact {
            serde_json::to_string(&analysis)
        } else {
            serde_json::to_string_pretty(&analysis)
        }
        .map_err(|e| e.to_string())?;
        return write_output(args.output.as_deref(), format!("{json}\n").as_bytes());
    }
    match format {
        InputFormat::Rgroup => {
            require_first_record(args.record, "RGfile")?;
            let doc = generator.inspect_rgroup(&input.content)?;
            let json = if args.compact {
                serde_json::to_string(&doc)
            } else {
                serde_json::to_string_pretty(&doc)
            }
            .map_err(|e| e.to_string())?;
            return write_output(args.output.as_deref(), format!("{json}\n").as_bytes());
        }
        InputFormat::Mol => require_first_record(args.record, "Molfile")?,
        InputFormat::Sdf => {}
        InputFormat::Mol2 => {
            let sdf = generator.mol2_to_sdf(&input.content, args.record.get())?;
            let json = generator.inspect_sdf_record_json(&sdf, 1, !args.compact)?;
            return write_output(args.output.as_deref(), format!("{json}\n").as_bytes());
        }
        InputFormat::Rxn | InputFormat::ReactionSmiles => unreachable!(),
        InputFormat::Smiles => {
            return Err("inspect currently accepts Molfile or SDF input, not SMILES".to_string())
        }
        InputFormat::Auto => unreachable!("auto format is resolved before inspection"),
    }

    let mut generator = Generator::new()
        .map_err(|error| format!("could not initialize the inspection engine: {error}"))?;
    let mut json = generator
        .inspect_sdf_record_json(&input.content, args.record.get(), !args.compact)
        .map_err(|error| format!("failed to inspect {format}: {error}"))?;
    json.push('\n');
    write_output(args.output.as_deref(), json.as_bytes())
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

fn read_inspect_input(args: &InspectArgs) -> Result<Input, String> {
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

fn apply_fidelity_policy(inspection: &str, policy: FidelityMode) -> Result<(), String> {
    let value: serde_json::Value = serde_json::from_str(inspection)
        .map_err(|error| format!("inspection engine returned invalid JSON: {error}"))?;
    let diagnostics = value
        .get("diagnostics")
        .and_then(serde_json::Value::as_array)
        .ok_or_else(|| "inspection engine omitted diagnostics".to_string())?;
    if diagnostics.is_empty() || policy == FidelityMode::Ignore {
        return Ok(());
    }

    let messages = diagnostics
        .iter()
        .map(|diagnostic| {
            let code = diagnostic
                .get("code")
                .and_then(serde_json::Value::as_str)
                .unwrap_or("fidelity");
            let message = diagnostic
                .get("message")
                .and_then(serde_json::Value::as_str)
                .unwrap_or("the input contains a feature that is not depicted faithfully");
            (code, message)
        })
        .collect::<Vec<_>>();

    match policy {
        FidelityMode::Ignore => Ok(()),
        FidelityMode::Warn => {
            for (code, message) in messages {
                eprintln!("warning[{code}]: {message}");
            }
            Ok(())
        }
        FidelityMode::Strict => Err(format!(
            "faithful depiction is not available: {}",
            messages
                .iter()
                .map(|(_, message)| *message)
                .collect::<Vec<_>>()
                .join("; ")
        )),
    }
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
    if input.content.contains("M  V30 BEGIN RGROUP ")
        || input.content.trim_start().starts_with("$MDL")
    {
        return Ok(InputFormat::Rgroup);
    }

    if let Some(extension) = input
        .path
        .as_deref()
        .and_then(Path::extension)
        .and_then(|extension| extension.to_str())
        .map(str::to_ascii_lowercase)
    {
        match extension.as_str() {
            "mol2" => return Ok(InputFormat::Mol2),
            "rxn" => return Ok(InputFormat::Rxn),
            "rsmi" => return Ok(InputFormat::ReactionSmiles),
            "mol" => return Ok(InputFormat::Mol),
            "sdf" => return Ok(InputFormat::Sdf),
            "smi" | "smiles" => return Ok(InputFormat::Smiles),
            _ => {}
        }
    }

    if input.content.trim_start().starts_with("$RXN") {
        return Ok(InputFormat::Rxn);
    }
    if input.content.contains("@<TRIPOS>MOLECULE") {
        return Ok(InputFormat::Mol2);
    }
    if input.content.trim().split('>').count() == 3 && input.content.lines().count() <= 1 {
        return Ok(InputFormat::ReactionSmiles);
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

    Err("could not detect input format; pass --format mol, sdf, smiles, mol2, rxn, reaction-smiles, or rgroup".to_string())
}

fn require_first_record(record: NonZeroUsize, format: &str) -> Result<(), String> {
    if record.get() == 1 {
        Ok(())
    } else {
        Err(format!(
            "--record can only exceed 1 for SDF or MOL2 input, not {format}"
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

impl std::fmt::Display for InputFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{} input",
            match self {
                Self::Smiles => "SMILES",
                Self::Sdf => "SDF",
                Self::Mol => "Molfile",
                Self::Mol2 => "MOL2",
                Self::Rxn => "RXN",
                Self::ReactionSmiles => "Reaction SMILES",
                Self::Rgroup => "RGfile",
                Self::Auto => "automatic",
            }
        )
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
    fn validates_typst_lengths() {
        assert_eq!(parse_typst_length("2.5mm").unwrap(), "2.5mm");
        assert!(parse_typst_length("calc(2mm)").is_err());
        assert!(parse_typst_length("-1em").is_err());
    }
}
