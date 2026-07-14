use std::fs;
use std::io::Write;
use std::process::{Command, Stdio};
use std::time::{SystemTime, UNIX_EPOCH};

const CID_241: &str = include_str!("fixtures/Structure2D_COMPOUND_CID_241.sdf");
const CID_93406: &str = include_str!("fixtures/Structure2D_COMPOUND_CID_93406.sdf");

fn molchemist() -> Command {
    Command::new(env!("CARGO_BIN_EXE_molchemist"))
}

fn temp_path(name: &str) -> std::path::PathBuf {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    std::env::temp_dir().join(format!("molchemist-{name}-{nonce}"))
}

#[test]
fn dumps_smiles_to_stdout_without_diagnostics() {
    let output = molchemist()
        .args(["dump", "--smiles", "c1ccccc1", "--mode", "skeletal"])
        .output()
        .unwrap();

    assert!(output.status.success());
    assert!(output.stderr.is_empty());
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.starts_with("#let base-sep = 3em\n#skeletize({\n"));
    assert!(stdout.contains("double(absolute:"));
    assert!(stdout.ends_with("})"));
}

#[test]
fn dumps_sdf_file_selected_by_extension() {
    let fixture = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/Structure2D_COMPOUND_CID_241.sdf");
    let output = molchemist()
        .args(["dump", fixture.to_str().unwrap()])
        .output()
        .unwrap();

    assert!(output.status.success());
    assert!(output.stderr.is_empty());
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.starts_with("#let base-sep = 3em\n#skeletize({\n"));
    assert!(stdout.contains("name: \"a0\""));
}

#[test]
fn reads_direct_text_with_an_explicit_format() {
    let output = molchemist()
        .args([
            "dump",
            "--text",
            "c1ccncc1",
            "--format",
            "smiles",
            "--mode",
            "abbreviate",
        ])
        .output()
        .unwrap();

    assert!(output.status.success());
    assert!(output.stderr.is_empty());
    assert!(String::from_utf8(output.stdout)
        .unwrap()
        .contains("fragment(\"N\""));
}

#[test]
fn selects_a_one_based_record_from_sdf() {
    let path = temp_path("records.sdf");
    fs::write(&path, format!("{CID_241}{CID_93406}")).unwrap();

    let first = molchemist()
        .args([
            "dump",
            path.to_str().unwrap(),
            "--record",
            "1",
            "--mode",
            "skeletal",
        ])
        .output()
        .unwrap();
    let second = molchemist()
        .args([
            "dump",
            path.to_str().unwrap(),
            "--record",
            "2",
            "--mode",
            "skeletal",
        ])
        .output()
        .unwrap();

    assert!(
        first.status.success(),
        "first record failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second record failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_ne!(first.stdout, second.stdout);
    fs::remove_file(path).unwrap();
}

#[test]
fn applies_atom_separation_and_indentation_options() {
    let output = molchemist()
        .args([
            "dump",
            "--smiles",
            "CCO",
            "--mode",
            "abbreviate",
            "--atom-sep",
            "4.5em",
            "--indent",
            "4",
        ])
        .output()
        .unwrap();

    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.starts_with("#let base-sep = 4.5em\n#skeletize({\n    "));
}

#[test]
fn reads_smiles_from_stdin_when_no_source_is_given() {
    let mut child = molchemist()
        .args(["dump", "--format", "smiles", "--mode", "abbreviate"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child.stdin.take().unwrap().write_all(b"CCO\n").unwrap();
    let output = child.wait_with_output().unwrap();

    assert!(output.status.success());
    assert!(output.stderr.is_empty());
    assert!(String::from_utf8(output.stdout)
        .unwrap()
        .contains("fragment(\"OH\""));
}

#[test]
fn writes_a_standalone_document_to_a_file() {
    let path = temp_path("standalone.typ");
    let output = molchemist()
        .args([
            "dump",
            "--smiles",
            "CC(=O)O",
            "--mode",
            "skeletal",
            "--standalone",
            "--output",
            path.to_str().unwrap(),
        ])
        .output()
        .unwrap();

    assert!(output.status.success());
    assert!(output.stdout.is_empty());
    assert!(output.stderr.is_empty());
    let document = fs::read_to_string(&path).unwrap();
    assert!(document.starts_with("#import \"@preview/alchemist:0.2.0\": *\n"));
    assert!(document.contains("#set page(width: auto, height: auto, margin: 3mm)"));
    assert!(document.ends_with("})"));
    fs::remove_file(path).unwrap();
}

#[test]
fn rejects_invalid_smiles_without_polluting_stdout() {
    let output = molchemist()
        .args(["dump", "--smiles", "C("])
        .output()
        .unwrap();

    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
    assert!(String::from_utf8(output.stderr)
        .unwrap()
        .starts_with("error: failed to convert SMILES input:"));
}

#[test]
#[ignore = "requires Typst and the alchemist package"]
fn generated_standalone_document_compiles_with_typst() {
    let source = temp_path("compile.typ");
    let pdf = source.with_extension("pdf");
    let generated = molchemist()
        .args([
            "dump",
            "--smiles",
            "CC(=O)O",
            "--mode",
            "skeletal",
            "--standalone",
            "--output",
            source.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(generated.status.success());

    let compiled = Command::new("typst")
        .args(["compile", source.to_str().unwrap(), pdf.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(
        compiled.status.success(),
        "Typst failed: {}",
        String::from_utf8_lossy(&compiled.stderr)
    );
    assert!(pdf.exists());
    fs::remove_file(source).unwrap();
    fs::remove_file(pdf).unwrap();
}
