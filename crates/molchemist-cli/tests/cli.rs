use std::fs;
use std::io::Write;
use std::path::PathBuf;
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
fn dumps_disconnected_smiles_as_separate_components() {
    let output = molchemist()
        .args(["dump", "--smiles", "[Na+].[Cl-]", "--mode", "abbreviate"])
        .output()
        .unwrap();

    assert!(output.status.success());
    assert!(output.stderr.is_empty());
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.contains("fragment(\"Na^+\", name: \"a0\")"));
    assert!(stdout.contains("operator(none, margin: base-sep * 0.5)"));
    assert!(stdout.contains("fragment(\"Cl^-\", name: \"a1\")"));
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
fn dumps_extended_sdf_bond_semantics_without_collapsing_to_single() {
    let fixture = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/bond-semantics.sdf");
    let output = molchemist()
        .args(["dump", fixture.to_str().unwrap()])
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "bond semantics dump failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output.stderr.is_empty());
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.contains("#import \"@preview/cetz:0.5.2\""));
    for function in [
        "_molchemist-aromatic(",
        "_molchemist-single-or-double(",
        "_molchemist-single-or-aromatic(",
        "_molchemist-double-or-aromatic(",
        "_molchemist-wavy(",
        "_molchemist-coordination-right(",
        "_molchemist-hydrogen(",
    ] {
        assert!(stdout.contains(function), "missing {function}");
    }
}

#[test]
fn dumps_stereochemistry_without_folding_or_dropping_semantics() {
    let fixture = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/stereochemistry.sdf");
    let output = molchemist()
        .args(["dump", fixture.to_str().unwrap(), "--mode", "skeletal"])
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stereochemistry dump failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output.stderr.is_empty());
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.contains("_molchemist-crossed-double("));
    assert!(stdout.contains("cram-filled-left("));
    assert!(stdout.contains("fragment(\"H\", name: \"a1\")"));
    assert!(stdout.contains("Stereo annotations: C CFG=1; OR7 (a0)"));
}

#[test]
fn dumps_l_and_d_alanine_with_the_expected_absolute_bond_styles() {
    let cases = [
        ("skeletal", "cram-filled-right(", "cram-dashed-right("),
        ("abbreviate", "cram-filled-right(", "cram-dashed-right("),
        ("full", "cram-filled-left(", "cram-dashed-left("),
    ];

    for (mode, l_style, d_style) in cases {
        let l_alanine = molchemist()
            .args(["dump", "--smiles", "N[C@@H](C)C(=O)O", "--mode", mode])
            .output()
            .unwrap();
        let d_alanine = molchemist()
            .args(["dump", "--smiles", "N[C@H](C)C(=O)O", "--mode", mode])
            .output()
            .unwrap();
        assert!(l_alanine.status.success(), "L-alanine failed in {mode}");
        assert!(d_alanine.status.success(), "D-alanine failed in {mode}");

        let l_source = String::from_utf8(l_alanine.stdout).unwrap();
        let d_source = String::from_utf8(d_alanine.stdout).unwrap();
        assert!(l_source.contains(l_style), "L-alanine {mode}: {l_source}");
        assert!(d_source.contains(d_style), "D-alanine {mode}: {d_source}");
        assert_ne!(l_source, d_source, "enantiomers collapsed in {mode}");
    }
}

#[test]
fn implicit_and_explicit_hydrogen_alanine_dump_the_same_center_orientation() {
    let sources = ["N[C@@H](C)C(=O)O", "N[C@@]([H])(C)C(=O)O"].map(|smiles| {
        let output = molchemist()
            .args(["dump", "--smiles", smiles, "--mode", "skeletal"])
            .output()
            .unwrap();
        assert!(output.status.success(), "{smiles}");
        String::from_utf8(output.stdout).unwrap()
    });

    assert!(sources
        .iter()
        .all(|source| source.contains("cram-filled-right(")));
}

#[test]
fn relayouts_sdf_records_with_collapsed_coordinates() {
    let fixture = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/layout-robustness.sdf");
    let output = molchemist()
        .args(["dump", fixture.to_str().unwrap(), "--mode", "skeletal"])
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "collapsed-coordinate dump failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.contains("cram-filled-"));
    assert!(!stdout.contains("atom-sep: base-sep * 0, name:"));
    assert!(!stdout.contains("NaN"));
    assert!(!stdout.contains("inf"));
}

#[test]
fn dumps_extended_chirality_as_native_geometry() {
    let cases = [
        ("[Pt@SP2](F)(Cl)(Br)I", false),
        ("[As@TB5](F)(Cl)(Br)(N)S", true),
        ("[Co@OH5](F)(Cl)(Br)(I)(N)S", true),
        ("NC(Br)=[C@AL1]=C(O)C", true),
    ];

    for (smiles, expects_stereo_bond) in cases {
        let output = molchemist()
            .args(["dump", "--smiles", smiles, "--mode", "skeletal"])
            .output()
            .unwrap();

        assert!(
            output.status.success(),
            "extended chirality dump failed for {smiles}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        let stdout = String::from_utf8(output.stdout).unwrap();
        assert!(!stdout.contains("Stereo annotations:"), "{smiles}");
        assert_eq!(stdout.contains("cram-"), expects_stereo_bond, "{smiles}");
    }
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

    let missing = molchemist()
        .args(["dump", path.to_str().unwrap(), "--record", "3"])
        .output()
        .unwrap();
    assert!(!missing.status.success());
    assert!(String::from_utf8(missing.stderr)
        .unwrap()
        .contains("SDF record 3 does not exist; input contains 2 record(s)"));
    fs::remove_file(path).unwrap();
}

#[test]
fn auto_detects_and_dumps_v3000_input() {
    let v3000 = concat!(
        "charged\n",
        "  molchemist\n",
        "\n",
        "  0  0  0     0  0            999 V3000\n",
        "M  V30 BEGIN CTAB\n",
        "M  V30 COUNTS 2 1 0 0 0\n",
        "M  V30 BEGIN ATOM\n",
        "M  V30 1 N 0.0000 0.0000 0.0000 0 CHG=1\n",
        "M  V30 2 O 1.5000 0.0000 0.0000 0 CHG=-1\n",
        "M  V30 END ATOM\n",
        "M  V30 BEGIN BOND\n",
        "M  V30 1 1 1 2\n",
        "M  V30 END BOND\n",
        "M  V30 END CTAB\n",
        "M  END\n",
    );
    let output = molchemist()
        .args(["dump", "--text", v3000, "--mode", "abbreviate"])
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "V3000 input failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.contains("fragment(\"N^+\""));
    assert!(stdout.contains("fragment(\"O^-\""));
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
    let fixture =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/stereochemistry.sdf");
    let generated = molchemist()
        .args([
            "dump",
            fixture.to_str().unwrap(),
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

#[test]
#[ignore = "requires Typst and the alchemist package"]
fn local_typst_package_renders_atom_metadata() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    let fixture = root.join("crates/molchemist-cli/tests/fixtures/atom-metadata.typ");
    let pdf = temp_path("atom-metadata").with_extension("pdf");

    let compiled = Command::new("typst")
        .current_dir(&root)
        .args([
            "compile",
            fixture.to_str().unwrap(),
            pdf.to_str().unwrap(),
            "--root",
            root.to_str().unwrap(),
        ])
        .output()
        .unwrap();

    assert!(
        compiled.status.success(),
        "Typst failed: {}",
        String::from_utf8_lossy(&compiled.stderr)
    );
    assert!(pdf.exists());
    fs::remove_file(pdf).unwrap();
}

#[test]
#[ignore = "requires Typst and the alchemist package"]
fn local_typst_package_selects_sdf_records_and_renders_v3000() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    let fixture = root.join("crates/molchemist-cli/tests/fixtures/sdf-input-fidelity.typ");
    let pdf = temp_path("sdf-input-fidelity").with_extension("pdf");

    let compiled = Command::new("typst")
        .current_dir(&root)
        .args([
            "compile",
            fixture.to_str().unwrap(),
            pdf.to_str().unwrap(),
            "--root",
            root.to_str().unwrap(),
        ])
        .output()
        .unwrap();

    assert!(
        compiled.status.success(),
        "Typst failed: {}",
        String::from_utf8_lossy(&compiled.stderr)
    );
    assert!(pdf.exists());
    fs::remove_file(pdf).unwrap();
}

#[test]
#[ignore = "requires Typst and the alchemist package"]
fn local_typst_package_renders_and_annotates_multiple_components() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    let fixture = root.join("crates/molchemist-cli/tests/fixtures/multicomponent-rendering.typ");
    let pdf = temp_path("multicomponent-rendering").with_extension("pdf");

    let compiled = Command::new("typst")
        .current_dir(&root)
        .args([
            "compile",
            fixture.to_str().unwrap(),
            pdf.to_str().unwrap(),
            "--root",
            root.to_str().unwrap(),
        ])
        .output()
        .unwrap();

    assert!(
        compiled.status.success(),
        "Typst failed: {}",
        String::from_utf8_lossy(&compiled.stderr)
    );
    assert!(pdf.exists());
    fs::remove_file(pdf).unwrap();
}

#[test]
#[ignore = "requires Typst and the alchemist package"]
fn local_typst_package_renders_extended_bond_semantics() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    let fixture = root.join("crates/molchemist-cli/tests/fixtures/bond-semantics-fidelity.typ");
    let pdf = temp_path("bond-semantics-fidelity").with_extension("pdf");

    let compiled = Command::new("typst")
        .current_dir(&root)
        .args([
            "compile",
            fixture.to_str().unwrap(),
            pdf.to_str().unwrap(),
            "--root",
            root.to_str().unwrap(),
        ])
        .output()
        .unwrap();

    assert!(
        compiled.status.success(),
        "Typst failed: {}",
        String::from_utf8_lossy(&compiled.stderr)
    );
    assert!(pdf.exists());
    fs::remove_file(pdf).unwrap();
}

#[test]
#[ignore = "requires Typst and the alchemist package"]
fn local_typst_package_renders_stereochemistry_fidelity_cases() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    let fixture = root.join("crates/molchemist-cli/tests/fixtures/stereochemistry-fidelity.typ");
    let pdf = temp_path("stereochemistry-fidelity").with_extension("pdf");

    let compiled = Command::new("typst")
        .current_dir(&root)
        .args([
            "compile",
            fixture.to_str().unwrap(),
            pdf.to_str().unwrap(),
            "--root",
            root.to_str().unwrap(),
        ])
        .output()
        .unwrap();

    assert!(
        compiled.status.success(),
        "Typst failed: {}",
        String::from_utf8_lossy(&compiled.stderr)
    );
    assert!(pdf.exists());
    fs::remove_file(pdf).unwrap();
}

#[test]
#[ignore = "requires Typst and the alchemist package"]
fn local_typst_package_recovers_collapsed_sdf_layouts() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    let fixture = root.join("crates/molchemist-cli/tests/fixtures/layout-robustness.typ");
    let pdf = temp_path("layout-robustness").with_extension("pdf");

    let compiled = Command::new("typst")
        .current_dir(&root)
        .args([
            "compile",
            fixture.to_str().unwrap(),
            pdf.to_str().unwrap(),
            "--root",
            root.to_str().unwrap(),
        ])
        .output()
        .unwrap();

    assert!(
        compiled.status.success(),
        "Typst failed: {}",
        String::from_utf8_lossy(&compiled.stderr)
    );
    assert!(pdf.exists());
    fs::remove_file(pdf).unwrap();
}
