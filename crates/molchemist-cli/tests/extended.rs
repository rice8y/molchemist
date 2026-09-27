use std::path::PathBuf;
use std::process::Command;

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/extended")
        .join(name)
}
fn cli() -> Command {
    Command::new(env!("CARGO_BIN_EXE_molchemist"))
}

#[test]
fn mol2_inspection_preserves_source_and_partial_charge() {
    let output = cli()
        .arg("inspect")
        .arg(fixture("water.mol2"))
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let json: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(json["atoms"][0]["sourceId"], 10);
    assert_eq!(json["atoms"][0]["formalCharge"], 0);
}

#[test]
fn rxn_and_reaction_smiles_report_bond_edits_and_ambiguity() {
    let output = cli()
        .arg("inspect")
        .arg(fixture("oxidation.rxn"))
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let json: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(json["changes"].as_array().unwrap().len(), 1);
    let output = cli()
        .args([
            "dump",
            "--text",
            "CC>>CC",
            "--infer-mapping",
            "--fidelity",
            "strict",
        ])
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("ambiguous"));
    assert!(output.stdout.is_empty());
}

#[test]
fn rgroup_members_keep_attachment_context() {
    let output = cli()
        .arg("inspect")
        .arg(fixture("alternatives.mol"))
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let json: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(json["members"][0]["attachments"][0][1][0], 1);
    assert_eq!(
        json["rawSource"],
        std::fs::read_to_string(fixture("alternatives.mol")).unwrap()
    );
    let output = cli()
        .arg("dump")
        .arg(fixture("alternatives.mol"))
        .args(["--fidelity", "strict"])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let invalid = std::fs::read_to_string(fixture("alternatives.mol"))
        .unwrap()
        .replacen("RGROUPS=(1 1)", "RGROUPS=(1 1) ATTCHPT=1", 1);
    let output = cli()
        .args([
            "dump",
            "--format",
            "rgroup",
            "--text",
            &invalid,
            "--fidelity",
            "strict",
        ])
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
}

#[test]
#[ignore = "requires Typst and Alchemist 0.2.0"]
fn coordinate_layout_contract_and_source_parity() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    let output = Command::new("typst")
        .current_dir(&root)
        .args([
            "query",
            "--root",
            ".",
            "crates/molchemist-cli/tests/fixtures/coordinate-contract.typ",
            "<coordinate-parity>",
            "--field",
            "value",
            "--one",
        ])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let package: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    let mut generator = molchemist_cli::Generator::new().unwrap();
    let code = generator
        .coordinate_to_code(
            "[13CH3:7]C(=O)O",
            "smiles",
            molchemist_cli::RenderMode::Skeletal,
            1,
            "3em",
            &serde_json::json!({"layout": "avoid"}),
        )
        .unwrap();
    assert_eq!(package["smiles"].as_str().unwrap(), code);
}

#[test]
#[ignore = "requires Typst and Alchemist 0.2.0"]
fn extended_package_and_standalone_documents_render() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    let target = std::env::temp_dir().join(format!("molchemist-extended-{}", std::process::id()));
    std::fs::create_dir_all(&target).unwrap();
    let output = Command::new("typst")
        .current_dir(&root)
        .args([
            "compile",
            "--root",
            ".",
            "crates/molchemist-cli/tests/fixtures/extended-rendering.typ",
        ])
        .arg(target.join("package.pdf"))
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let cases = [
        vec![
            "--smiles".into(),
            "[13CH3:7]C(=O)O".into(),
            "--layout".into(),
            "avoid".into(),
        ],
        vec![
            fixture("water.mol2").to_string_lossy().into_owned(),
            "--layout".into(),
            "coordinates".into(),
        ],
        vec![
            fixture("tetra-3d.mol").to_string_lossy().into_owned(),
            "--layout".into(),
            "reflow".into(),
            "--infer-stereo".into(),
        ],
        vec![fixture("oxidation.rxn").to_string_lossy().into_owned()],
        vec![
            "--text".into(),
            "[CH3:1][CH2:2][CH2:3]Br>>[CH3:1][CH2:2][CH2:3]O".into(),
        ],
        vec![fixture("alternatives.mol").to_string_lossy().into_owned()],
    ];
    for (i, args) in cases.iter().enumerate() {
        let output = cli()
            .arg("dump")
            .args(args)
            .arg("--standalone")
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "case {i}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        let input = target.join(format!("{i}.typ"));
        std::fs::write(&input, &output.stdout).unwrap();
        let rendered = Command::new("typst")
            .arg("compile")
            .arg(&input)
            .arg(target.join(format!("{i}.svg")))
            .output()
            .unwrap();
        assert!(
            rendered.status.success(),
            "case {i}: {}",
            String::from_utf8_lossy(&rendered.stderr)
        );
    }
}
