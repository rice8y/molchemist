use std::collections::HashMap;
use std::path::PathBuf;
use std::process::Command;

use molchemist_cli::{Generator, RenderMode};

#[test]
#[ignore = "requires Typst and the embedded WASM modules"]
fn typst_and_cli_outputs_match_exactly() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let output = Command::new("typst")
        .current_dir(&root)
        .args([
            "eval",
            "query(<parity>).first().value",
            "--in",
            "tests/fixtures/typst-parity.typ",
            "--root",
            ".",
        ])
        .output()
        .expect("Typst must be installed for the parity test");
    assert!(
        output.status.success(),
        "Typst parity fixture failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let wasm: HashMap<String, String> = serde_json::from_slice(&output.stdout).unwrap();

    let cases = [
        ("benzene", "c1ccccc1", RenderMode::Skeletal),
        (
            "charged",
            "OCCc1c(C)[n+](=cs1)Cc2cnc(C)nc(N)2",
            RenderMode::Abbreviate,
        ),
        ("chiral", "N[C@@H](C)C(=O)O", RenderMode::Full),
        ("ez", r"F/C=C\F", RenderMode::Skeletal),
        (
            "complex",
            "CC[C@@H]([C@@H]1[C@H](C[C@@](O1)(CC)[C@H]2CC[C@@]([C@@H](O2)C)(CC)O)C)C(=O)[C@@H](C)[C@H]([C@H](C)CCC3=C(C=C(C(=C3C(=O)O)O)C)Br)O",
            RenderMode::Abbreviate,
        ),
    ];

    let mut generator = Generator::new().unwrap();
    let sdf = include_str!("fixtures/Structure2D_COMPOUND_CID_241.sdf");
    let cli = generator
        .sdf_to_code(sdf, RenderMode::Abbreviate, "3em", 2)
        .unwrap();
    assert_eq!(wasm.get("sdf").unwrap(), &cli, "SDF parity failed");

    for (name, smiles, mode) in cases {
        let cli = generator.smiles_to_code(smiles, mode, "3em", 2).unwrap();
        assert_eq!(wasm.get(name).unwrap(), &cli, "parity failed for {name}");
    }
}

#[test]
#[ignore = "requires Typst, local WASM modules, and generated PubChem fixtures"]
fn first_hundred_pubchem_compounds_match_exactly() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    let output = Command::new("typst")
        .current_dir(&root)
        .args([
            "eval",
            "query(<pubchem-parity>).first().value",
            "--in",
            "test-results/cli/wasm-parity.typ",
            "--root",
            ".",
        ])
        .output()
        .expect("Typst must be installed for the parity test");
    assert!(
        output.status.success(),
        "Typst PubChem parity fixture failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let wasm: HashMap<String, String> = serde_json::from_slice(&output.stdout).unwrap();

    let source =
        std::fs::read_to_string(root.join("test-results/smiles/pubchem-cid-1-100.json")).unwrap();
    let pubchem: serde_json::Value = serde_json::from_str(&source).unwrap();
    let compounds = pubchem["PropertyTable"]["Properties"].as_array().unwrap();
    let mut generator = Generator::new().unwrap();
    let mut mismatches = Vec::new();
    for compound in compounds {
        let cid = compound["CID"].as_u64().unwrap().to_string();
        let smiles = compound["SMILES"].as_str().unwrap();
        let cli = generator
            .smiles_to_code(smiles, RenderMode::Abbreviate, "3em", 2)
            .unwrap();
        if let Some(wasm) = wasm.get(&cid) {
            if wasm != &cli {
                let (wasm_line, cli_line) = first_different_line(wasm, &cli);
                mismatches.push(format!("CID {cid}: Typst `{wasm_line}` / CLI `{cli_line}`"));
            }
        } else {
            mismatches.push(format!("CID {cid}: missing WASM output"));
        }
    }

    assert!(
        mismatches.is_empty(),
        "WASM/native output differed:\n{}",
        mismatches.join("\n")
    );
}

fn first_different_line<'a>(left: &'a str, right: &'a str) -> (&'a str, &'a str) {
    left.lines()
        .zip(right.lines())
        .find(|(left, right)| left != right)
        .unwrap_or(("<different length>", "<different length>"))
}
