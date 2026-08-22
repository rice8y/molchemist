use std::fs;
use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Stdio};
use std::time::{SystemTime, UNIX_EPOCH};

const CID_241: &str = include_str!("fixtures/Structure2D_COMPOUND_CID_241.sdf");
const CID_93406: &str = include_str!("fixtures/Structure2D_COMPOUND_CID_93406.sdf");
const CTFILE_V3000: &str = include_str!("fixtures/ctfile-fidelity.sdf");
const FIDELITY_V3000: &str = concat!(
    "foundation\n",
    "  molchemist\n",
    "semantic record\n",
    "  0  0  0     0  0            999 V3000\n",
    "M  V30 BEGIN CTAB\n",
    "M  V30 COUNTS 2 1 1 0 0\n",
    "M  V30 BEGIN ATOM\n",
    "M  V30 10 C 0.0000 0.0000 0.0000 7 MASS=13\n",
    "M  V30 20 O 1.5000 0.0000 0.0000 0\n",
    "M  V30 END ATOM\n",
    "M  V30 BEGIN BOND\n",
    "M  V30 7 1 10 20 RXCTR=1\n",
    "M  V30 END BOND\n",
    "M  V30 BEGIN SGROUP\n",
    "M  V30 3 SUP ATOMS=(1 10) LABEL=\"Me\"\n",
    "M  V30 END SGROUP\n",
    "M  V30 END CTAB\n",
    "M  END\n",
    "> <NAME>\nfirst\nline two\n\n",
    "> <NAME>\nsecond\n\n",
    "$$$$\n",
);

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

fn fixture_path(relative: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(relative)
}

fn inspect_fixture(relative: &str, record: Option<usize>) -> serde_json::Value {
    let path = fixture_path(relative);
    let mut command = molchemist();
    command.arg("inspect").arg(&path);
    if let Some(record) = record {
        command.arg("--record").arg(record.to_string());
    }
    let output = command.output().unwrap();
    assert!(
        output.status.success(),
        "inspection failed for {relative}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        output.stderr.is_empty(),
        "unexpected warning for {relative}"
    );
    serde_json::from_slice(&output.stdout).unwrap()
}

fn svg_path_tag_with_fill<'a>(svg: &'a str, fill: &str) -> &'a str {
    let fill_attribute = format!("fill=\"#{fill}");
    for (tag_start, _) in svg.match_indices("<path") {
        let Some(relative_end) = svg[tag_start..].find("/>") else {
            continue;
        };
        let tag_end = tag_start + relative_end;
        let tag = &svg[tag_start..tag_end];
        if tag.contains(&fill_attribute) {
            return tag;
        }
    }
    panic!("SVG path with fill #{fill} was not found")
}

fn svg_path_data_with_fill<'a>(svg: &'a str, fill: &str) -> &'a str {
    let tag = svg_path_tag_with_fill(svg, fill);
    let data_start = tag
        .find(" d=\"")
        .unwrap_or_else(|| panic!("SVG path with fill #{fill} has no path data"))
        + 4;
    let data_end = data_start
        + tag[data_start..]
            .find('"')
            .unwrap_or_else(|| panic!("SVG path with fill #{fill} has invalid path data"));
    &tag[data_start..data_end]
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
fn dumps_atom_classes_as_inline_label_suffixes() {
    for mode in ["full", "abbreviate", "skeletal"] {
        let output = molchemist()
            .args(["dump", "--smiles", "[CH3:1]O", "--mode", mode])
            .output()
            .unwrap();

        assert!(output.status.success(), "atom class failed in {mode}");
        assert!(output.stderr.is_empty());
        let stdout = String::from_utf8(output.stdout).unwrap();
        assert!(stdout.contains(" + [:1]))"), "atom class missing in {mode}");
        assert!(
            !stdout.contains("br: [:1]"),
            "atom class is subscripted in {mode}"
        );
        if mode != "full" {
            assert!(
                stdout.contains(
                    "math.attach(math.attach([C#math.attach([H], b: [3], t: std.hide([3]))]) + [:1])"
                ),
                "hydrogen count or atom class order changed in {mode}"
            );
        }
    }
}

#[test]
fn dumps_smiles_quadruple_bonds_in_every_mode() {
    for mode in ["full", "abbreviate", "skeletal"] {
        let output = molchemist()
            .args(["dump", "--smiles", "[Cr]$[Cr]", "--mode", mode])
            .output()
            .unwrap();

        assert!(
            output.status.success(),
            "quadruple bond failed in {mode}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(output.stderr.is_empty());
        let stdout = String::from_utf8(output.stdout).unwrap();
        assert!(stdout.contains("#let _molchemist-quadruple = build-link"));
        assert!(stdout.contains("_molchemist-quadruple(absolute:"));
        assert!(stdout.contains("fragment(\"Cr\", name: \"a0\")"));
        assert!(stdout.contains("fragment(\"Cr\", name: \"a1\")"));
    }
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
    assert!(stdout.contains("math.attach([Na], tr: [+], br: std.hide([+]))"));
    assert!(stdout.contains("operator(none, margin: base-sep * 0.5)"));
    assert!(stdout.contains("math.attach([Cl], tr: [−], br: std.hide([−]))"));
}

#[test]
fn dumps_component_labels_with_balanced_attachments() {
    let output = molchemist()
        .args(["dump", "--smiles", "[H+].C.[Cl-]", "--mode", "skeletal"])
        .output()
        .unwrap();

    assert!(output.status.success());
    assert!(output.stderr.is_empty());
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.contains("math.attach([H], tr: [+], br: std.hide([+]))"));
    assert!(stdout.contains("[C#math.attach([H], b: [4], t: std.hide([4]))]"));
    assert!(stdout.contains("math.attach([Cl], tr: [−], br: std.hide([−]))"));
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
    assert!(stdout.contains("mark: (end: \">\", fill: black)"));
    assert!(stdout.contains("mark: (start: \">\", fill: black)"));
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
    assert!(stdout.contains("math.attach([N], tr: [+], br: std.hide([+]))"));
    assert!(stdout.contains("math.attach([O], tr: [−], br: std.hide([−]))"));
}

#[test]
fn inspects_semantic_records_as_ordered_json() {
    let output = molchemist()
        .args(["inspect", "--text", FIDELITY_V3000])
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "inspection failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output.stderr.is_empty());
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(value["schemaVersion"], 1);
    assert_eq!(value["format"], "molfile-v3000");
    assert_eq!(value["atoms"][0]["index"], 0);
    assert_eq!(value["atoms"][0]["sourceId"], 10);
    assert_eq!(value["atoms"][0]["atomMap"], 7);
    assert_eq!(value["bonds"][0]["sourceId"], 7);
    assert_eq!(value["sgroups"][0]["kind"], "superatom");
    assert_eq!(value["properties"][0]["name"], "NAME");
    assert_eq!(value["properties"][0]["value"], "first\nline two");
    assert_eq!(value["properties"][1]["name"], "NAME");
    assert_eq!(value["properties"][1]["value"], "second");
    assert_eq!(
        value["diagnostics"][0]["code"],
        "reaction-center-not-depicted"
    );
    assert!(value["rawRecord"]
        .as_str()
        .unwrap()
        .starts_with("foundation\n"));
}

#[test]
fn fidelity_policy_warns_rejects_or_ignores_unsupported_depiction() {
    let warning = molchemist()
        .args(["dump", "--text", FIDELITY_V3000])
        .output()
        .unwrap();
    assert!(warning.status.success());
    let stderr = String::from_utf8(warning.stderr).unwrap();
    assert!(stderr.contains("warning[reaction-center-not-depicted]"));

    let strict = molchemist()
        .args(["dump", "--text", FIDELITY_V3000, "--fidelity", "strict"])
        .output()
        .unwrap();
    assert!(!strict.status.success());
    assert!(String::from_utf8(strict.stderr)
        .unwrap()
        .contains("faithful depiction is not available"));

    let ignored = molchemist()
        .args(["dump", "--text", FIDELITY_V3000, "--fidelity", "ignore"])
        .output()
        .unwrap();
    assert!(ignored.status.success());
    assert!(ignored.stderr.is_empty());
}

#[test]
fn dumps_ctfile_fidelity_overlays_without_warnings() {
    let output = molchemist()
        .args(["dump", "--text", CTFILE_V3000, "--fidelity", "strict"])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "CTfile dump failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output.stderr.is_empty());
    let source = String::from_utf8(output.stdout).unwrap();
    assert!(source.contains("draw-skeleton(name: \"molchemist-structure\""));
    assert!(source.contains("cetz.draw.on-layer(−1"));
    assert!(source.contains("_molchemist-highlight-path"));
    assert!(source.contains("molchemist-sgroup-3-right"));
    assert!(source.contains("!\\[C,N\\]"));
    assert!(source.contains("(H1) (s3) (u) (r2)"));
    assert!(source.contains("rel: (0pt, base-sep * −0.36)"));
    assert!(source.contains("text(size: 0.44em, fill: luma(38%))"));
    assert!(!source.contains("implicit H >= 1"));
    assert!(source.contains("[rn]"));
    assert!(source.contains("text(size: 0.56em, fill: luma(32%))"));
}

#[test]
fn real_pubchem_fixtures_preserve_sdf_metadata() {
    let benzene = inspect_fixture("Structure2D_COMPOUND_CID_241.sdf", None);
    assert_eq!(benzene["name"], "241");
    assert_eq!(benzene["atoms"].as_array().unwrap().len(), 12);
    assert_eq!(benzene["bonds"].as_array().unwrap().len(), 12);
    let benzene_properties = benzene["properties"]
        .as_array()
        .unwrap()
        .iter()
        .map(|property| property["name"].as_str().unwrap())
        .collect::<Vec<_>>();
    assert_eq!(
        benzene_properties,
        [
            "PUBCHEM_COMPOUND_CID",
            "PUBCHEM_IUPAC_NAME",
            "PUBCHEM_MOLECULAR_FORMULA",
            "PUBCHEM_SMILES",
            "PUBCHEM_COORDINATE_TYPE",
        ]
    );

    let complex = inspect_fixture("Structure2D_COMPOUND_CID_93406.sdf", None);
    assert_eq!(complex["name"], "93406");
    assert_eq!(complex["atoms"].as_array().unwrap().len(), 29);
    assert_eq!(complex["bonds"].as_array().unwrap().len(), 30);
    let complex_properties = complex["properties"].as_array().unwrap();
    assert_eq!(complex_properties.len(), 34);
    assert_eq!(complex_properties[0]["name"], "PUBCHEM_COMPOUND_CID");
    assert_eq!(
        complex_properties.last().unwrap()["name"],
        "PUBCHEM_BONDANNOTATIONS"
    );
    assert_eq!(complex["diagnostics"], serde_json::json!([]));
}

#[test]
fn real_rdkit_query_fixtures_preserve_atom_lists_and_constraints() {
    let expected_lists: [(usize, &[&str]); 3] =
        [(1, &["C", "N", "O"]), (2, &["N", "C", "O"]), (3, &["C"])];
    for (record, expected) in expected_lists {
        let value = inspect_fixture("rdkit/github8823.sdf", Some(record));
        let query = &value["atoms"][0]["query"];
        let elements = query["elements"].as_array().unwrap();
        assert_eq!(
            elements
                .iter()
                .map(|element| element.as_str().unwrap())
                .collect::<Vec<_>>(),
            expected,
            "atom-list order changed in record {record}"
        );
        assert_eq!(query["isNotList"], false);
        assert_eq!(value["diagnostics"], serde_json::json!([]));
    }

    let atom_query = inspect_fixture("rdkit/AtomQuery1.mol", None);
    assert_eq!(atom_query["format"], "molfile-v3000");
    assert_eq!(atom_query["atoms"][6]["sourceId"], 7);
    assert_eq!(atom_query["atoms"][6]["hydrogenCount"], 3);
    assert_eq!(atom_query["atoms"][6]["query"]["hydrogenCount"], 3);
}

#[test]
fn real_rdkit_ctfiles_preserve_bond_topology_and_sgroups() {
    for (fixture, expected_topology, marker) in [
        ("rdkit/RingBondQuery.mol", 1, "[rn]"),
        ("rdkit/ChainBondQuery.mol", 2, "[ch]"),
    ] {
        let value = inspect_fixture(fixture, None);
        assert_eq!(value["bonds"][4]["topology"], expected_topology);

        let output = molchemist()
            .args([
                "dump",
                fixture_path(fixture).to_str().unwrap(),
                "--fidelity",
                "strict",
            ])
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "dump failed for {fixture}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(output.stderr.is_empty());
        assert!(String::from_utf8(output.stdout).unwrap().contains(marker));
    }

    let v2000_sru = inspect_fixture("rdkit/Sgroups_SRU_01.mol", None);
    let sru = &v2000_sru["sgroups"][0];
    assert_eq!(sru["kind"], "structure-repeat-unit");
    assert_eq!(sru["connectivity"], "HT");
    assert_eq!(sru["subscript"], "n");
    assert_eq!(sru["brackets"].as_array().unwrap().len(), 2);
    assert_eq!(sru["crossingBondSourceIds"], serde_json::json!([10, 9]));

    let v3000_sru = inspect_fixture("rdkit/repeat_groups_query1.mol", None);
    let repeat = &v3000_sru["sgroups"][0];
    assert_eq!(repeat["kind"], "structure-repeat-unit");
    assert_eq!(repeat["connectivity"], "HT");
    assert_eq!(repeat["subscript"], "1-3");
    assert_eq!(repeat["brackets"].as_array().unwrap().len(), 2);
    assert_eq!(repeat["crossingBondSourceIds"], serde_json::json!([1, 6]));

    let data = inspect_fixture("rdkit/Sgroups_Data_01.mol", None);
    assert_eq!(data["sgroups"].as_array().unwrap().len(), 2);
    assert_eq!(data["sgroups"][0]["kind"], "data");
    assert_eq!(data["sgroups"][0]["fieldName"], "pH");
    assert_eq!(
        data["sgroups"][0]["fieldValue"].as_str().unwrap().trim(),
        "4.6"
    );
    assert_eq!(data["sgroups"][1]["fieldName"], "Stereo");
    assert_eq!(
        data["sgroups"][1]["fieldValue"].as_str().unwrap().trim(),
        "E/Z unknown"
    );
}

#[test]
fn real_rdkit_highlight_fixture_preserves_separate_atom_and_bond_collections() {
    let value = inspect_fixture("rdkit/v3k.crash1.mol", None);
    assert_eq!(value["format"], "molfile-v3000");
    assert_eq!(value["diagnostics"], serde_json::json!([]));
    assert_eq!(value["collections"].as_array().unwrap().len(), 2);
    assert_eq!(value["collections"][0]["kind"], "highlight");
    assert_eq!(
        value["collections"][0]["atomSourceIds"],
        serde_json::json!([])
    );
    assert_eq!(
        value["collections"][0]["bondSourceIds"],
        serde_json::json!([7])
    );
    assert_eq!(value["collections"][1]["kind"], "highlight");
    assert_eq!(
        value["collections"][1]["atomSourceIds"],
        serde_json::json!([7])
    );
    assert_eq!(
        value["collections"][1]["bondSourceIds"],
        serde_json::json!([])
    );
}

#[test]
fn synthetic_highlight_shape_corpus_keeps_atom_and_bond_membership_distinct() {
    let expected = [
        (vec![2], vec![]),
        (vec![2], vec![]),
        (vec![], vec![1]),
        (vec![1, 2, 3], vec![1, 2]),
        (vec![1, 4], vec![2]),
    ];

    for (record, (atom_ids, bond_ids)) in expected.into_iter().enumerate() {
        let value = inspect_fixture("highlight-shapes.sdf", Some(record + 1));
        assert_eq!(value["diagnostics"], serde_json::json!([]));
        assert_eq!(value["collections"].as_array().unwrap().len(), 1);
        assert_eq!(value["collections"][0]["kind"], "highlight");
        assert_eq!(
            value["collections"][0]["atomSourceIds"],
            serde_json::json!(atom_ids),
            "atom membership changed in highlight record {}",
            record + 1
        );
        assert_eq!(
            value["collections"][0]["bondSourceIds"],
            serde_json::json!(bond_ids),
            "bond membership changed in highlight record {}",
            record + 1
        );
    }
}

#[test]
fn documentation_highlight_fixture_matches_the_regression_corpus() {
    assert_eq!(
        include_str!("fixtures/highlight-shapes.sdf"),
        include_str!("../../../package/docs/assets/highlight-shapes.sdf")
    );
}

#[test]
fn documentation_superatom_fixture_matches_the_real_regression_corpus() {
    assert_eq!(
        include_str!("fixtures/rdkit/Sgroups_Abbreviations.mol"),
        include_str!("../../../package/docs/assets/Sgroups_Abbreviations.mol")
    );
    let documentation = include_str!("../../../package/docs/documentation.typ");
    assert!(documentation.contains("=== Multi-atom Superatom Contraction"));
    assert!(documentation.contains("#render-mol(data, skeletal: true, fidelity: \"strict\")"));
}

#[test]
fn documentation_attachment_fixtures_match_the_regression_corpus() {
    assert_eq!(
        include_str!("fixtures/ctfile-attachments-collections.sdf"),
        include_str!("../../../package/docs/assets/ctfile-attachments-collections.sdf")
    );
    assert_eq!(
        include_str!("fixtures/rdkit/Sgroups_Link_01.mol"),
        include_str!("../../../package/docs/assets/Sgroups_Link_01.mol")
    );
    let documentation = include_str!("../../../package/docs/documentation.typ");
    assert!(documentation.contains("=== Attachment Points, Link Nodes, and Collections"));
    assert!(documentation
        .contains("#let data = read(\"ctfile-attachments-collections.sdf\", encoding: none)"));
    assert!(documentation
        .contains("#let example = example.with(side-by-side: false, breakable: false)"));
    let section_start = documentation
        .find("=== Attachment Points, Link Nodes, and Collections")
        .unwrap();
    let section_end = documentation[section_start..]
        .find("=== Highlight Depiction")
        .map(|offset| section_start + offset)
        .unwrap();
    let section = &documentation[section_start..section_end];
    assert_eq!(section.matches("#example[").count(), 4);
    assert!(section.contains("dotted branches are the `ATTACH=ANY` endpoint alternatives"));
    assert!(section.contains("the link repeat range is 1–3"));
    assert!(section.contains("inspect-mol(data, record: 3)"));
    assert!(section.contains("collection.name"));
}

#[test]
fn documentation_places_fidelity_details_after_the_introduction() {
    let documentation = include_str!("../../../package/docs/documentation.typ");
    let getting_started = documentation.find("= Getting Started").unwrap();
    let choosing_input = documentation.find("== Choosing an Input Format").unwrap();
    let compatibility = documentation
        .find("= Compatibility and Test Coverage")
        .unwrap();
    let fidelity = documentation.find("= Input and Semantic Fidelity").unwrap();
    let inspection = documentation
        .find("== Semantic Inspection and Fidelity Policy")
        .unwrap();
    let ctfile = documentation
        .find("== CTfile Queries, SGroups, Collections, and Attachments")
        .unwrap();
    let bonds = documentation.find("== Bond Semantics").unwrap();

    assert!(getting_started < choosing_input);
    assert!(choosing_input < compatibility);
    assert!(compatibility < fidelity);
    assert!(fidelity < inspection);
    assert!(inspection < ctfile);
    assert!(ctfile < bonds);

    let introduction = &documentation[getting_started..compatibility];
    assert!(!introduction.contains("#inspect-mol("));
    assert!(!introduction.contains("ctfile: ("));
    assert!(!introduction.contains("fidelity: \"strict\""));
}

#[test]
fn real_rdkit_superatom_fixture_contracts_in_strict_mode() {
    let fixture = "rdkit/Sgroups_Abbreviations.mol";
    let value = inspect_fixture(fixture, None);
    assert_eq!(value["sgroups"].as_array().unwrap().len(), 2);
    assert_eq!(value["sgroups"][0]["kind"], "superatom");
    assert_eq!(value["sgroups"][0]["label"], "NO2");
    assert_eq!(value["sgroups"][1]["label"], "COOH");
    assert_eq!(value["diagnostics"], serde_json::json!([]));

    let output = molchemist()
        .args([
            "dump",
            fixture_path(fixture).to_str().unwrap(),
            "--fidelity",
            "strict",
        ])
        .output()
        .unwrap();
    assert!(output.status.success());
    assert!(output.stderr.is_empty());
    let source = String::from_utf8(output.stdout).unwrap();
    assert_eq!(source.matches("fragment(").count(), 8);
    assert!(source.contains("[NO₂]"));
    assert!(source.contains("[COOH]"));
}

#[test]
fn attachment_link_node_and_collection_fixtures_preserve_standard_semantics() {
    let link = inspect_fixture("rdkit/Sgroups_Link_01.mol", None);
    assert_eq!(link["linkNodes"][0]["minRepeat"], 1);
    assert_eq!(link["linkNodes"][0]["maxRepeat"], 3);
    assert_eq!(
        link["linkNodes"][0]["connections"]
            .as_array()
            .unwrap()
            .len(),
        2
    );
    assert_eq!(link["diagnostics"], serde_json::json!([]));

    let sap = inspect_fixture("rdkit/sgroup_ap_bug.mol", None);
    assert_eq!(sap["sgroups"][0]["attachmentPoints"][0]["atomSourceId"], 1);
    assert_eq!(sap["sgroups"][2]["attachmentPoints"][0]["id"], "Al");
    assert_eq!(sap["diagnostics"], serde_json::json!([]));

    let attachments = inspect_fixture("ctfile-attachments-collections.sdf", Some(1));
    assert_eq!(
        attachments["atoms"][0]["attachmentPoints"],
        serde_json::json!([])
    );
    assert_eq!(
        attachments["bonds"][6]["endpointSourceIds"],
        serde_json::json!([20, 30, 40])
    );
    assert_eq!(attachments["bonds"][6]["attachmentMode"], "ANY");
    assert_eq!(attachments["diagnostics"], serde_json::json!([]));

    let atom_attributes = inspect_fixture("ctfile-attachments-collections.sdf", Some(2));
    assert_eq!(
        atom_attributes["atoms"][0]["query"]["elements"],
        serde_json::json!(["C", "N"])
    );
    assert_eq!(atom_attributes["atoms"][0]["query"]["isNotList"], true);
    assert_eq!(
        atom_attributes["atoms"][1]["rgroupLabels"],
        serde_json::json!([7])
    );
    assert_eq!(atom_attributes["collections"], serde_json::json!([]));
    assert_eq!(atom_attributes["diagnostics"], serde_json::json!([]));

    let collections = inspect_fixture("ctfile-attachments-collections.sdf", Some(3));
    assert_eq!(collections["collections"].as_array().unwrap().len(), 1);
    assert_eq!(
        collections["collections"][0]["name"],
        "MOLCHEMIST/SELECTION"
    );
    assert_eq!(collections["collections"][0]["kind"], "user-defined");
    assert_eq!(
        collections["collections"][0]["atomSourceIds"],
        serde_json::json!([10])
    );
    assert_eq!(
        collections["collections"][0]["bondSourceIds"],
        serde_json::json!([7])
    );
    assert_eq!(
        collections["diagnostics"][0]["code"],
        "user-defined-collection-not-depicted"
    );

    for (fixture, record) in [
        ("rdkit/Sgroups_Link_01.mol", None),
        ("ctfile-attachments-collections.sdf", Some("1")),
        ("ctfile-attachments-collections.sdf", Some("2")),
    ] {
        let mut command = molchemist();
        command.args([
            "dump",
            fixture_path(fixture).to_str().unwrap(),
            "--fidelity",
            "strict",
        ]);
        if let Some(record) = record {
            command.args(["--record", record]);
        }
        let output = command.output().unwrap();
        assert!(
            output.status.success(),
            "strict rendering failed for {fixture}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(output.stderr.is_empty());
    }
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
    for smiles in ["C(", ""] {
        let output = molchemist()
            .args(["dump", "--smiles", smiles])
            .output()
            .unwrap();

        assert!(!output.status.success(), "{smiles:?}");
        assert!(output.stdout.is_empty(), "{smiles:?}");
        let stderr = String::from_utf8(output.stderr).unwrap();
        if smiles.is_empty() {
            assert_eq!(stderr, "error: SMILES input is empty\n");
        } else {
            assert!(
                stderr.starts_with("error: failed to convert SMILES input:"),
                "{smiles:?}"
            );
        }
    }
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
fn local_typst_package_exposes_fidelity_inspection() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    let fixture = root.join("crates/molchemist-cli/tests/fixtures/fidelity-foundation.typ");
    let pdf = temp_path("fidelity-foundation").with_extension("pdf");

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
fn local_typst_package_depicts_ctfile_fidelity_features() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    let fixture = root.join("crates/molchemist-cli/tests/fixtures/ctfile-fidelity.typ");
    let pdf = temp_path("ctfile-fidelity").with_extension("pdf");

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
fn local_typst_package_depicts_variable_attachments_link_nodes_and_atom_attributes() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    let fixture =
        root.join("crates/molchemist-cli/tests/fixtures/ctfile-attachments-collections.typ");
    let pdf = temp_path("ctfile-attachments-collections").with_extension("pdf");

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
fn local_typst_package_renders_expanded_query_details_on_request() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    let fixture = root.join("crates/molchemist-cli/tests/fixtures/ctfile-fidelity-details.typ");
    let pdf = temp_path("ctfile-fidelity-details").with_extension("pdf");

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
fn local_typst_package_renders_real_world_ctfile_corpus() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    let fixture = root.join("crates/molchemist-cli/tests/fixtures/ctfile-real-world-corpus.typ");
    let pdf = temp_path("ctfile-real-world-corpus").with_extension("pdf");

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
fn local_typst_package_renders_highlight_shape_corpus_to_svg() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    let fixture = root.join("crates/molchemist-cli/tests/fixtures/ctfile-highlight-corpus.typ");
    let svg = temp_path("ctfile-highlight-corpus").with_extension("svg");

    let compiled = Command::new("typst")
        .current_dir(&root)
        .args([
            "compile",
            fixture.to_str().unwrap(),
            svg.to_str().unwrap(),
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
    let svg_source = fs::read_to_string(&svg).unwrap();

    let cases = [
        ("ff6b6b", 4),
        ("74c0fc", 1),
        ("51cf66", 1),
        ("ffd43b", 3),
        ("cc5de8", 9),
        ("ff922b", 5),
    ];
    for (paint, expected_subpaths) in cases {
        let tag = svg_path_tag_with_fill(&svg_source, paint);
        assert!(tag.contains("fill-rule=\"nonzero\""));
        let path = svg_path_data_with_fill(&svg_source, paint);
        assert_eq!(
            path.matches('Z').count(),
            expected_subpaths,
            "highlight topology changed for #{paint}"
        );
    }

    let skeletal_atom = svg_path_data_with_fill(&svg_source, "74c0fc");
    assert!(!skeletal_atom.contains("h "));
    assert!(!skeletal_atom.contains("v "));
    let query_glyph = svg_path_data_with_fill(&svg_source, "51cf66");
    assert!(query_glyph.contains("h "));
    assert!(query_glyph.contains("v "));

    fs::remove_file(svg).unwrap();
}

#[test]
#[ignore = "requires Typst and the alchemist package"]
fn local_typst_package_strict_fidelity_rejects_unsupported_depiction() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    let fixture = root.join("crates/molchemist-cli/tests/fixtures/fidelity-strict.typ");
    let pdf = temp_path("fidelity-strict").with_extension("pdf");

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

    assert!(!compiled.status.success());
    assert!(
        String::from_utf8_lossy(&compiled.stderr).contains("faithful depiction is not available")
    );
    assert!(!pdf.exists());
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
