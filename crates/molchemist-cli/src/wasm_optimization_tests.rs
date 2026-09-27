use super::*;
use std::path::PathBuf;

fn compare<T: PartialEq>(
    before: &mut Generator,
    after: &mut Generator,
    label: &str,
    call: impl Fn(&mut Generator) -> Result<T, String>,
) {
    let expected = call(before).unwrap_or_else(|e| panic!("Original {label}: {e}"));
    let actual = call(after).unwrap_or_else(|e| panic!("Optimized {label}: {e}"));
    assert!(expected == actual, "Optimization changed {label}");
}

fn interface(bytes: &[u8]) -> Vec<String> {
    let mut config = Config::default();
    config.wasm_relaxed_simd(false);
    let engine = Engine::new(&config);
    let module = Module::new(&engine, bytes).unwrap();
    let mut entries = module
        .imports()
        .map(|i| format!("import {}::{} {:?}", i.module(), i.name(), i.ty()))
        .chain(
            module
                .exports()
                .map(|e| format!("export {} {:?}", e.name(), e.ty())),
        )
        .collect::<Vec<_>>();
    entries.sort();
    entries
}

#[test]
#[ignore = "requires candidate plugins from just optimize-wasm"]
fn optimized_wasm_preserves_behavior() {
    let directory = PathBuf::from(
        std::env::var_os("MOLCHEMIST_OPTIMIZED_WASM_DIR")
            .expect("MOLCHEMIST_OPTIMIZED_WASM_DIR must contain both candidate plugins"),
    );
    let core = std::fs::read(directory.join("molchemist_plugin.wasm")).unwrap();
    let layout = std::fs::read(directory.join("molchemist_smiles_plugin.wasm")).unwrap();
    assert_eq!(interface(CORE_WASM), interface(&core), "Core interface");
    assert_eq!(
        interface(LAYOUT_WASM),
        interface(&layout),
        "Layout interface"
    );
    let mut before = Generator::new().unwrap();
    let mut after = Generator {
        core: Plugin::new(&core).unwrap(),
        layout: Some(Plugin::new(&layout).unwrap()),
    };
    let modes = [
        RenderMode::Full,
        RenderMode::Abbreviate,
        RenderMode::Skeletal,
    ];
    for smiles in [
        "c1ccccc1",
        "[13CH3:7]C(=O)O",
        "N[C@@H](C)C(=O)O",
        r"F/C=C\F",
        "[CH2-]C",
        "*",
        "[Cr]$[Cr]",
        "[Na+].[Cl-]",
        "[Pt@SP1]1(F)(Cl)CCC1",
    ] {
        for mode in modes {
            compare(
                &mut before,
                &mut after,
                &format!("{smiles} {mode:?}"),
                |g| g.smiles_to_code(smiles, mode, "3em"),
            );
        }
    }
    let fixtures = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let read = |name: &str| std::fs::read_to_string(fixtures.join(name)).unwrap();
    for name in [
        "Structure2D_COMPOUND_CID_241.sdf",
        "bond-semantics.sdf",
        "stereochemistry.sdf",
        "layout-robustness.sdf",
        "ctfile-fidelity.sdf",
    ] {
        let data = read(name);
        compare(&mut before, &mut after, name, |g| {
            g.inspect_sdf_record_json(&data, 1, false)
        });
        for mode in modes {
            compare(&mut before, &mut after, &format!("{name} {mode:?}"), |g| {
                g.sdf_to_code(&data, mode, "3em")
            });
        }
    }
    for config in [
        serde_json::json!({"layout": "reflow"}),
        serde_json::json!({"layout": "reflow", "infer-stereo": true}),
    ] {
        compare(&mut before, &mut after, "3D coordinate reflow", |g| {
            g.coordinate_to_code(
                &read("extended/tetra-3d.mol"),
                "sdf",
                RenderMode::Skeletal,
                1,
                "3em",
                &config,
            )
        });
    }
    let config = serde_json::json!({});
    compare(&mut before, &mut after, "MOL2 normalization", |g| {
        g.mol2_to_sdf(&read("extended/water.mol2"), 1)
    });
    compare(&mut before, &mut after, "R-group inspection", |g| {
        g.inspect_rgroup(&read("extended/alternatives.mol"))
    });
    compare(&mut before, &mut after, "R-group rendering", |g| {
        g.rgroup_to_code(
            &read("extended/alternatives.mol"),
            RenderMode::Skeletal,
            "3em",
            &config,
        )
    });
    compare(&mut before, &mut after, "SUP expansion", |g| {
        g.expand_superatoms(&read("rdkit/Sgroups_Abbreviations.mol"), 1)
    });
    for (data, format, infer) in [
        (
            "[CH3:1][CH2:2][CH2:3]Br>>[CH3:1][CH2:2][CH2:3]O".into(),
            "reaction-smiles",
            false,
        ),
        ("CCO>>CC=O".into(), "reaction-smiles", true),
        ("CC>>CC".into(), "reaction-smiles", true),
        (read("extended/oxidation.rxn"), "rxn", false),
    ] {
        compare(&mut before, &mut after, "Reaction inspection", |g| {
            g.inspect_reaction(&data, format, infer, 200000)
        });
        let analysis = before
            .inspect_reaction(&data, format, infer, 200000)
            .unwrap();
        compare(&mut before, &mut after, "Reaction rendering", |g| {
            g.reaction_to_code(&analysis, RenderMode::Skeletal, "3em", &config, "reaction")
        });
    }
    for invalid in ["C(", "[NH4", "C1CC"] {
        compare(&mut before, &mut after, "Invalid SMILES diagnostics", |g| {
            g.smiles_to_code(invalid, RenderMode::Skeletal, "3em")
                .err()
                .ok_or_else(|| "Invalid SMILES was accepted".into())
        });
    }
    compare(&mut before, &mut after, "Invalid SDF diagnostics", |g| {
        g.inspect_sdf_record_json("invalid", 1, false)
            .err()
            .ok_or_else(|| "Invalid SDF was accepted".into())
    });
}
