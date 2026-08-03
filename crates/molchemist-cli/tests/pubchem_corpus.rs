use std::collections::HashSet;
use std::fs;
use std::path::PathBuf;
use std::process::{Command, ExitStatus, Stdio};
use std::thread;
use std::time::{Duration, Instant};

use serde_json::Value;

fn molchemist() -> Command {
    Command::new(env!("CARGO_BIN_EXE_molchemist"))
}

fn corpus_root() -> PathBuf {
    std::env::var_os("MOLCHEMIST_PUBCHEM_CORPUS")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("../..")
                .join(".local-tests/pubchem-visual")
        })
}

fn requested_modes() -> Vec<String> {
    std::env::var("MOLCHEMIST_PUBCHEM_MODES")
        .unwrap_or_else(|_| "skeletal".to_string())
        .split(',')
        .map(str::trim)
        .filter(|mode| !mode.is_empty())
        .map(str::to_string)
        .collect()
}

fn requested_cids() -> Option<HashSet<u64>> {
    std::env::var("MOLCHEMIST_PUBCHEM_CIDS").ok().map(|value| {
        value
            .split(',')
            .map(str::trim)
            .filter(|cid| !cid.is_empty())
            .map(|cid| {
                cid.parse()
                    .expect("MOLCHEMIST_PUBCHEM_CIDS must contain CIDs")
            })
            .collect()
    })
}

fn png_dimensions(bytes: &[u8]) -> Option<(u32, u32)> {
    if !bytes.starts_with(b"\x89PNG\r\n\x1a\n") || bytes.get(12..16) != Some(b"IHDR") {
        return None;
    }
    let width = u32::from_be_bytes(bytes.get(16..20)?.try_into().ok()?);
    let height = u32::from_be_bytes(bytes.get(20..24)?.try_into().ok()?);
    Some((width, height))
}

fn render_with_timeout(
    cid: u64,
    smiles: &str,
    mode: &str,
    timeout: Duration,
) -> (bool, ExitStatus, String, String) {
    let output_path = std::env::temp_dir().join(format!(
        "molchemist-pubchem-{cid}-{mode}-{}.typ",
        std::process::id()
    ));
    let mut child = molchemist()
        .args([
            "dump",
            "--smiles",
            smiles,
            "--mode",
            mode,
            "--output",
            output_path.to_str().unwrap(),
        ])
        .stdout(Stdio::null())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let started = Instant::now();
    let mut timed_out = false;

    loop {
        if child.try_wait().unwrap().is_some() {
            break;
        }
        if started.elapsed() >= timeout {
            timed_out = true;
            child.kill().unwrap();
            break;
        }
        thread::sleep(Duration::from_millis(25));
    }

    let output = child.wait_with_output().unwrap();
    let source = fs::read_to_string(&output_path).unwrap_or_default();
    let _ = fs::remove_file(output_path);
    (
        timed_out,
        output.status,
        String::from_utf8_lossy(&output.stderr).into_owned(),
        source,
    )
}

#[test]
#[ignore = "requires scripts/fetch-pubchem-visual-corpus.py"]
fn downloaded_pubchem_smiles_render_and_keep_their_reference_images() {
    let root = corpus_root();
    let manifest_path = root.join("manifest.json");
    let manifest = fs::read_to_string(&manifest_path).unwrap_or_else(|error| {
        panic!(
            "failed to read {}: {error}; run python3 scripts/fetch-pubchem-visual-corpus.py",
            manifest_path.display()
        )
    });
    let requested_cids = requested_cids();
    let include_stress = requested_cids.is_some()
        || std::env::var("MOLCHEMIST_PUBCHEM_INCLUDE_STRESS").is_ok_and(|value| value == "1");
    let cases = serde_json::from_str::<Vec<Value>>(&manifest)
        .unwrap()
        .into_iter()
        .filter(|case| {
            let cid_matches = requested_cids
                .as_ref()
                .is_none_or(|cids| cids.contains(&case["cid"].as_u64().unwrap_or_default()));
            let is_stress = case["stress"].as_bool().unwrap_or(false);
            cid_matches && (include_stress || !is_stress)
        })
        .collect::<Vec<_>>();
    assert!(!cases.is_empty(), "the PubChem corpus is empty");
    let modes = requested_modes();
    assert!(!modes.is_empty(), "no PubChem corpus modes were requested");
    let timeout = Duration::from_secs(
        std::env::var("MOLCHEMIST_PUBCHEM_TIMEOUT_SECS")
            .ok()
            .and_then(|value| value.parse().ok())
            .unwrap_or(30),
    );

    let mut failures = Vec::new();
    for (position, case) in cases.iter().enumerate() {
        let cid = case["cid"].as_u64().expect("CID must be an integer");
        let smiles = case["smiles"].as_str().expect("SMILES must be a string");
        let image = case["image"].as_str().expect("image must be a string");
        let image_bytes = fs::read(root.join(image)).unwrap_or_else(|error| {
            panic!("CID {cid} reference image is missing: {error}");
        });
        let (width, height) = png_dimensions(&image_bytes)
            .unwrap_or_else(|| panic!("CID {cid} reference image is not a valid PNG header"));
        assert!(
            width >= 100 && height >= 100 && image_bytes.len() >= 1024,
            "CID {cid} reference image is unexpectedly small: {width}x{height}, {} bytes",
            image_bytes.len()
        );

        for mode in &modes {
            eprintln!("[{}/{}] CID {cid} ({mode})", position + 1, cases.len());
            let (timed_out, status, stderr, source) =
                render_with_timeout(cid, smiles, mode, timeout);
            if timed_out
                || !status.success()
                || source.contains("NaN")
                || source.contains("inf")
                || !source.contains("\"a0\"")
            {
                failures.push(format!(
                    "CID {cid} ({mode}): timed-out={timed_out} status={status} stderr={} output-prefix={}",
                    stderr.trim(),
                    source.chars().take(160).collect::<String>()
                ));
            }
        }
    }

    assert!(
        failures.is_empty(),
        "PubChem corpus failures:\n{}",
        failures.join("\n")
    );
}
