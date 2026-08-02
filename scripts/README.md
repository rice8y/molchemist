# Development Scripts

## Local PubChem visual corpus

`fetch-pubchem-visual-corpus.py` builds a local-only rendering corpus from a diverse CID catalog. It retrieves isomeric SMILES through PubChem PUG REST and the corresponding PubChem 2D PNG, then writes a manifest and a Typst comparison sheet under `.local-tests/pubchem-visual/`.

```sh
python3 scripts/fetch-pubchem-visual-corpus.py
cargo test -p molchemist-cli --test pubchem_corpus -- --ignored --nocapture
typst compile --root . \
  .local-tests/pubchem-visual/comparison.typ \
  .local-tests/pubchem-visual/comparison.pdf
```

Use `--limit 3` for a quick downloader smoke test, `--refresh` to replace cached responses, or `--catalog path/to/cids.tsv` for a larger private catalog. The catalog format is `CID<TAB>category`; comment and blank lines are ignored. The ignored Rust test renders skeletal mode by default, enforces a 30-second per-case timeout, and accepts `MOLCHEMIST_PUBCHEM_MODES=skeletal,abbreviate,full` and `MOLCHEMIST_PUBCHEM_TIMEOUT_SECS=60` for a broader or slower local run. Set `MOLCHEMIST_PUBCHEM_CIDS=14969,392622` to rerun selected records only. Catalog categories prefixed with `stress-` remain in the manifest and visual sheet source but are skipped by default. Set `MOLCHEMIST_PUBCHEM_INCLUDE_STRESS=1` for the Rust run or pass `--input include-stress=true` to `typst compile` to include them.

The generated manifest, PNG files, Typst source, and PDF are ignored by Git and must not be committed. PubChem records incorporate data from many contributors, so users remain responsible for the provenance and licensing restrictions of downloaded content. The script stays below PubChem's documented request-rate limit and retries temporary throttling responses.
