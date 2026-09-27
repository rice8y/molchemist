# Development Scripts

## Local PubChem visual corpus

`fetch-pubchem-visual-corpus.py` builds a local-only rendering corpus from a diverse CID catalog. It retrieves isomeric SMILES through PubChem PUG REST and the corresponding PubChem 2D PNG, then writes a manifest and a Typst comparison sheet under `.local-tests/pubchem-visual/`.

```sh
python3 scripts/fetch-pubchem-visual-corpus.py
cargo test -p molchemist-cli --test pubchem_corpus -- --ignored --nocapture
typst compile --root . \
  .local-tests/pubchem-visual/comparison.typ \
  .local-tests/pubchem-visual/comparison.pdf
python3 scripts/check-pubchem-visual-regression.py
```

Use `--limit 3` for a quick downloader smoke test, `--refresh` to replace cached responses, or `--catalog path/to/cids.tsv` for a larger private catalog. The catalog format is `CID<TAB>category`; comment and blank lines are ignored. The ignored Rust test renders skeletal mode by default, enforces a 30-second per-case timeout, and accepts `MOLCHEMIST_PUBCHEM_MODES=skeletal,abbreviate,full` and `MOLCHEMIST_PUBCHEM_TIMEOUT_SECS=60` for a broader or slower local run. Set `MOLCHEMIST_PUBCHEM_CIDS=14969,392622` to rerun selected records only. Catalog categories prefixed with `stress-` remain in the manifest and visual sheet source but are skipped by default. Set `MOLCHEMIST_PUBCHEM_INCLUDE_STRESS=1` for the Rust run or pass `--input include-stress=true` to `typst compile` to include them.

The generated manifest, PNG files, Typst source, and PDF are ignored by Git and must not be committed. PubChem records incorporate data from many contributors, so users remain responsible for the provenance and licensing restrictions of downloaded content. The script stays below PubChem's documented request-rate limit and retries temporary throttling responses.

`check-pubchem-visual-regression.py` compiles the side-by-side PubChem/molchemist comparison sheet, rasterizes it with `pdftoppm`, and compares exact page hashes against a local ignored baseline. On the first run, inspect `.local-tests/pubchem-visual/visual-regression/current.pdf`, then run `python3 scripts/check-pubchem-visual-regression.py --accept`; subsequent runs fail when pages are added, removed, or changed. Baselines are intentionally machine-local because Typst, Poppler, and font updates can alter raster output without changing molecular semantics.

## Rust plugin and figure regression checks

When only Rust or shared Typst renderer code has changed, rebuild and synchronize the core WASM module without rebuilding the unchanged Coordgen C++ module:

```sh
bash scripts/build-core-wasm.sh
```

Use `scripts/build-wasm.sh` with its pinned Emscripten version for C++ changes.

To optimize the existing distributed plugins without rebuilding them, run `just --justfile package/justfile optimize-wasm` from the repository root. This requires Python 3, Cargo, Typst with Alchemist 0.2.0, and Binaryen's `wasm-opt`. The task uses `wasm-opt -Os` with the same WebAssembly features as the build scripts, compares plugin interfaces, generated source, inspection results, error messages and rendered figures, then updates both the Typst and CLI copies. Failed checks leave the original plugins intact; larger optimized files are not installed. These regression checks cover the bundled examples rather than proving equivalence for every possible input.

`python3 scripts/check-rendering.py` checks committed, reviewed figures with Typst 0.15.1 and bundled fonts. Inspect changed images under `target/visual-regression` before running `--accept` to update the baseline. The script compares decoded pixel samples with a small antialiasing tolerance.

## Documentation comparison figures

`package/docs/examples/` contains the executable sources shared by the README comparisons and manual. To regenerate a PNG with Typst 0.15.1 and bundled fonts, run the following command from the repository root; use `reaction`, `superatoms`, or `rgroups` for the other figures. The matching code blocks in both READMEs use the same source, with asset paths relative to each README directory.

```sh
typst compile --root . --ignore-system-fonts --ppi 192 \
  --input figure=layout package/docs/render-comparison.typ \
  package/images/comparison-layout.png
```
