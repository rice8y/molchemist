# molchemist-cli

`molchemist-cli` converts Molfile, SDF, SMILES, MOL2, RXN, Reaction SMILES, and RGfile input into formatted [`alchemist`](https://typst.app/universe/package/alchemist/) source. The installed executable is named `molchemist`.

The CLI embeds the same Rust parser and Coordgen WebAssembly modules shipped with the molchemist Typst package. The package and CLI share one coordinate renderer, and equivalent configurations produce identical `dump` output.

## Install

Rust 1.86 or later is required.

```sh
cargo install --locked molchemist-cli
```

No JavaScript runtime or system chemistry library is required.

## Compatibility

| Concern | Requirement or CI coverage |
| --- | --- |
| Building and installing the CLI | Rust 1.86 or later |
| CLI platforms | Rust 1.86 tests on Ubuntu, macOS, and Windows |
| Compiling generated Typst source | Typst 0.14.0, 0.14.1, 0.14.2, 0.15.0, and 0.15.1 |
| Package/CLI output parity | Byte-for-byte checked on Ubuntu for every listed Typst version |

The CLI itself does not invoke Typst unless you separately compile its generated source. Compatibility with Typst 0.14.0 and 0.14.1 is tested, but upstream recommends 0.14.2 or later because those earlier releases contain a [WebAssembly runtime security issue](https://github.com/typst/typst/releases/tag/v0.14.2).

## Usage

Dump a Molfile or SDF file to standard output:

```sh
molchemist dump molecule.sdf
molchemist dump molecule.mol --mode abbreviate
```

Convert a SMILES string:

```sh
molchemist dump --smiles 'CC(=O)Oc1ccccc1C(=O)O' --mode skeletal
```

SMILES input is parsed strictly. Malformed branch, dot, bond, bracket-property, charge, isotope, atom-class, directional-bond, and aromatic notation returns an error instead of being normalized silently.

Input can also come from standard input or `--text`:

```sh
printf '%s\n' 'c1ccccc1' | molchemist dump --format smiles
molchemist dump --text 'c1ccncc1' --format smiles
```

Use `--output` to write the generated source to a file. Add `--standalone` to create a directly compilable Typst document with the current Alchemist import and an auto-sized page:

```sh
molchemist dump \
  --smiles 'CC(=O)O' \
  --mode skeletal \
  --standalone \
  --output acetic-acid.typ

typst compile acetic-acid.typ
```

Generated source imports `@preview/alchemist:0.2.0`. The standalone wrapper sets a `3mm` page margin. Use `--page-margin` and `--atom-sep` to change the page margin and bond-length unit.

For a multi-record SDF, select a one-based record with `--record`:

```sh
molchemist dump compounds.sdf --record 3
```

Inspect the selected record as semantic JSON when downstream processing needs data that is not part of the drawing:

```sh
molchemist inspect compounds.sdf --record 3
```

Inspection preserves the header, stable source atom/bond IDs, enhanced stereo groups, SGroups with original type codes, complete COLLECTION names and member lists, normalized V2000/V3000 query metadata, duplicate and multiline SDF properties in source order, the exact selected record, and diagnostics for remaining depiction gaps. V3000 `sourceId` values are the original CTAB IDs; V2000 uses one-based source positions. Pass `--compact` for single-line JSON.

`dump` uses `--fidelity warn` by default. Warnings go to standard error without contaminating generated source on standard output. Use `--fidelity strict` to reject any known unsupported depiction feature or `--fidelity ignore` to omit those diagnostics.

Generated source depicts CTfile atom lists and query constraints, R-group labels, ring/chain bond topology, SGroup brackets and labels, contracted multi-atom superatoms, variable-attachment bonds, link nodes, and V3000 atom/bond HILITE collections. Contracted superatoms reconnect crossing bonds at a labelled graph node, use explicit SAP atoms for placement, and project highlight membership onto the contracted glyph. Known SGroup HILITE membership projects onto its atoms and bonds. User-defined collections and unresolved 3D-object, external R-group or generic members remain diagnostic in strict mode. The same overlay code is included in `--standalone` output. Standalone-molecule reaction-center query flags remain diagnostic; reaction rendering highlights graph differences across mapped reactants and products.

Highlight behavior has dedicated semantic and SVG-shape regressions. The corpus separates atom-only, long-query-glyph, bond-only, connected, and disconnected selections, and also includes RDKit/Bingo's real `v3k.crash1.mol`. Run them with:

```sh
cargo test -p molchemist-cli --test cli highlight_
cargo test -p molchemist-cli --test cli local_typst_package_renders_highlight_shape_corpus_to_svg -- --ignored --exact
```

The same five cases are shown with their source directly in `package/docs/documentation.typ`.

Each selected Molfile/SDF record is detected as V2000 or V3000. Empty structures, malformed records, non-finite coordinates, and out-of-range record numbers are reported as conversion errors.

Extended SDF bond orders are preserved in generated source: aromatic and query bonds use distinct dashed/dotted forms, any and `either` bonds are wavy, coordination bonds retain their arrow direction, hydrogen bonds are dotted, and undefined double-bond geometry is crossed. Wedge/dash bonds to explicit hydrogen remain visible in abbreviated and skeletal modes. SDF atom parity, enhanced stereo groups, and extended OpenSMILES chirality classes are retained as annotations in generated source. The generated helpers are included automatically, including in `--standalone` output.

Disconnected Molfile/SDF graphs and dot-separated SMILES retain every component. The renderer packs disconnected components side by side; `--components preserve` retains their source positions:

```sh
molchemist dump --smiles '[Na+].[Cl-]' --mode abbreviate
```

The three rendering modes match the Typst API:

- `full` draws every atom represented by the conversion pipeline.
- `abbreviate` folds common hydrogens and terminal groups into labels.
- `skeletal` hides the carbon backbone and attached hydrogens.

`--format auto` first uses an explicit input kind, then the file extension, then the content. Supported extensions are `.mol`, `.sdf`, `.smi`, and `.smiles`. Pass `--format` when piped or extensionless input is ambiguous.

Run `molchemist dump --help` or `molchemist inspect --help` for the complete option list. Generated source and inspection JSON are written exclusively to standard output; diagnostics are written to standard error, so shell redirection is safe.

## Layout, reactions, and R-groups

```sh
molchemist dump molecule.sdf --layout avoid --strict-collisions --standalone
molchemist dump molecule.mol --layout reflow --infer-stereo --standalone
molchemist dump molecule.mol2 --layout coordinates --components preserve
molchemist dump alternatives.mol --format rgroup --fidelity strict
molchemist dump reaction.rxn --conditions oxidation --standalone
molchemist inspect --text 'CCO>>CC=O' --infer-mapping
```

`--layout` accepts `avoid` (default), `coordinates`, and `reflow`. The coordinate renderer measures actual font bounds and clips bonds at labels. `avoid` increases spacing; `reflow` also regenerates CTfile coordinates. `--strict-collisions` rejects remaining measured atom-label collisions when the output is compiled; the default layout is `avoid`. `--expand-superatoms` depicts the source atoms of SUP groups.

MOL2 source text is preserved as a JSON string in the `MOL2_SOURCE` property. Partial charges are retained without rounding them into formal charges. RGfile inspection preserves `rawSource` and all member alternatives. Reaction mapping inference is optional, respects supplied maps, and reports ambiguity and search-limit exhaustion. `--fidelity strict` rejects those cases and unsupported CTfile features in reaction and R-group components. Use `--mapping-search-limit` to change the mapping search budget.

The [manual](../../package/docs/documentation.pdf) documents layout configuration, supported stereochemistry, input formats, and inference limits.

## License

Molchemist-authored CLI source is MIT-licensed. `src/runtime.rs`, which adapts Typst's WebAssembly plugin host, is Apache-2.0-licensed. The crate also embeds the same precompiled WASM components as the Typst package, so its Cargo metadata uses the aggregate SPDX expression `MIT AND BSD-3-Clause AND Apache-2.0 AND (Apache-2.0 WITH LLVM-exception)`.

`wasm-minimal-protocol` is released under the Unlicense and is recorded in the notices rather than the Cargo license expression. See [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md) for the complete file-to-license mapping and the included license files.
