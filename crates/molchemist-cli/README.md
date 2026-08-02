# molchemist-cli

`molchemist-cli` converts Molfile, SDF, and SMILES input into formatted [`alchemist`](https://typst.app/universe/package/alchemist/) source. The installed executable is named `molchemist`.

The CLI embeds the same Rust parser and Coordgen WebAssembly modules shipped with the molchemist Typst package. With the default `3em` atom separation and two-space indentation, its output is byte-for-byte identical to the package's `dump: true` output.

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

The standalone wrapper defaults to `@preview/alchemist:0.2.0`, a `3mm` page margin, and `3em` atom separation. Override these with `--alchemist-import`, `--page-margin`, and `--atom-sep`.

For a multi-record SDF, select a one-based record with `--record`:

```sh
molchemist dump compounds.sdf --record 3
```

Each selected Molfile/SDF record is detected as V2000 or V3000. Empty structures, malformed records, non-finite coordinates, and out-of-range record numbers are reported as conversion errors.

Extended SDF bond orders are preserved in generated source: aromatic and query bonds use distinct dashed/dotted forms, any and `either` bonds are wavy, coordination bonds retain their arrow direction, hydrogen bonds are dotted, and undefined double-bond geometry is crossed. Wedge/dash bonds to explicit hydrogen remain visible in abbreviated and skeletal modes. SDF atom parity, enhanced stereo groups, and extended OpenSMILES chirality classes are retained as annotations in generated source. The generated helpers are included automatically, including in `--standalone` output.

Disconnected Molfile/SDF graphs and dot-separated SMILES retain every component. Generated Alchemist source places components side by side with a neutral `operator(none, ...)` boundary:

```sh
molchemist dump --smiles '[Na+].[Cl-]' --mode abbreviate
```

The three rendering modes match the Typst API:

- `full` draws every atom represented by the conversion pipeline.
- `abbreviate` folds common hydrogens and terminal groups into labels.
- `skeletal` hides the carbon backbone and attached hydrogens.

`--format auto` first uses an explicit input kind, then the file extension, then the content. Supported extensions are `.mol`, `.sdf`, `.smi`, and `.smiles`. Pass `--format` when piped or extensionless input is ambiguous.

Run `molchemist dump --help` for the complete option list. Generated source is written exclusively to standard output; diagnostics are written to standard error, so shell redirection is safe.

## License

Molchemist-authored CLI source is MIT-licensed. `src/runtime.rs`, which adapts Typst's WebAssembly plugin host, is Apache-2.0-licensed. The crate also embeds the same precompiled WASM components as the Typst package, so its Cargo metadata uses the aggregate SPDX expression `MIT AND BSD-3-Clause AND Apache-2.0 AND (Apache-2.0 WITH LLVM-exception)`.

`wasm-minimal-protocol` is released under the Unlicense and is recorded in the notices rather than the Cargo license expression. See [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md) for the complete file-to-license mapping and the included license files.
