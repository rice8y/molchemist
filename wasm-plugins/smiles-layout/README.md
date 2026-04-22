# SMILES Layout WASM Plugin

This directory contains the C++/WASM companion plugin used by `render-smiles(...)`.

Its role is intentionally small:

`SMILES (parsed in Rust) -> atom/bond payload -> CoordgenLibs -> 2D coordinates`

The Rust plugin in `../core` handles OpenSMILES parsing, aromatic bond kekulization, and AST generation.
This plugin only reads the compact layout payload, runs CoordgenLibs, and returns packed 2D coordinates.

## Build notes

- `vendor/coordgenlibs` is copied from the `coordgen` crate source so the Typst plugin build stays self-contained.
- The output module is `molchemist_smiles_plugin.wasm`.
- The module uses Typst's `wasm_minimal_protocol` imports and is intended to be called from Typst, not from JavaScript.
