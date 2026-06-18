# WASM Plugins

This directory contains the WebAssembly modules and local plugin-specific dependencies used by `molchemist`.

## Layout

- `core/`
  - Rust/WASM plugin for SDF/Molfile parsing, SMILES parsing, AST generation, and Typst-facing exports.
- `smiles-layout/`
  - C++/WASM plugin that runs CoordgenLibs to generate 2D coordinates for SMILES input.
- `vendor/`
  - Vendored plugin dependencies that are maintained locally for the Typst/WASM workflow.

## Notes

- The Typst package still ships the same output modules:
  - `molchemist_plugin.wasm`
  - `molchemist_smiles_plugin.wasm`
- `vendor/opensmiles` is intentionally kept in-repo because `molchemist` relies on local SMILES parser behavior and fixes. Treat it as a maintained vendored dependency for this package rather than a temporary patch queue awaiting upstream synchronization. Local changes are documented in `vendor/opensmiles/PATCHES.md`.
