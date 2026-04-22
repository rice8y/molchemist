# WASM Plugins

This directory contains the WebAssembly modules and local plugin-specific dependencies used by `molchemist`.

## Layout

- `core/`
  - Rust/WASM plugin for SDF/Molfile parsing, SMILES parsing, AST generation, and Typst-facing exports.
- `smiles-layout/`
  - C++/WASM plugin that runs CoordgenLibs to generate 2D coordinates for SMILES input.
- `vendor/`
  - Vendored plugin dependencies that are patched locally when needed.

## Notes

- The Typst package still ships the same output modules:
  - `molchemist_plugin.wasm`
  - `molchemist_smiles_plugin.wasm`
- `vendor/opensmiles` is intentionally kept in-repo because `molchemist` currently relies on a local parser fix for nested branch and ring-closure handling.
