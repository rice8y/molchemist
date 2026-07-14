# Core WASM Plugin

This Rust crate builds `molchemist_plugin.wasm`.

It exposes the shared `molchemist-core` engine through Typst's minimal plugin protocol. The exports cover:

- parsing Molfile/SDF input
- parsing SMILES input
- preparing layout payloads for the companion SMILES layout plugin
- converting molecule graphs into the CBOR AST consumed by the Typst package
- formatting the same graphs as editable Alchemist source

The parsing and formatting implementation is in `../../crates/molchemist-core`. This crate deliberately stays a small ABI wrapper so the Typst package and CLI cannot drift into separate conversion implementations.

Build from the repository root with:

```sh
cargo build -p molchemist_plugin --target wasm32-unknown-unknown --release
```
