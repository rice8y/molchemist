# Core WASM Plugin

This Rust crate builds `molchemist_plugin.wasm`.

It is responsible for:

- parsing Molfile/SDF input
- parsing SMILES input
- preparing layout payloads for the companion SMILES layout plugin
- converting molecule graphs into the CBOR AST consumed by the Typst package

The crate depends on the vendored `opensmiles` copy in `../vendor/opensmiles` because `molchemist` currently needs a local parser patch that is not available in the upstream crate as-is.
