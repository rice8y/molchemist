# WASM Plugins

This directory contains the thin WebAssembly adapters used by `molchemist`. Shared conversion code and vendored chemistry dependencies live in `../crates/molchemist-core`.

## Layout

- `core/`
  - Rust/WASM exports for the shared `molchemist-core` crate.
- `smiles-layout/`
  - C++/WASM adapter that runs the shared Coordgen engine.

## Notes

- The Typst package still ships the same output modules:
  - `molchemist_plugin.wasm`
  - `molchemist_smiles_plugin.wasm`
- `../crates/molchemist-core/vendor/opensmiles` is intentionally maintained in-repo because molchemist relies on local parser behavior and fixes. Local changes are documented in its `PATCHES.md`.
- `../crates/molchemist-core/vendor/coordgenlibs` is shared by the native test bridge and the C++/WASM layout adapter.
- The CLI embeds byte-identical copies of both package WASM files under `../crates/molchemist-cli/wasm`.

## Build and synchronize

Install the Rust `wasm32-unknown-unknown` target, Emscripten 5.0.5, CMake, Ninja, and `wasi-stub`, then run from the repository root:

```sh
./scripts/build-wasm.sh
```

The script rebuilds both plugins, copies each byte-identical result to the Typst package and CLI crate, and verifies that the two destinations match.
