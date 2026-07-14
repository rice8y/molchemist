#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
layout_build="${MOLCHEMIST_LAYOUT_BUILD_DIR:-$root/target/wasm-plugins/smiles-layout}"
emcc_version="$(emcc --version)"
emcc_version="${emcc_version%%$'\n'*}"

if [[ "$emcc_version" != *" 5.0.5"* ]]; then
  echo "Expected Emscripten 5.0.5, found: $emcc_version" >&2
  echo "Update the bundled runtime notices before changing toolchain versions." >&2
  exit 1
fi

cd "$root"

cargo build -p molchemist_plugin --target wasm32-unknown-unknown --release --locked
cp target/wasm32-unknown-unknown/release/molchemist_plugin.wasm \
  package/molchemist_plugin.wasm
cp target/wasm32-unknown-unknown/release/molchemist_plugin.wasm \
  crates/molchemist-cli/wasm/molchemist_plugin.wasm

emcmake cmake \
  -S wasm-plugins/smiles-layout \
  -B "$layout_build" \
  -G Ninja \
  -DCMAKE_BUILD_TYPE=Release
cmake --build "$layout_build" --target molchemist_smiles_plugin

cp "$layout_build/dist/molchemist_smiles_plugin.wasm" \
  package/molchemist_smiles_plugin.wasm
cp "$layout_build/dist/molchemist_smiles_plugin.wasm" \
  crates/molchemist-cli/wasm/molchemist_smiles_plugin.wasm

"$root/scripts/check-wasm-sync.sh"
