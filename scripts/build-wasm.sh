#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
layout_build="${MOLCHEMIST_LAYOUT_BUILD_DIR:-$root/target/wasm-plugins/smiles-layout}"
optimized_build="$root/target/wasm-plugins/optimized"
wasm_features=(
  --enable-mutable-globals
  --enable-sign-ext
  --enable-bulk-memory
  --enable-nontrapping-float-to-int
)
emcc_version="$(emcc --version)"
emcc_version="${emcc_version%%$'\n'*}"

if ! command -v wasm-opt >/dev/null 2>&1; then
  echo "wasm-opt is required to build the WASM plugins." >&2
  exit 1
fi

if [[ "$emcc_version" != *" 5.0.5"* ]]; then
  echo "Expected Emscripten 5.0.5, found: $emcc_version" >&2
  echo "Update the bundled runtime notices before changing toolchain versions." >&2
  exit 1
fi

cd "$root"
mkdir -p "$optimized_build"

cargo build -p molchemist_plugin --target wasm32-unknown-unknown --release --locked
wasm-opt \
  -Os \
  "${wasm_features[@]}" \
  target/wasm32-unknown-unknown/release/molchemist_plugin.wasm \
  -o "$optimized_build/molchemist_plugin.wasm"
cp "$optimized_build/molchemist_plugin.wasm" \
  package/molchemist_plugin.wasm
cp "$optimized_build/molchemist_plugin.wasm" \
  crates/molchemist-cli/wasm/molchemist_plugin.wasm

emcmake cmake \
  -S wasm-plugins/smiles-layout \
  -B "$layout_build" \
  -G Ninja \
  -DCMAKE_BUILD_TYPE=Release
cmake --build "$layout_build" --target molchemist_smiles_plugin

wasm-opt \
  -Os \
  "${wasm_features[@]}" \
  "$layout_build/dist/molchemist_smiles_plugin.wasm" \
  -o "$optimized_build/molchemist_smiles_plugin.wasm"
cp "$optimized_build/molchemist_smiles_plugin.wasm" \
  package/molchemist_smiles_plugin.wasm
cp "$optimized_build/molchemist_smiles_plugin.wasm" \
  crates/molchemist-cli/wasm/molchemist_smiles_plugin.wasm

"$root/scripts/check-wasm-sync.sh"
