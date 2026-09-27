#!/usr/bin/env bash
set -euo pipefail
# Build the Rust plugin and synchronize the package and CLI distributions.
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
cargo build -p molchemist_plugin --target wasm32-unknown-unknown --release --locked
mkdir -p target/wasm-plugins/optimized
wasm-opt -Os --enable-mutable-globals --enable-sign-ext --enable-bulk-memory --enable-nontrapping-float-to-int \
  target/wasm32-unknown-unknown/release/molchemist_plugin.wasm \
  -o target/wasm-plugins/optimized/molchemist_plugin.wasm
cp target/wasm-plugins/optimized/molchemist_plugin.wasm package/molchemist_plugin.wasm
cp target/wasm-plugins/optimized/molchemist_plugin.wasm crates/molchemist-cli/wasm/molchemist_plugin.wasm
bash scripts/check-wasm-sync.sh
