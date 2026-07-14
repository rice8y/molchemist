#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

cmp "$root/package/molchemist_plugin.wasm" \
  "$root/crates/molchemist-cli/wasm/molchemist_plugin.wasm"
cmp "$root/package/molchemist_smiles_plugin.wasm" \
  "$root/crates/molchemist-cli/wasm/molchemist_smiles_plugin.wasm"

echo "Typst package and CLI WASM binaries are synchronized."
