#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
expected='MIT AND BSD-3-Clause AND Apache-2.0 AND (Apache-2.0 WITH LLVM-exception)'

package_license="$(sed -n 's/^license = "\(.*\)"$/\1/p' "$root/package/typst.toml")"
cli_license="$(sed -n 's/^license = "\(.*\)"$/\1/p' "$root/crates/molchemist-cli/Cargo.toml")"

if [[ "$package_license" != "$expected" || "$cli_license" != "$expected" ]]; then
  echo "Typst package and CLI must use: $expected" >&2
  echo "Typst package: $package_license" >&2
  echo "CLI: $cli_license" >&2
  exit 1
fi

cmp "$root/LICENSE" "$root/package/LICENSE"
cmp "$root/LICENSE" "$root/crates/molchemist-cli/LICENSE"

for license in \
  LICENSE-APACHE-2.0 \
  LICENSE-EMSCRIPTEN \
  LICENSE-LIBCXX \
  LICENSE-LIBCXXABI \
  LICENSE-MUSL
do
  cmp "$root/package/$license" "$root/crates/molchemist-cli/$license"
done

grep -Fq 'SPDX-License-Identifier: Apache-2.0' \
  "$root/crates/molchemist-cli/src/runtime.rs"

echo "Typst package and CLI license metadata are synchronized."
