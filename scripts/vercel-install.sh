#!/usr/bin/env bash
set -euo pipefail

typst_version="${TYPST_VERSION:-0.15.0}"
typage_version="${TYPAGE_VERSION:-0.1.5}"
starter_version="${TYPST_DOCS_STARTER_VERSION:-0.1.3}"
starter="${TYPST_DOCS_STARTER:-github:rice8y/typage-starter-typst-docs#v${starter_version}}"
docs_root="${TYPAGE_DOCS_ROOT:-molchemist-docs}"

case "$docs_root" in
  "" | "." | "/" | /* | *".."*)
    echo "Unsafe TYPAGE_DOCS_ROOT: ${docs_root}" >&2
    exit 1
    ;;
esac

platform="$(uname -s):$(uname -m)"
case "$platform" in
  Linux:x86_64 | Linux:amd64)
    typst_target="x86_64-unknown-linux-musl"
    ;;
  Linux:aarch64 | Linux:arm64)
    typst_target="aarch64-unknown-linux-musl"
    ;;
  Darwin:x86_64)
    typst_target="x86_64-apple-darwin"
    ;;
  Darwin:aarch64 | Darwin:arm64)
    typst_target="aarch64-apple-darwin"
    ;;
  *)
    echo "Unsupported platform: ${platform}" >&2
    exit 1
    ;;
esac

typst_dir=".cache/typst/v${typst_version}"
typst_archive="${typst_dir}/typst-${typst_target}.tar.xz"
typst_url="https://github.com/typst/typst/releases/download/v${typst_version}/typst-${typst_target}.tar.xz"

mkdir -p .bin "$typst_dir"

echo "[setup] installing Typst ${typst_version} for ${typst_target}"
curl -fsSL "$typst_url" -o "$typst_archive"
tar -xJf "$typst_archive" -C "$typst_dir"
cp "${typst_dir}/typst-${typst_target}/typst" .bin/typst
chmod +x .bin/typst

echo "[setup] installing Typage ${typage_version}"
cargo install typage --version "$typage_version" --locked

export PATH="$PWD/.bin:$HOME/.cargo/bin:/rust/bin:$PATH"

echo "[setup] initializing Typage docs starter from ${starter}"
rm -rf "$docs_root"
typage init "$docs_root" \
  --starter "$starter" \
  --install-plugins
