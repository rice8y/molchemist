#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

docs_root="${TYPAGE_DOCS_ROOT:-molchemist-docs}"
jobs="${TYPAGE_BUILD_JOBS:-0}"

case "$docs_root" in
  "" | "." | "/" | /* | *".."*)
    echo "Unsafe TYPAGE_DOCS_ROOT: ${docs_root}" >&2
    exit 1
    ;;
esac

if [ ! -d "$docs_root" ]; then
  echo "Missing ${docs_root}. Run scripts/vercel-install.sh first." >&2
  exit 1
fi

echo "[build] replacing starter package with top-level files from molchemist package"
rm -rf "${docs_root}/package"
mkdir -p "${docs_root}/package"

find package \
  -mindepth 1 \
  -maxdepth 1 \
  -type f \
  -exec cp -p {} "${docs_root}/package/" \;

export PATH="$repo_root/.bin:$HOME/.cargo/bin:/rust/bin:$PATH"

echo "[build] building Typage docs in ${docs_root}"
typage build --root "$docs_root" --force --jobs "$jobs"