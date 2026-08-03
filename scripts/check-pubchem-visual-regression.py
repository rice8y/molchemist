#!/usr/bin/env python3
"""Compile and compare the local PubChem side-by-side visual corpus."""

from __future__ import annotations

import argparse
import hashlib
import shutil
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    repo_root = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--corpus",
        type=Path,
        default=repo_root / ".local-tests/pubchem-visual",
        help="directory produced by fetch-pubchem-visual-corpus.py",
    )
    parser.add_argument(
        "--accept",
        action="store_true",
        help="replace the local raster baseline with the current rendering",
    )
    parser.add_argument(
        "--include-stress",
        action="store_true",
        help="include stress-large-molecule records",
    )
    parser.add_argument("--dpi", type=int, default=144, help="rasterization resolution")
    args = parser.parse_args()
    if args.dpi < 72:
        parser.error("--dpi must be at least 72")
    return args


def require_program(name: str) -> str:
    path = shutil.which(name)
    if path is None:
        raise RuntimeError(f"required program not found: {name}")
    return path


def clear_pngs(directory: Path) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    for path in directory.glob("*.png"):
        path.unlink()


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def render(repo_root: Path, corpus: Path, dpi: int, include_stress: bool) -> list[Path]:
    source = corpus / "comparison.typ"
    manifest = corpus / "manifest.json"
    if not source.is_file() or not manifest.is_file():
        raise RuntimeError(
            f"missing local corpus under {corpus}; run fetch-pubchem-visual-corpus.py"
        )

    visual_root = corpus / "visual-regression"
    current_pdf = visual_root / "current.pdf"
    current_pages = visual_root / "current-pages"
    visual_root.mkdir(parents=True, exist_ok=True)
    clear_pngs(current_pages)

    typst = require_program("typst")
    pdftoppm = require_program("pdftoppm")
    subprocess.run(
        [
            typst,
            "compile",
            "--root",
            str(repo_root),
            "--input",
            f"include-stress={'true' if include_stress else 'false'}",
            str(source),
            str(current_pdf),
        ],
        cwd=repo_root,
        check=True,
    )
    subprocess.run(
        [
            pdftoppm,
            "-png",
            "-r",
            str(dpi),
            str(current_pdf),
            str(current_pages / "page"),
        ],
        check=True,
    )
    pages = sorted(current_pages.glob("page-*.png"))
    if not pages:
        raise RuntimeError("comparison PDF produced no raster pages")
    print(f"rendered {len(pages)} pages to {current_pdf}")
    return pages


def main() -> int:
    args = parse_args()
    repo_root = Path(__file__).resolve().parent.parent
    corpus = args.corpus.resolve()
    current_pages = render(repo_root, corpus, args.dpi, args.include_stress)
    baseline_dir = corpus / "visual-regression" / "baseline-pages"

    if args.accept:
        clear_pngs(baseline_dir)
        for page in current_pages:
            shutil.copy2(page, baseline_dir / page.name)
        print(f"accepted {len(current_pages)} baseline pages under {baseline_dir}")
        return 0

    baseline_pages = sorted(baseline_dir.glob("page-*.png"))
    if not baseline_pages:
        raise RuntimeError(
            "visual baseline is missing; inspect visual-regression/current.pdf, then rerun with --accept"
        )

    current_by_name = {page.name: page for page in current_pages}
    baseline_by_name = {page.name: page for page in baseline_pages}
    missing = sorted(baseline_by_name.keys() - current_by_name.keys())
    added = sorted(current_by_name.keys() - baseline_by_name.keys())
    changed = sorted(
        name
        for name in current_by_name.keys() & baseline_by_name.keys()
        if digest(current_by_name[name]) != digest(baseline_by_name[name])
    )
    if missing or added or changed:
        details = []
        if missing:
            details.append(f"missing pages: {', '.join(missing)}")
        if added:
            details.append(f"added pages: {', '.join(added)}")
        if changed:
            details.append(f"changed pages: {', '.join(changed)}")
        raise RuntimeError("visual regression detected; " + "; ".join(details))

    print(f"visual regression passed for {len(current_pages)} pages")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
