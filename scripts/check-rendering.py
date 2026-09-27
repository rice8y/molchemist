#!/usr/bin/env python3
"""Pixel regression for the distributed synthetic/attributed figure corpus.

Uses Typst's bundled fonts and PNG rasterizer. No OS font or PDF rasterizer is
involved. Update a baseline only after viewing all generated pages.
"""
import argparse
import json
from pathlib import Path
import shutil
import struct
import subprocess
import zlib

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "crates/molchemist-cli/tests/visual"
CURRENT = ROOT / "target/visual-regression"
VERSION = "0.15.1"


def pixels(path):
    data = path.read_bytes()
    if data[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError(f"Not a PNG: {path}")
    offset, payload = 8, bytearray()
    width = height = channels = None
    while offset < len(data):
        size = struct.unpack(">I", data[offset:offset + 4])[0]
        kind = data[offset + 4:offset + 8]
        body = data[offset + 8:offset + 8 + size]
        offset += size + 12
        if kind == b"IHDR":
            width, height, depth, color, _, _, interlace = struct.unpack(">IIBBBBB", body)
            if depth != 8 or color not in (2, 6) or interlace:
                raise ValueError("Expected noninterlaced RGB/RGBA8 PNG")
            channels = 4 if color == 6 else 3
        elif kind == b"IDAT":
            payload.extend(body)
    raw = zlib.decompress(payload)
    stride = width * channels
    previous = bytearray(stride)
    output = bytearray()
    for row in range(height):
        start = row * (stride + 1)
        mode = raw[start]
        current = bytearray(raw[start + 1:start + 1 + stride])
        for i in range(stride):
            left = current[i - channels] if i >= channels else 0
            up = previous[i]
            diagonal = previous[i - channels] if i >= channels else 0
            if mode == 1:
                prediction = left
            elif mode == 2:
                prediction = up
            elif mode == 3:
                prediction = (left + up) // 2
            elif mode == 4:
                p = left + up - diagonal
                distances = (abs(p - left), abs(p - up), abs(p - diagonal))
                prediction = (left, up, diagonal)[distances.index(min(distances))]
            elif mode == 0:
                prediction = 0
            else:
                raise ValueError("Invalid PNG filter")
            current[i] = (current[i] + prediction) & 255
        output.extend(current)
        previous = current
    return (width, height, channels), output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accept", action="store_true")
    args = parser.parse_args()
    version = subprocess.check_output(["typst", "--version"], text=True).split()[1]
    if version != VERSION:
        raise SystemExit(f"Pixel baselines require Typst {VERSION}; found {version}")
    CURRENT.mkdir(parents=True, exist_ok=True)
    # Isolate each run so stale pages cannot satisfy the page-count check.
    import tempfile
    with tempfile.TemporaryDirectory(prefix="run-", dir=CURRENT) as directory:
        generated = Path(directory)
        subprocess.run([
            "typst", "compile", "--root", str(ROOT), "--ignore-system-fonts", "--ppi", "144",
            str(ROOT / "crates/molchemist-cli/tests/fixtures/extended-rendering.typ"),
            str(generated / "page-{p}.png"),
        ], check=True, cwd=ROOT)
        pages = sorted(generated.glob("page-*.png"))
        if not pages:
            raise SystemExit("No rendered pages")
        if args.accept:
            BASE.mkdir(parents=True, exist_ok=True)
            for old in BASE.glob("page-*.png"):
                old.unlink()
            for page in pages:
                shutil.copyfile(page, BASE / page.name)
            (BASE / "manifest.json").write_text(json.dumps({
                "typst": VERSION, "ppi": 144, "systemFonts": False,
                "pages": [page.name for page in pages],
            }, indent=2) + "\n")
            print(f"Accepted {len(pages)} pages")
            return
        manifest = json.loads((BASE / "manifest.json").read_text())
        if [page.name for page in pages] != manifest["pages"]:
            raise SystemExit("Rendered page set differs from the reviewed baseline")
        failures = []
        for page in pages:
            shape, actual = pixels(page)
            expected_shape, expected = pixels(BASE / page.name)
            if shape != expected_shape:
                failures.append(f"{page.name}: dimensions {expected_shape} -> {shape}")
            else:
                differences = [abs(a-b) for a, b in zip(actual, expected)]
                # Permit subpixel rounding on a tiny fraction of samples;
                # moved labels, missing bonds and glyph changes still fail.
                significant = sum(d > 3 for d in differences)
                if significant > len(differences) * 0.0001 or sum(differences)/len(differences) > 0.01:
                    failures.append(f"{page.name}: {significant} changed color samples")
            shutil.copyfile(page, CURRENT / page.name)
        if failures:
            raise SystemExit("\n".join(failures) + f"\nReview new pages in {CURRENT}")
        print(f"{len(pages)} pixel regressions passed")


if __name__ == "__main__":
    main()
