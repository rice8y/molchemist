#!/usr/bin/env python3
"""Optimize distributed plugins after comparing their APIs, output and rendering."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
PLUGINS = ("molchemist_plugin.wasm", "molchemist_smiles_plugin.wasm")
FIXTURES = Path("crates/molchemist-cli/tests/fixtures")
FEATURES = (
    "--enable-mutable-globals",
    "--enable-sign-ext",
    "--enable-bulk-memory",
    "--enable-nontrapping-float-to-int",
)


def run(*args, **kwargs):
    subprocess.run(args, cwd=ROOT, check=True, **kwargs)


def replace(path, data):
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as stream:
        temporary = Path(stream.name)
        stream.write(data)
    try:
        shutil.copymode(path, temporary)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def main():
    for tool in ("wasm-opt", "cargo", "typst"):
        if shutil.which(tool) is None:
            raise SystemExit(f"{tool} must be installed to optimize and verify the plugins")
    originals = {}
    for name in PLUGINS:
        package = ROOT / "package" / name
        cli = ROOT / "crates/molchemist-cli/wasm" / name
        originals[package] = package.read_bytes()
        originals[cli] = cli.read_bytes()
        if originals[package] != originals[cli]:
            raise SystemExit(f"Package and CLI copies differ: {name}; synchronize them first")

    with tempfile.TemporaryDirectory(prefix="molchemist-wasm-") as directory:
        work = Path(directory)
        for variant in ("before", "after"):
            root = work / variant
            (root / "package").mkdir(parents=True)
            shutil.copy2(ROOT / "package/lib.typ", root / "package/lib.typ")
            shutil.copytree(ROOT / "package/src", root / "package/src")
            shutil.copytree(ROOT / FIXTURES, root / FIXTURES)
            for name in PLUGINS:
                (root / "package" / name).write_bytes(originals[ROOT / "package" / name])
        for name in PLUGINS:
            run("wasm-opt", "-Os", *FEATURES,
                str(work / "before/package" / name), "-o", str(work / "after/package" / name))
            candidate = work / "after/package" / name
            original = originals[ROOT / "package" / name]
            if candidate.stat().st_size > len(original):
                print(f"{name}: optimized output is larger; keeping the original plugin", flush=True)
                candidate.write_bytes(original)

        print("Checking plugin interfaces, generated source, inspection and errors...", flush=True)
        environment = os.environ.copy()
        environment["MOLCHEMIST_OPTIMIZED_WASM_DIR"] = str(work / "after/package")
        run("cargo", "test", "--locked", "-p", "molchemist-cli", "--lib",
            "runtime::optimization_tests::optimized_wasm_preserves_behavior", "--", "--ignored", "--exact",
            env=environment)

        print("Comparing rendered figures before and after optimization...", flush=True)
        for variant in ("before", "after"):
            root = work / variant
            (root / "images").mkdir()
            run("typst", "compile", "--root", str(root), "--ignore-system-fonts", "--ppi", "144",
                str(root / FIXTURES / "extended-rendering.typ"), str(root / "images/page-{p}.png"))
        before = {p.name: p.read_bytes() for p in (work / "before/images").glob("*.png")}
        after = {p.name: p.read_bytes() for p in (work / "after/images").glob("*.png")}
        if not before or before != after:
            raise SystemExit("Optimization changed the rendered figures; original plugins were retained")
        print(f"All {len(before)} rendered pages match exactly.", flush=True)

        # Do not overwrite a plugin that changed while verification was running.
        for path, data in originals.items():
            if path.read_bytes() != data:
                raise SystemExit(f"{path} changed during verification; no plugins were replaced")
        written = []
        try:
            for name in PLUGINS:
                data = (work / "after/package" / name).read_bytes()
                original = originals[ROOT / "package" / name]
                for path in (ROOT / "package" / name, ROOT / "crates/molchemist-cli/wasm" / name):
                    if data != originals[path]:
                        written.append(path)
                        replace(path, data)
                print(f"{name}: {len(original)} -> {len(data)} bytes")
        except BaseException:
            for path in written:
                replace(path, originals[path])
            raise
    print("Optimized plugins verified; package and CLI copies are synchronized.")


if __name__ == "__main__":
    main()
