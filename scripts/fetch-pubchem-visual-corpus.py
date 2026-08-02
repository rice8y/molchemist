#!/usr/bin/env python3
"""Build a local-only PubChem SMILES and 2D-image visual corpus."""

from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path


PUBCHEM_PROPERTY_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
    "{cids}/property/Title,IsomericSMILES,IUPACName/JSON"
)
PUBCHEM_IMAGE_URL = "https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={cid}&t=l"
PUBCHEM_COMPOUND_URL = "https://pubchem.ncbi.nlm.nih.gov/compound/{cid}#section=SMILES"
PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"
USER_AGENT = "molchemist-local-visual-corpus/0.1 (https://github.com/rice8y/molchemist)"


class PubChemClient:
    def __init__(self, timeout: float = 30.0, minimum_interval: float = 0.26) -> None:
        self.timeout = timeout
        self.minimum_interval = minimum_interval
        self.last_request = 0.0

    def get(self, url: str) -> bytes:
        for attempt in range(5):
            wait = self.minimum_interval - (time.monotonic() - self.last_request)
            if wait > 0:
                time.sleep(wait)
            request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            self.last_request = time.monotonic()
            try:
                with urllib.request.urlopen(request, timeout=self.timeout) as response:
                    return response.read()
            except urllib.error.HTTPError as error:
                if error.code not in {429, 500, 502, 503, 504} or attempt == 4:
                    raise
            except urllib.error.URLError:
                if attempt == 4:
                    raise
            time.sleep(2**attempt)
        raise RuntimeError(f"failed to retrieve {url}")


def parse_args() -> argparse.Namespace:
    repo_root = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--catalog",
        type=Path,
        default=repo_root / "scripts/pubchem-visual-cids.tsv",
        help="CID/category TSV catalog",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=repo_root / ".local-tests/pubchem-visual",
        help="local output directory",
    )
    parser.add_argument(
        "--limit",
        type=int,
        help="only process the first N catalog entries",
    )
    parser.add_argument(
        "--refresh",
        action="store_true",
        help="replace already downloaded PNG files",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=30.0,
        help="HTTP timeout in seconds",
    )
    args = parser.parse_args()
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be positive")
    return args


def read_catalog(path: Path) -> list[tuple[int, str]]:
    entries: list[tuple[int, str]] = []
    seen: set[int] = set()
    for line_number, raw_line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 2:
            raise ValueError(f"{path}:{line_number}: expected CID<TAB>category")
        cid = int(fields[0])
        category = fields[1].strip()
        if cid <= 0 or not category:
            raise ValueError(f"{path}:{line_number}: invalid CID or category")
        if cid in seen:
            raise ValueError(f"{path}:{line_number}: duplicate CID {cid}")
        seen.add(cid)
        entries.append((cid, category))
    if not entries:
        raise ValueError(f"{path}: catalog is empty")
    return entries


def chunks(values: list[int], size: int = 100) -> list[list[int]]:
    return [values[offset : offset + size] for offset in range(0, len(values), size)]


def fetch_properties(client: PubChemClient, cids: list[int]) -> dict[int, dict[str, object]]:
    properties: dict[int, dict[str, object]] = {}
    for cid_chunk in chunks(cids):
        url = PUBCHEM_PROPERTY_URL.format(cids=",".join(map(str, cid_chunk)))
        payload = json.loads(client.get(url))
        for record in payload.get("PropertyTable", {}).get("Properties", []):
            properties[int(record["CID"])] = record
    missing = [cid for cid in cids if cid not in properties]
    if missing:
        raise RuntimeError(f"PubChem did not return metadata for CIDs: {missing}")
    return properties


def write_comparison_source(path: Path) -> None:
    path.write_text(
        '''#import "../../package/lib.typ": render-smiles

#set page(paper: "a4", margin: 10mm)
#set text(size: 8pt)

#let include-stress = sys.inputs.at("include-stress", default: "false") == "true"
#let cases = json("manifest.json").filter(case => include-stress or not case.stress)
#let card(case) = block(
  breakable: false,
  inset: 6pt,
  radius: 3pt,
  stroke: 0.5pt + luma(180),
)[
  #text(weight: "bold")[CID #case.cid · #case.title]
  #h(0.6em)
  #text(fill: luma(90))[#case.category]
  #v(3pt)
  #grid(
    columns: (1fr, 1fr),
    gutter: 8pt,
    align: center + horizon,
    [
      #image(case.image, width: 100%, height: 46mm, fit: "contain")
      #align(center)[#text(size: 6pt, fill: luma(90))[PubChem reference]]
    ],
    [
      #box(
        width: 100%,
        height: 46mm,
        align(center + horizon, scale(
          55%,
          origin: center,
          render-smiles(case.smiles, skeletal: true),
        )),
      )
      #align(center)[#text(size: 6pt, fill: luma(90))[molchemist skeletal]]
    ],
  )
  #v(2pt)
  #text(size: 6pt)[#case.smiles]
]

#for case in cases {
  card(case)
  v(5pt)
}
''',
        encoding="utf-8",
    )


def main() -> int:
    args = parse_args()
    entries = read_catalog(args.catalog)
    if args.limit is not None:
        entries = entries[: args.limit]

    output = args.output.resolve()
    image_dir = output / "images"
    image_dir.mkdir(parents=True, exist_ok=True)

    client = PubChemClient(timeout=args.timeout)
    properties = fetch_properties(client, [cid for cid, _ in entries])
    manifest = []

    for position, (cid, category) in enumerate(entries, 1):
        record = properties[cid]
        smiles = record.get("SMILES") or record.get("IsomericSMILES")
        if not isinstance(smiles, str) or not smiles:
            raise RuntimeError(f"CID {cid} did not provide an isomeric SMILES")

        image_path = image_dir / f"{cid}.png"
        if args.refresh or not image_path.exists():
            image = client.get(PUBCHEM_IMAGE_URL.format(cid=cid))
            if not image.startswith(PNG_SIGNATURE):
                raise RuntimeError(f"CID {cid} image response is not a PNG")
            image_path.write_bytes(image)

        title = record.get("Title") or record.get("IUPACName") or f"CID {cid}"
        manifest.append(
            {
                "cid": cid,
                "category": category,
                "stress": category.startswith("stress-"),
                "title": title,
                "iupac_name": record.get("IUPACName", ""),
                "smiles": smiles,
                "compound_url": PUBCHEM_COMPOUND_URL.format(cid=cid),
                "image_url": PUBCHEM_IMAGE_URL.format(cid=cid),
                "image": f"images/{cid}.png",
            }
        )
        print(f"[{position:>3}/{len(entries)}] CID {cid}: {title}")

    (output / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    write_comparison_source(output / "comparison.typ")
    print(f"wrote {len(manifest)} cases to {output}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError, urllib.error.URLError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
