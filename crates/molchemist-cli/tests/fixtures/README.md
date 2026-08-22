# CLI Test Structures

These SDF files are copies of the documentation fixtures used to test file input and multi-record selection:

- `Structure2D_COMPOUND_CID_241.sdf`: PubChem Compound CID 241
- `Structure2D_COMPOUND_CID_93406.sdf`: PubChem Compound CID 93406
- `bond-semantics.sdf`: synthetic V3000 coverage for aromatic, query, coordination, hydrogen, and either-stereo bonds
- `stereochemistry.sdf`: synthetic V3000 coverage for atom parity, enhanced stereo groups, stereochemical hydrogen bonds, and undefined double bonds
- `layout-robustness.sdf`: synthetic V2000 coverage for collapsed XY and 3D-only coordinates while retaining source wedge semantics
- `ctfile-fidelity.sdf`: synthetic V3000 coverage for NOT atom lists, atom query constraints, R-groups, bond topology, SGroup brackets, ordered SDF properties, and atom/bond highlights
- `ctfile-attachments-collections.sdf`: synthetic V3000 coverage for position-variation bond endpoints/modes, standard atom-list and R-group atom attributes, and preservation of a user-defined COLLECTION
- `highlight-shapes.sdf`: synthetic V3000 boundary matrix for atom-only, bond-only, connected, disconnected, skeletal-carbon, and long-query-glyph highlights; `ctfile-highlight-corpus.typ` renders the matrix as an SVG shape regression test
- `rdkit/github8823.sdf`: real RDKit regression corpus with three V3000 atom lists whose element ordering differs
- `rdkit/AtomQuery1.mol`: real Marvin V3000 query structure with raw `HCOUNT=3`, normalized to the MDL `H3` query
- `rdkit/RingBondQuery.mol` and `rdkit/ChainBondQuery.mol`: real ISIS V2000 ring/chain bond-topology queries
- `rdkit/Sgroups_SRU_01.mol`: real ACD/Labs V2000 SRU with explicit brackets, crossing bonds, head-to-tail connectivity, and an `n` subscript
- `rdkit/repeat_groups_query1.mol`: real Marvin V3000 SRU with a `1-3` repeat range
- `rdkit/Sgroups_Data_01.mol`: real ACD/Labs V2000 DAT SGroups with `pH` and `Stereo` fields
- `rdkit/Sgroups_Abbreviations.mol`: real ACD/Labs V2000 SUP SGroups for the `NO2` and `COOH` abbreviations; this exercises multi-atom contraction, crossing-bond reconnection, and strict-fidelity rendering
- `rdkit/Sgroups_Link_01.mol`: real ACD/Labs V2000 link-node structure with a `1–3` repeat range and two substituent connections
- `rdkit/sgroup_ap_bug.mol`: real Marvin V2000 SAP regression structure with single- and multi-attachment SUP SGroups and named `Al`/`Br` connections
- `rdkit/v3k.crash1.mol`: real Bingo V3000 regression fixture with separate `MDLV30/HILITE` atom and bond collections

The files under `rdkit/` are copied from RDKit commit `b421f19c9f564d0cb66148c4e614c59abadf5413` (2026-08-21). Their original paths are `Code/GraphMol/FileParsers/test_data/` and `Code/GraphMol/FileParsers/sgroup_test_data/`; only CRLF line endings were normalized to LF. `github8823.sdf` is the committed regression data associated with RDKit issues 8820/8823 and pull request 8824.

Source and data-usage details are recorded in `THIRD_PARTY_NOTICES.md`.

## Highlight regression checks

Run the semantic membership checks with:

```sh
cargo test -p molchemist-cli --test cli highlight_
```

When Typst is installed, run the SVG shape regression (including the real RDKit/Bingo case) with:

```sh
cargo test -p molchemist-cli --test cli local_typst_package_renders_highlight_shape_corpus_to_svg -- --ignored --exact
```

For manual visual inspection, compile the same six-case corpus directly:

```sh
typst compile --root . crates/molchemist-cli/tests/fixtures/ctfile-highlight-corpus.typ /tmp/ctfile-highlight-corpus.svg
```
