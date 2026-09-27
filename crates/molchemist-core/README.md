# molchemist-core

Shared parsing, layout-payload preparation, AST generation, and Alchemist source formatting for molchemist's WebAssembly plugins.

The crate also defines the versioned `ChemicalRecord` semantic IR used before SDF depiction. `inspect_sdf_record` preserves source IDs, headers, normalized V2000/V3000 atom queries and SGroups, original SGroup type codes, SAP attachment points, variable-bond endpoints/modes, link nodes, full COLLECTION names and member lists, ordered/duplicate SDF properties, the exact selected record, and diagnostics for remaining unsupported features. SDF rendering and CTfile overlays are generated from that same semantic record, including query labels, R-groups, bond topology, SGroup brackets, variable bonds, link nodes, and V3000 atom/bond highlights. User-defined and unknown internal collections are preserved without inventing a depiction; strict fidelity reports them. Highlight formatting keeps atom and bond membership distinct, resolves all skeletal vertices through atom-center anchors, uses rendered glyph anchors for visible labels, and uses consistently wound non-zero compound paths for connected unions.

This is an unpublished implementation crate. End users should install `molchemist-cli` and use the `molchemist` executable instead.

The `native-layout` feature compiles the same Coordgen engine for regression testing and development. Production CLI output intentionally executes the packaged WebAssembly modules so it remains byte-for-byte identical to Typst.

`format_coordinate_code` serializes the AST and configuration, then embeds the ordered drawing modules under `package/src`. These are the same modules used by the package entry point; the Rust formatter contains no separate renderer. Register drawing modules in the formatter’s `RENDER_MODULES` list in dependency order and keep their imports on one line. Local imports are resolved within the bundle; the generated source retains one Alchemist import. Rebuild the plugin after changing shared rendering code; see [the development scripts](../../scripts/README.md).
