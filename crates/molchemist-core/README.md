# molchemist-core

Shared parsing, layout-payload preparation, AST generation, and Alchemist source formatting for molchemist's WebAssembly plugins.

This is an unpublished implementation crate. End users should install `molchemist-cli` and use the `molchemist` executable instead.

The `native-layout` feature compiles the same Coordgen engine for regression testing and development. Production CLI output intentionally executes the packaged WebAssembly modules so it remains byte-for-byte identical to Typst.
