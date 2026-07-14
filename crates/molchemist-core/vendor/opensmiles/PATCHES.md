# Local Patches for molchemist

This directory contains a vendored copy of `opensmiles` used by `molchemist`.
The local copy is intentionally maintained for `molchemist` and is not treated as a temporary patch queue waiting for upstream synchronization.

## Policy

- Keep local changes focused on behavior required by `molchemist`.
- Document each local change here when it affects parsing behavior or public package output.
- Do not block `molchemist` releases on upstreaming these changes.
- Preserve upstream license notices in `THIRD_PARTY_NOTICES.md`.

## Patches

### Branch-leading bond type must not leak into the internal branch chain

- Files:
  - `src/parser.rs`
  - `tests/branches.rs`
- Reason:
  - `molchemist` renders SMILES through a parsed molecular graph. For inputs such as `C(=CC=O)`, the branch-leading `=` describes the bond from the parent atom into the branch. It must not also become the implicit bond inside the branch.
- User-visible bug:
  - CID 30, whose SMILES contains `C(=CC=O)`, could be rendered as if the branch contained `C=C=C=O` rather than `C=C-C=O`, making the carbonyl region look chemically wrong.
- Local change:
  - After assigning a branch-leading bond to `branch_bond_type`, clear `next_bond_type` so it is not reused for the next internal branch bond.
- Regression:
  - `parse_branch_leading_bond_does_not_leak_to_internal_chain` asserts that `C(=CC=O)` parses as double, single, double.
- Verification:
  - From the repository root: `cargo test --manifest-path crates/molchemist-core/vendor/opensmiles/Cargo.toml --test branches parse_branch_leading_bond_does_not_leak_to_internal_chain -- --nocapture`
  - From the repository root: `cargo test -p molchemist-core branch_leading_bond_type_does_not_leak_into_smiles_layout -- --nocapture`

### Keep the spanning-tree helper independent of `Molecule`

- File: `src/ast/molecule.rs`
- Reason: The recursive helper operates entirely on its `SpanningTreeState` argument and does not read the `Molecule` instance.
- Local change: Remove the unused `&self` receiver and call the helper as an associated function. This keeps strict workspace Clippy checks clean without suppressing `clippy::only_used_in_recursion`.
- Behavior: No parsing or canonicalization behavior changes.
- Verification: From the repository root: `cargo clippy --workspace --all-targets -- -D warnings`
