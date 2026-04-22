# opensmiles

A fast, correct SMILES parser for Rust, following the [OpenSMILES specification](http://opensmiles.org/opensmiles.html).

[![Crates.io](https://img.shields.io/crates/v/opensmiles.svg)](https://crates.io/crates/opensmiles)
[![docs.rs](https://docs.rs/opensmiles/badge.svg)](https://docs.rs/opensmiles)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![CI](https://github.com/Peariforme/bigsmiles-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/Peariforme/bigsmiles-rs/actions/workflows/ci.yml)
[![Benchmarks](https://img.shields.io/badge/benchmarks-view-blue)](https://peariforme.github.io/bigsmiles-rs/dev/bench/)

## Features

- Full [OpenSMILES](http://opensmiles.org/opensmiles.html) compliance — all 118 elements, organic subset, bracket atoms, rings, branches, stereochemistry
- Canonical SMILES output via `Display` (round-trip)
- Detailed parse errors with character position
- Optional parallel batch parsing with Rayon
- Optional Hückel's rule aromaticity validation (4n+2 π-electron check)
- Zero unsafe code, no C dependencies

## Installation

```toml
[dependencies]
opensmiles = "0.1"
```

With optional features:

```toml
[dependencies]
opensmiles = { version = "0.1", features = ["parallel", "huckel-validation"] }
```

## Usage

### Basic parsing

```rust
use opensmiles::parse;

let mol = parse("CCO").unwrap();            // ethanol
let mol = parse("c1ccccc1").unwrap();       // benzene
let mol = parse("[C@H](F)(Cl)Br").unwrap(); // chiral center

// Access atoms
for node in mol.nodes() {
    println!(
        "{} — aromatic: {}, H: {}, charge: {}",
        node.atom().element(),
        node.aromatic(),
        node.hydrogens(),
        node.atom().charge(),
    );
}

// Access bonds
for bond in mol.bonds() {
    println!("{} → {} ({:?})", bond.source(), bond.target(), bond.kind());
}
```

### Round-trip to canonical SMILES

`Molecule` implements `Display`, which serializes back to a canonical SMILES string:

```rust
use opensmiles::parse;

let mol = parse("OCC").unwrap();
println!("{}", mol); // OCC
```

### Error handling

```rust
use opensmiles::{parse, ParserError};

match parse("C(C") {
    Ok(mol) => println!("parsed: {}", mol),
    Err(ParserError::UnclosedParenthesis) => eprintln!("missing closing )"),
    Err(e) => eprintln!("parse error: {}", e),
}
```

### Parallel batch parsing

Enable the `parallel` feature for multi-threaded parsing of large datasets:

```rust
use opensmiles::parse_batch;

let dataset = vec!["CCO", "c1ccccc1", "CC(=O)O", /* ... */];
let results = parse_batch(&dataset); // Vec<Result<Molecule, ParserError>>
```

**Benchmark results (4-core CPU):**

| Batch size | Sequential | Parallel | Speedup |
|------------|-----------|----------|---------|
| 100        | 76 µs     | 169 µs   | 0.45×   |
| 1 000      | 877 µs    | 396 µs   | **2.2×**|
| 10 000     | 8.6 ms    | 2.2 ms   | **3.9×**|

> For batches smaller than ~500 molecules, sequential is faster due to thread overhead.

See the full [benchmark dashboard](https://peariforme.github.io/bigsmiles-rs/dev/bench/) and [sequential vs parallel comparison](https://peariforme.github.io/bigsmiles-rs/dev/bench/compare.html).

### Aromaticity validation (Hückel's rule)

Enable the `huckel-validation` feature to have `parse()` reject chemically invalid aromatic rings:

```rust
// With the feature enabled, this returns Err(MoleculeError::HuckelViolation)
// for rings that don't satisfy 4n+2 π electrons.
let mol = parse("c1ccccc1").unwrap(); // benzene: 6 π-electrons ✓
```

The validation API is also available explicitly, without the feature flag:

```rust
use opensmiles::{parse, ast::aromaticity::validate_aromaticity};

let mol = parse("c1ccccc1").unwrap();
let checks = validate_aromaticity(&mol);
assert!(checks[0].is_valid);
assert_eq!(checks[0].pi_electrons, Some(6));
```

## Supported SMILES features

| Feature | Status |
|---------|--------|
| All 118 elements in bracket atoms | ✅ |
| Organic subset (B C N O P S F Cl Br I) with implicit H | ✅ |
| Wildcard `*` | ✅ |
| Isotopes `[13C]` | ✅ |
| Formal charges `[NH4+]`, `[Fe-3]` | ✅ |
| Explicit hydrogen count `[CH3]` | ✅ |
| Atom class `[C:1]` | ✅ |
| Single, double, triple, quadruple bonds | ✅ |
| Aromatic bonds `:` | ✅ |
| Directional bonds `/` `\` (E/Z stereochemistry) | ✅ |
| Branches `()` with arbitrary nesting | ✅ |
| Ring closures 0–9 and `%10`–`%99` | ✅ |
| Disconnected structures `.` | ✅ |
| Tetrahedral chirality `@` `@@` | ✅ |
| Extended chirality `@TH`, `@AL`, `@SP`, `@TB`, `@OH` | ✅ |
| Kekule aromatic forms | ✅ |
| Aromatic bracket symbols `[se]`, `[as]` | ✅ |
| Whitespace terminator | ✅ |

## Feature flags

| Flag | Default | Description |
|------|---------|-------------|
| `parallel` | off | Multi-threaded batch parsing via [Rayon](https://crates.io/crates/rayon) |
| `huckel-validation` | off | Reject aromatic rings violating Hückel's 4n+2 rule in `parse()` |

## Part of the bigsmiles-rs ecosystem

`opensmiles` is the SMILES foundation of [bigsmiles-rs](https://github.com/Peariforme/bigsmiles-rs).
The [`bigsmiles`](https://crates.io/crates/bigsmiles) crate extends it with support for polymer notation.

## References

- [OpenSMILES Specification](http://opensmiles.org/opensmiles.html)
- [SMILES Formal Grammar (LL(1))](https://depth-first.com/articles/2020/12/21/smiles-formal-grammar-revisited/)
- Weininger, D. *J. Chem. Inf. Comput. Sci.* **1988**, 28, 31–36.

## License

MIT — see [LICENSE](https://github.com/Peariforme/bigsmiles-rs/blob/master/LICENSE).
