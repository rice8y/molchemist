//! Benchmarks for SMILES parsing performance.
//!
//! Run with:
//!   cargo bench -p opensmiles                     # without parallel
//!   cargo bench -p opensmiles --features parallel  # with parallel
//!
//! Benchmark groups:
//! - `reference`: 5 representative molecules for cross-language comparison
//! - `seq_vs_parallel`: find the batch size where parallelism pays off (parallel feature only)
//! - `scaling`: parser performance vs molecule size, tracks memory footprint
//! - `huckel`: overhead of Hückel aromaticity validation on aromatic molecules

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use opensmiles::ast::aromaticity::validate_aromaticity;
use opensmiles::{parse, Molecule};

#[cfg(feature = "parallel")]
use opensmiles::parser_parallel::parse_batch;

#[cfg(feature = "parallel")]
/// Diverse set of molecules for batch benchmarks.
const BATCH_MOLECULES: &[&str] = &[
    "C",                            // methane
    "CCO",                          // ethanol
    "C=O",                          // formaldehyde
    "C1CCCCC1",                     // cyclohexane
    "c1ccccc1",                     // benzene
    "CC(C)C",                       // isobutane
    "CC(=O)O",                      // acetic acid
    "c1ccc2ccccc2c1",               // naphthalene
    "CC(C)Cc1ccc(C(C)C(=O)O)cc1",   // ibuprofen-like
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", // caffeine
];

/// PEG (polyethylene glycol): OCCOCCOCC...O — linear chain with O and C
fn generate_peg(n: usize) -> String {
    "OCC".repeat(n) + "O"
}

/// Teflon (PTFE): C(F)(F)C(F)(F)... — branching stress test with fluorine
fn generate_teflon(n: usize) -> String {
    "C(F)(F)".repeat(n)
}

#[cfg(feature = "parallel")]
fn generate_batch_dataset(size: usize) -> Vec<&'static str> {
    (0..size)
        .map(|i| BATCH_MOLECULES[i % BATCH_MOLECULES.len()])
        .collect()
}

/// Estimate heap memory used by a Molecule (nodes + bonds vectors).
fn estimate_heap_bytes(mol: &Molecule) -> usize {
    std::mem::size_of_val(mol.nodes()) + std::mem::size_of_val(mol.bonds())
}

/// Reference benchmarks: 5 representative molecules for cross-language comparison.
///
/// Covers: simple chain, ring closure, aromaticity, branching, complex real molecule.
/// Use these to compare against other SMILES parsers (RDKit, OpenBabel, etc.).
fn bench_reference(c: &mut Criterion) {
    let mut group = c.benchmark_group("reference");

    let molecules = [
        ("ethanol", "CCO"),
        ("cyclohexane", "C1CCCCC1"),
        ("benzene", "c1ccccc1"),
        ("ibuprofen", "CC(C)Cc1ccc(C(C)C(=O)O)cc1"),
        ("caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"),
    ];

    for (name, smiles) in &molecules {
        group.bench_with_input(BenchmarkId::from_parameter(name), smiles, |b, s| {
            b.iter(|| parse(black_box(s)))
        });
    }

    group.finish();
}

/// Sequential vs parallel crossover: find the batch size where parallelism pays off.
///
/// Tests 5 batch sizes from 10 to 10000. Look for the crossover point
/// where parallel throughput (elements/sec) exceeds sequential.
#[cfg(feature = "parallel")]
fn bench_seq_vs_parallel(c: &mut Criterion) {
    let mut group = c.benchmark_group("seq_vs_parallel");

    for size in [10, 100, 1000, 5000, 10000, 50000] {
        let dataset = generate_batch_dataset(size);
        group.throughput(Throughput::Elements(size as u64));

        group.bench_with_input(BenchmarkId::new("sequential", size), &dataset, |b, data| {
            b.iter(|| {
                let results: Vec<_> = data.iter().map(|s| parse(black_box(s))).collect();
                black_box(results)
            })
        });

        group.bench_with_input(BenchmarkId::new("parallel", size), &dataset, |b, data| {
            b.iter(|| black_box(parse_batch(black_box(data))))
        });
    }

    group.finish();
}

/// Parser scaling: performance vs molecule size, with memory footprint tracking.
///
/// Uses Throughput::Bytes so Criterion reports bytes/sec — if structs grow bigger
/// after a refactor, the throughput numbers will shift and Criterion will flag it.
fn bench_scaling(c: &mut Criterion) {
    let mut group = c.benchmark_group("scaling");

    // PEG chains: detect O(n^2) algorithmic regressions (n = repeat units)
    for n in [100, 500, 1000, 5000] {
        let smiles = generate_peg(n);
        let mol = parse(&smiles).expect("should parse");
        let heap_bytes = estimate_heap_bytes(&mol);

        group.throughput(Throughput::Bytes(heap_bytes as u64));
        group.bench_with_input(BenchmarkId::new("peg", n), &smiles, |b, s| {
            b.iter(|| parse(black_box(s)))
        });
    }

    // Teflon (PTFE): 2 branches per carbon, tests branch parsing with fluorine
    for n in [100, 500, 1000, 5000] {
        let smiles = generate_teflon(n);
        let mol = parse(&smiles).expect("should parse");
        let heap_bytes = estimate_heap_bytes(&mol);

        group.throughput(Throughput::Bytes(heap_bytes as u64));
        group.bench_with_input(BenchmarkId::new("teflon", n), &smiles, |b, s| {
            b.iter(|| parse(black_box(s)))
        });
    }

    group.finish();
}

/// Hückel validation overhead: measures the cost of aromaticity checking
/// after parsing, for molecules with different aromatic complexity.
///
/// Benchmarks `parse()` alone vs `parse() + validate_aromaticity()` to
/// isolate the validation overhead. Useful for deciding whether to enable
/// the `huckel-validation` feature flag.
fn bench_huckel(c: &mut Criterion) {
    let mut group = c.benchmark_group("huckel");

    let molecules = [
        ("ethanol", "CCO"),                           // no aromatic rings (baseline)
        ("benzene", "c1ccccc1"),                      // 1 ring, 6 atoms
        ("naphthalene", "c1ccc2ccccc2c1"),            // 2 fused rings
        ("ibuprofen", "CC(C)Cc1ccc(C(C)C(=O)O)cc1"),  // 1 ring, complex
        ("caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"), // Kekulé form
    ];

    for (name, smiles) in &molecules {
        group.bench_with_input(BenchmarkId::new("parse_only", name), smiles, |b, s| {
            b.iter(|| parse(black_box(s)))
        });

        group.bench_with_input(
            BenchmarkId::new("parse_and_validate", name),
            smiles,
            |b, s| {
                b.iter(|| {
                    let mol = parse(black_box(s)).unwrap();
                    validate_aromaticity(black_box(&mol));
                    mol
                })
            },
        );
    }

    group.finish();
}

#[cfg(feature = "parallel")]
criterion_group!(
    benches,
    bench_reference,
    bench_seq_vs_parallel,
    bench_scaling,
    bench_huckel,
);

#[cfg(not(feature = "parallel"))]
criterion_group!(benches, bench_reference, bench_scaling, bench_huckel,);

criterion_main!(benches);
