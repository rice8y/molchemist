# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.3](https://github.com/Peariforme/bigsmiles-rs/compare/opensmiles-v0.1.2...opensmiles-v0.1.3) - 2026-03-03

### Added

- *(opensmiles)* widen node index from u16 to u32 ([#32](https://github.com/Peariforme/bigsmiles-rs/pull/32))

## [0.1.2](https://github.com/Peariforme/bigsmiles-rs/compare/opensmiles-v0.1.1...opensmiles-v0.1.2) - 2026-02-25

### Fixed

- *(opensmiles)* correct ring-closure index in branch parser ([#26](https://github.com/Peariforme/bigsmiles-rs/pull/26))

## [0.1.1](https://github.com/Peariforme/bigsmiles-rs/compare/opensmiles-v0.1.0...opensmiles-v0.1.1) - 2026-02-25

### Other

- Bigsmiles ([#23](https://github.com/Peariforme/bigsmiles-rs/pull/23))

## [0.1.0](https://github.com/Peariforme/bigsmiles-rs/releases/tag/opensmiles-v0.1.0) - 2026-02-25

### Added

- OpenSMILES compliance audit and fixes ([#17](https://github.com/Peariforme/bigsmiles-rs/pull/17))
- add stereochemistry parsing (chirality + directional bonds) ([#16](https://github.com/Peariforme/bigsmiles-rs/pull/16))
- add parallel batch parsing with rayon and criterion benchmarks ([#7](https://github.com/Peariforme/bigsmiles-rs/pull/7))
- add wildcard atom (*) parsing support ([#4](https://github.com/Peariforme/bigsmiles-rs/pull/4))
- add ring closure support with aromatic bond detection
- add bracket atom parsing with isotope, hydrogen, charge and class
- add explicit bonds and branch parsing support
- initial SMILES parser implementation

### Fixed

- simplify CI by merging feature-powerset testing into single job
- fix comparison charts x-axis labels, align benchmark sizes, remove empty graph ([#14](https://github.com/Peariforme/bigsmiles-rs/pull/14))
- change error messages from French to English

### Other

- Prepare to publish ([#21](https://github.com/Peariforme/bigsmiles-rs/pull/21))
- Display molecule ([#20](https://github.com/Peariforme/bigsmiles-rs/pull/20))
- refine benchmarks ([#12](https://github.com/Peariforme/bigsmiles-rs/pull/12))
- optimize benchmark workflow and add GitHub Pages integration ([#11](https://github.com/Peariforme/bigsmiles-rs/pull/11))
- introduce Parser struct with position tracking
- add comprehensive unit tests for SMILES chemistry features ([#2](https://github.com/Peariforme/bigsmiles-rs/pull/2))
- split ast.rs into modular structure with builder pattern
- make molecule tests comprehensive by checking all node attributes ([#1](https://github.com/Peariforme/bigsmiles-rs/pull/1))
