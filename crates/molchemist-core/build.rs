use std::env;
use std::path::Path;

const COORDGEN_SOURCES: &[&str] = &[
    "CoordgenFragmentBuilder.cpp",
    "CoordgenFragmenter.cpp",
    "CoordgenMacrocycleBuilder.cpp",
    "CoordgenMinimizer.cpp",
    "CoordgenTemplates.cpp",
    "sketcherMinimizer.cpp",
    "sketcherMinimizerAtom.cpp",
    "sketcherMinimizerBond.cpp",
    "sketcherMinimizerFragment.cpp",
    "sketcherMinimizerMarchingSquares.cpp",
    "sketcherMinimizerMolecule.cpp",
    "sketcherMinimizerResidue.cpp",
    "sketcherMinimizerResidueInteraction.cpp",
    "sketcherMinimizerRing.cpp",
];

fn main() {
    if env::var_os("CARGO_FEATURE_NATIVE_LAYOUT").is_none() {
        return;
    }

    let vendor = Path::new("vendor/coordgenlibs");
    let mut build = cc::Build::new();
    build
        .cpp(true)
        .std("c++17")
        .define("STATIC_COORDGEN", None)
        .include(vendor)
        .include("native")
        .file("native/coordgen_engine.cpp")
        .file("native/native_bridge.cpp")
        .warnings(false);

    for source in COORDGEN_SOURCES {
        build.file(vendor.join(source));
    }

    build.compile("molchemist_coordgen");

    println!("cargo:rerun-if-changed=native/coordgen_engine.cpp");
    println!("cargo:rerun-if-changed=native/coordgen_engine.h");
    println!("cargo:rerun-if-changed=native/native_bridge.cpp");
    println!("cargo:rerun-if-changed=vendor/coordgenlibs");
}
