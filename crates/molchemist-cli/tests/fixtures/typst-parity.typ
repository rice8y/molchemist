#let mol-plugin = plugin("../../wasm/molchemist_plugin.wasm")
#let smiles-plugin = plugin("../../wasm/molchemist_smiles_plugin.wasm")

#let dump-smiles(smiles, mode) = {
  let layout-input = if mode == "full" {
    mol-plugin.smiles_to_full_layout_input(bytes(smiles))
  } else {
    mol-plugin.smiles_to_layout_input(bytes(smiles))
  }
  let coords = smiles-plugin.layout_coordinates(layout-input)
  str(mol-plugin.smiles_to_code(
    bytes(smiles),
    coords,
    bytes(mode),
    bytes("3em"),
    bytes("2"),
  ))
}

#let dump-sdf(sdf, mode, record: "1") = {
  let record = bytes(record)
  let layout-input = mol-plugin.sdf_record_to_layout_input(sdf, record)
  if layout-input.len() == 0 {
    return str(mol-plugin.sdf_record_to_code(
      sdf,
      bytes(mode),
      record,
      bytes("3em"),
      bytes("2"),
    ))
  }
  let coords = smiles-plugin.layout_coordinates(layout-input)
  str(mol-plugin.sdf_record_to_code_with_coords(
    sdf,
    coords,
    bytes(mode),
    record,
    bytes("3em"),
    bytes("2"),
  ))
}

#metadata((
  sdf: dump-sdf(
    read("Structure2D_COMPOUND_CID_241.sdf", encoding: none),
    "abbreviate",
  ),
  bond-semantics: dump-sdf(
    read("bond-semantics.sdf", encoding: none),
    "full",
  ),
  stereochemistry: dump-sdf(
    read("stereochemistry.sdf", encoding: none),
    "skeletal",
  ),
  collapsed-sdf: dump-sdf(
    read("layout-robustness.sdf", encoding: none),
    "skeletal",
  ),
  benzene: dump-smiles("c1ccccc1", "skeletal"),
  charged: dump-smiles("OCCc1c(C)[n+](=cs1)Cc2cnc(C)nc(N)2", "abbreviate"),
  chiral: dump-smiles("N[C@@H](C)C(=O)O", "full"),
  ez: dump-smiles("F/C=C\\F", "skeletal"),
  isotope-map: dump-smiles("[13CH3:7]C", "abbreviate"),
  charged-carbon: dump-smiles("[CH2-]C", "skeletal"),
  wildcard: dump-smiles("*", "skeletal"),
  multicomponent: dump-smiles("[Na+].[Cl-]", "abbreviate"),
  complex: dump-smiles("CC[C@@H]([C@@H]1[C@H](C[C@@](O1)(CC)[C@H]2CC[C@@]([C@@H](O2)C)(CC)O)C)C(=O)[C@@H](C)[C@H]([C@H](C)CCC3=C(C=C(C(=C3C(=O)O)O)C)Br)O", "abbreviate"),
)) <parity>
