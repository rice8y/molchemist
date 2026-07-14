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

#metadata((
  sdf: str(mol-plugin.sdf_to_code(
    read("Structure2D_COMPOUND_CID_241.sdf", encoding: none),
    bytes("abbreviate"),
    bytes("3em"),
    bytes("2"),
  )),
  benzene: dump-smiles("c1ccccc1", "skeletal"),
  charged: dump-smiles("OCCc1c(C)[n+](=cs1)Cc2cnc(C)nc(N)2", "abbreviate"),
  chiral: dump-smiles("N[C@@H](C)C(=O)O", "full"),
  ez: dump-smiles("F/C=C\\F", "skeletal"),
  complex: dump-smiles("CC[C@@H]([C@@H]1[C@H](C[C@@](O1)(CC)[C@H]2CC[C@@]([C@@H](O2)C)(CC)O)C)C(=O)[C@@H](C)[C@H]([C@H](C)CCC3=C(C=C(C(=C3C(=O)O)O)C)Br)O", "abbreviate"),
)) <parity>
