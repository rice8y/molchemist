#import "../../../../package/lib.typ": render-mol, render-smiles

#let dump-smiles(smiles, mode) = render-smiles(smiles, skeletal: mode == "skeletal", abbreviate: mode == "abbreviate", dump: true).text
#let dump-sdf(sdf, mode, record: "1") = render-mol(sdf, record: int(record), skeletal: mode == "skeletal", abbreviate: mode == "abbreviate", dump: true).text

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
  ctfile-fidelity: dump-sdf(
    read("ctfile-fidelity.sdf", encoding: none),
    "full",
  ),
  benzene: dump-smiles("c1ccccc1", "skeletal"),
  charged: dump-smiles("OCCc1c(C)[n+](=cs1)Cc2cnc(C)nc(N)2", "abbreviate"),
  chiral: dump-smiles("N[C@@H](C)C(=O)O", "full"),
  ez: dump-smiles("F/C=C\\F", "skeletal"),
  isotope-map: dump-smiles("[13CH3:7]C", "abbreviate"),
  charged-carbon: dump-smiles("[CH2-]C", "skeletal"),
  wildcard: dump-smiles("*", "skeletal"),
  quadruple: dump-smiles("[Cr]$[Cr]", "skeletal"),
  multicomponent: dump-smiles("[Na+].[Cl-]", "abbreviate"),
  component-baselines: dump-smiles("[H+].C.[Cl-]", "skeletal"),
  complex: dump-smiles("CC[C@@H]([C@@H]1[C@H](C[C@@](O1)(CC)[C@H]2CC[C@@]([C@@H](O2)C)(CC)O)C)C(=O)[C@@H](C)[C@H]([C@H](C)CCC3=C(C=C(C(=C3C(=O)O)O)C)Br)O", "abbreviate"),
)) <parity>
