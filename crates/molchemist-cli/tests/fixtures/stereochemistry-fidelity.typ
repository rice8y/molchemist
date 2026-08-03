#import "../../../../package/lib.typ": render-mol, render-smiles

#set page(width: auto, height: auto, margin: 4mm)

#let sdf = read("stereochemistry.sdf", encoding: none)

#grid(
  columns: (8em, auto),
  row-gutter: 3mm,
  column-gutter: 5mm,
  align: (left + horizon, center + horizon),
  [SDF stereo], render-mol(sdf, skeletal: true),
  [D / L alanine], grid(
    columns: 2,
    gutter: 4mm,
    render-smiles("N[C@H](C)C(=O)O", abbreviate: true),
    render-smiles("N[C@@H](C)C(=O)O", abbreviate: true),
  ),
  [Implicit / explicit H], grid(
    columns: 2,
    gutter: 4mm,
    render-smiles("N[C@@H](C)C(=O)O", skeletal: true),
    render-smiles("N[C@@]([H])(C)C(=O)O", skeletal: true),
  ),
  [Leading center], grid(
    columns: 2,
    gutter: 4mm,
    render-smiles("[C@H](F)(Cl)Br", skeletal: true),
    render-smiles("[C@@H](F)(Cl)Br", skeletal: true),
  ),
  [Ring order], render-smiles("[C@]1(Br)(Cl)CCCC(F)C1", skeletal: true),
  [Double bond], grid(
    columns: 2,
    gutter: 4mm,
    render-smiles("F/C=C/F", skeletal: true),
    render-smiles("F/C=C\\F", skeletal: true),
  ),
  [Extended], grid(
    columns: 2,
    gutter: 4mm,
    render-smiles("NC(Br)=[C@AL1]=C(O)C", skeletal: true),
    render-smiles("[Pt@SP2](F)(Cl)(Br)I", skeletal: true),
    render-smiles("[As@TB5](F)(Cl)(Br)(N)S", skeletal: true),
    render-smiles("[Co@OH5](F)(Cl)(Br)(I)(N)S", skeletal: true),
  ),
)
