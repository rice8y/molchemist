#import "../../../../package/lib.typ": render-mol, render-smiles

#set page(width: auto, height: auto, margin: 3mm)

#let radical-sdf = "metadata\n  molchemist\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  3  0  0  0  0  0  0  0  7  0  0\nM  ISO  1   1  13\nM  RAD  1   1   2\nM  END\n"

#grid(
  columns: 4,
  gutter: 6mm,
  render-smiles("[13CH3:7]C", abbreviate: true),
  render-smiles("[CH2-]C", skeletal: true),
  render-smiles("*", skeletal: true),
  render-mol(radical-sdf),
)
