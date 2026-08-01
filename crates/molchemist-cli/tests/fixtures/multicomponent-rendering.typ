#import "../../../../package/lib.typ": *

#set page(width: auto, height: auto, margin: 4mm)

#let salt-sdf = "sodium chloride\n  molchemist\n\n  2  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 Na  0  0  0  0  0  0  0  0  0  0  0  0\n    3.0000    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\nM  CHG  2   1   1   2  -1\nM  END\n"

#grid(
  columns: 1,
  row-gutter: 7mm,
  [SDF ions: #render-mol(salt-sdf, abbreviate: true)],
  [SMILES ions: #render-smiles(
    "[Na+].[Cl-]",
    abbreviate: true,
    show-indices: true,
    annotations: callout-annotation(atom-anchor(1), [chloride], side: "south"),
  )],
  [Visible isolated components: #render-smiles("[H+].C.[Cl-]", skeletal: true)],
)
