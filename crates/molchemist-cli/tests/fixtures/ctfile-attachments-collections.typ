#import "../../../../package/lib.typ": *

#set page(width: 150mm, height: auto, margin: 8mm)

#let data = read("ctfile-attachments-collections.sdf", encoding: none)
#let real-link = read("rdkit/Sgroups_Link_01.mol", encoding: none)

#grid(
  columns: 2,
  gutter: 10mm,
  align: center,
  render-mol(data, record: 1, skeletal: true, fidelity: "strict"),
  render-mol(data, record: 2, skeletal: true, fidelity: "strict"),
)
#v(5mm)
#align(center)[
  #render-mol(real-link, skeletal: true, fidelity: "strict")
]
