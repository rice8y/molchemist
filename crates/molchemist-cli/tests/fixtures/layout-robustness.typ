#import "../../../../package/lib.typ": render-mol

#set page(width: auto, height: auto, margin: 4mm)

#let collapsed = read("layout-robustness.sdf", encoding: none)

#grid(
  columns: (auto, auto),
  gutter: 5mm,
  align: horizon,
  [Collapsed source coordinates],
  render-mol(collapsed, skeletal: true),
)
