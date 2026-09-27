#import "../../lib.typ": *

#let data = read("../assets/rgroup-alternatives.mol")
#let groups = inspect-rgroup(data)
#let root = render-mol(groups.root, skeletal: true)
#let alternatives = render-rgroup(data)
#grid(
  columns: 2, gutter: 10mm, align: center + top,
  [*Root only* #v(3mm) #root],
  [*Root and alternatives* #v(3mm) #alternatives],
)
