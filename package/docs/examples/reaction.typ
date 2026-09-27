#import "../../lib.typ": *

#let reaction = "[CH3:1][CH2:2][CH2:3]Br>>[CH3:1][CH2:2][CH2:3]O"
#grid(
  columns: 2, column-gutter: 5mm, row-gutter: 6mm,
  align: (left + horizon, left + horizon),
  [*Highlight off*], render-reaction(reaction, highlight-center: false),
  [*Highlight on*], render-reaction(reaction, highlight-center: true),
)
