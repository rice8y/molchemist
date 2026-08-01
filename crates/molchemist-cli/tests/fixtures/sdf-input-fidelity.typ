#import "../../../../package/lib.typ": render-mol

#set page(width: auto, height: auto, margin: 3mm)

#let v2000 = "carbon\n  molchemist\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n"
#let v3000 = "charged\n  molchemist\n\n  0  0  0     0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 N 0.0000 0.0000 0.0000 7 CHG=1 MASS=15 RAD=2\nM  V30 2 O 1.5000 0.0000 0.0000 0 CHG=-1\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2 CFG=1\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n"
#let records = v2000 + "$$$$\n" + v3000 + "$$$$\n"

#grid(
  columns: 2,
  gutter: 6mm,
  render-mol(records, record: 1),
  render-mol(records, record: 2),
)
