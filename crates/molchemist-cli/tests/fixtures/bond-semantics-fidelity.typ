#import "../../../../package/lib.typ": render-mol

#set page(width: auto, height: auto, margin: 4mm)

#let v2000(order, stereo: 0) = {
  let sdf = "bond semantics\n  molchemist\n\n"
  sdf += "  2  1  0  0  0  0  0  0  0  0999 V2000\n"
  sdf += "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
  sdf += "    1.5000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n"
  sdf += "  1  2  " + str(order) + "  " + str(stereo) + "  0  0  0\n"
  sdf += "M  END\n"
  sdf
}

#let v3000(order, atom1: 1, atom2: 2) = {
  let sdf = "bond semantics\n  molchemist\n\n"
  sdf += "  0  0  0     0  0            999 V3000\n"
  sdf += "M  V30 BEGIN CTAB\n"
  sdf += "M  V30 COUNTS 2 1 0 0 0\n"
  sdf += "M  V30 BEGIN ATOM\n"
  sdf += "M  V30 1 C 0.0000 0.0000 0.0000 0\n"
  sdf += "M  V30 2 N 1.5000 0.0000 0.0000 0\n"
  sdf += "M  V30 END ATOM\n"
  sdf += "M  V30 BEGIN BOND\n"
  sdf += "M  V30 1 " + str(order) + " " + str(atom1) + " " + str(atom2) + "\n"
  sdf += "M  V30 END BOND\n"
  sdf += "M  V30 END CTAB\n"
  sdf += "M  END\n"
  sdf
}

#grid(
  columns: (9.5em, auto),
  row-gutter: 2mm,
  column-gutter: 5mm,
  align: (left + horizon, center + horizon),
  [Aromatic], render-mol(v2000(4)),
  [Single or double], render-mol(v2000(5)),
  [Single or aromatic], render-mol(v2000(6)),
  [Double or aromatic], render-mol(v2000(7)),
  [Any], render-mol(v2000(8)),
  [Either stereo], render-mol(v2000(1, stereo: 4)),
  [Coordination →], render-mol(v3000(9)),
  [Coordination ←], render-mol(v3000(9, atom1: 2, atom2: 1)),
  [Hydrogen], render-mol(v3000(10)),
)
