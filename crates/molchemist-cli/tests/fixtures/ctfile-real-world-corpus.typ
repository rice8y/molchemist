#import "../../../../package/lib.typ": inspect-mol, render-mol

#set page(width: 210mm, height: auto, margin: 8mm)
#set text(size: 8pt)

#let rdkit(path) = read("rdkit/" + path, encoding: none)
#let case(title, body) = block(
  width: 100%,
  [
    #text(weight: "bold", size: 8pt)[#title]
    #v(2mm)
    #align(center, body)
  ],
)

#let atom-lists = rdkit("github8823.sdf")
#let second-list = inspect-mol(atom-lists, record: 2)
#assert.eq(second-list.atoms.first().query.elements, ("N", "C", "O"))

#let v2000-sru = inspect-mol(rdkit("Sgroups_SRU_01.mol"))
#assert.eq(v2000-sru.sgroups.first().kind, "structure-repeat-unit")
#assert.eq(v2000-sru.sgroups.first().subscript, "n")

#let data-sgroups = inspect-mol(rdkit("Sgroups_Data_01.mol"))
#assert.eq(data-sgroups.sgroups.map(group => group.fieldName), ("pH", "Stereo"))

#grid(
  columns: 2,
  gutter: 5mm,
  row-gutter: 5mm,
  case(
    "RDKit #8823 — atom list order",
    render-mol(atom-lists, record: 2, fidelity: "strict"),
  ),
  case(
    "Marvin V3000 — HCOUNT",
    render-mol(rdkit("AtomQuery1.mol"), fidelity: "strict"),
  ),
  case(
    "ISIS V2000 — ring topology",
    render-mol(rdkit("RingBondQuery.mol"), fidelity: "strict"),
  ),
  case(
    "ISIS V2000 — chain topology",
    render-mol(rdkit("ChainBondQuery.mol"), fidelity: "strict"),
  ),
  case(
    "ACD/Labs V2000 — SRU",
    render-mol(rdkit("Sgroups_SRU_01.mol"), fidelity: "strict"),
  ),
  case(
    "Marvin V3000 — ranged SRU",
    render-mol(rdkit("repeat_groups_query1.mol"), fidelity: "strict"),
  ),
  case(
    "ACD/Labs V2000 — DAT SGroups",
    render-mol(rdkit("Sgroups_Data_01.mol"), fidelity: "strict"),
  ),
  case(
    "ACD/Labs V2000 — contracted SUP",
    render-mol(rdkit("Sgroups_Abbreviations.mol"), fidelity: "strict"),
  ),
)
