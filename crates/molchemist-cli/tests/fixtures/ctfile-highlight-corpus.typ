#import "../../../../package/lib.typ": inspect-mol, render-mol

#set page(width: 210mm, height: auto, margin: 8mm)
#set text(size: 8pt)

#let shapes = read("highlight-shapes.sdf", encoding: none)
#let real = read("rdkit/v3k.crash1.mol", encoding: none)
#let highlight(hex) = (
  ctfile: (
    highlight-paint: rgb(hex).transparentize(35%),
  ),
)
#let case(title, body) = [
  #text(weight: "bold", size: 8pt)[#title]
  #v(2mm)
  #align(center, body)
]

#let real-record = inspect-mol(real)
#assert.eq(real-record.collections.map(collection => collection.kind), ("highlight", "highlight"))
#assert.eq(real-record.collections.first().bondSourceIds, (7,))
#assert.eq(real-record.collections.last().atomSourceIds, (7,))

#for record in range(1, 6) {
  let inspected = inspect-mol(shapes, record: record)
  assert.eq(inspected.collections.len(), 1)
  assert.eq(inspected.collections.first().kind, "highlight")
  assert.eq(inspected.diagnostics, ())
}

#grid(
  columns: 2,
  gutter: 5mm,
  row-gutter: 5mm,
  case(
    "Real RDKit/Bingo — atom + bond",
    render-mol(real, skeletal: true, fidelity: "strict", config: highlight("#ff6b6b")),
  ),
  case(
    "Atom only — skeletal vertex",
    render-mol(shapes, record: 1, skeletal: true, fidelity: "strict", config: highlight("#74c0fc")),
  ),
  case(
    "Atom only — long query glyph",
    render-mol(shapes, record: 2, skeletal: true, fidelity: "strict", config: highlight("#51cf66")),
  ),
  case(
    "Bond only — diagonal capsule",
    render-mol(shapes, record: 3, skeletal: true, fidelity: "strict", config: highlight("#ffd43b")),
  ),
  case(
    "Connected atom + bond union",
    render-mol(shapes, record: 4, skeletal: true, fidelity: "strict", config: highlight("#cc5de8")),
  ),
  case(
    "Disconnected highlight regions",
    render-mol(shapes, record: 5, skeletal: true, fidelity: "strict", config: highlight("#ff922b")),
  ),
)
