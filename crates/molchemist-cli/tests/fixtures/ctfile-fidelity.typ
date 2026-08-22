#import "../../../../package/lib.typ": inspect-mol, render-mol

#set page(width: auto, height: auto, margin: 3mm)

#let sdf = read("ctfile-fidelity.sdf")

#let record = inspect-mol(sdf)
#assert.eq(record.atoms.first().query.elements, ("C", "N"))
#assert.eq(record.atoms.first().query.isNotList, true)
#assert.eq(record.atoms.first().query.hydrogenCount, 1)
#assert.eq(record.atoms.first().query.substitutionCount, 3)
#assert.eq(record.atoms.at(1).rgroupLabels, (7,))
#assert.eq(record.sgroups.first().brackets.len(), 2)
#assert.eq(record.collections.first().atomSourceIds, (10, 20))
#assert.eq(record.collections.first().bondSourceIds, (5,))
#assert.eq(record.properties.map(property => property.name), ("SOURCE", "SOURCE"))
#assert.eq(record.properties.first().value, "first\nline two")
#assert.eq(record.diagnostics, ())

#render-mol(
  sdf,
  fidelity: "strict",
  config: (
    ctfile: (
      highlight-paint: rgb("#74c0fc").transparentize(45%),
    ),
  ),
)
