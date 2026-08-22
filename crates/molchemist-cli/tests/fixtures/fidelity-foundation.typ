#import "../../../../package/lib.typ": inspect-mol, render-mol

#set page(width: auto, height: auto, margin: 3mm)

#let sdf = "foundation\n  molchemist\nsemantic record\n  0  0  0     0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 1 0 0\nM  V30 BEGIN ATOM\nM  V30 10 C 0.0000 0.0000 0.0000 7 MASS=13\nM  V30 20 O 1.5000 0.0000 0.0000 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 7 1 10 20 RXCTR=1\nM  V30 END BOND\nM  V30 BEGIN SGROUP\nM  V30 3 SUP ATOMS=(1 10) LABEL=\"Me\"\nM  V30 END SGROUP\nM  V30 END CTAB\nM  END\n> <NAME>\nfirst\nline two\n\n> <NAME>\nsecond\n\n$$$$\n"

#let record = inspect-mol(sdf)
#assert.eq(record.schemaVersion, 1)
#assert.eq(record.format, "molfile-v3000")
#assert.eq(record.atoms.first().index, 0)
#assert.eq(record.atoms.first().sourceId, 10)
#assert.eq(record.atoms.first().atomMap, 7)
#assert.eq(record.bonds.first().sourceId, 7)
#assert.eq(record.sgroups.first().kind, "superatom")
#assert.eq(record.properties.map(property => property.name), ("NAME", "NAME"))
#assert.eq(record.properties.first().value, "first\nline two")
#assert.eq(record.diagnostics.map(diagnostic => diagnostic.code), (
  "reaction-center-not-depicted",
))

#render-mol(sdf, fidelity: "ignore")
