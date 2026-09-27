#import "../../../../package/lib.typ": *
#import "../../../../package/src/layout.typ": coordinate-plan, place-callouts, _label-cut
#import "../../../../package/src/labels.typ": _structured-atom-label
#import "../../../../package/src/highlights.typ": _highlight-band, _highlight-label-exclusions
#set page(width: auto, height: auto)

#metadata((
  smiles: render-smiles("[13CH3:7]C(=O)O", skeletal: true, dump: true, config: (layout: "avoid")).text,
)) <coordinate-parity>

#context {
  let atom(i) = (type: "fragment", name: "a" + str(i), element: "Cl", links: ())
  let commands = (atom(0), (type: "bond", name: "b0", bondType: "single", angle: 0, lengthScale: 0.3), atom(1))
  let original = coordinate-plan(commands, 1em, _ => [Chlorine], config: (layout: "coordinates"))
  assert(original.collisions.len() > 0)
  let fitted = coordinate-plan(commands, 1em, _ => [Chlorine], config: (layout: "avoid", collision-policy: "error", max-layout-scale: 64))
  assert(fitted.collisions.len() == 0 and fitted.scale > 1)
  let a = fitted.atoms.at(0)
  let b = fitted.atoms.at(1)
  let padding = measure(h(0.12em)).width / measure(h(1em)).width
  let visible = b.pos.at(0) - a.pos.at(0) - _label-cut(a, (1, 0), padding) - _label-cut(b, (-1, 0), padding)
  assert(visible >= 1, message: "Avoid must leave at least 1em of visible bond")
  let note = _structured-atom-label((symbol: "NO₂"))
  let nitrogen = measure(math.equation(math.upright([N]))).width
  assert(calc.abs(note.offset.at(0) - (measure(note.body).width - nitrogen)/2) < 0.01pt,
    message: "A superatom bond must target its first element, not the label midpoint")
  let carbon = _structured-atom-label((symbol: "C"))
  for left in (false, true) {
    let decorated = _structured-atom-label((symbol: "C", hydrogenCount: 3, isotope: 13, atomMap: 7), left: left)
    assert(decorated.symbol-size == carbon.symbol-size,
      message: "Reaction highlights must exclude hydrogen, isotope and mapping labels")
    assert(measure(decorated.body).width > decorated.symbol-size.at(0))
  }
  let oxygen = _structured-atom-label((symbol: "O"))
  let hydroxyl = _structured-atom-label((symbol: "O", hydrogenCount: 1), left: true)
  assert(hydroxyl.symbol-size == oxygen.symbol-size,
    message: "Only O, not the H in HO, belongs to the atom highlight")
  let contains(polygons, point) = polygons.any(polygon => {
    let sides = ()
    let previous = polygon.last()
    for vertex in polygon {
      sides.push((vertex.at(0)-previous.at(0))*(point.at(1)-previous.at(1))
        - (vertex.at(1)-previous.at(1))*(point.at(0)-previous.at(0)))
      previous = vertex
    }
    sides.all(v => v >= -0.000001) or sides.all(v => v <= 0.000001)
  })
  let exclusions = _highlight-label-exclusions((-2, -1, 2, 1), (0, 0), 1, 0.04,
    ((-2, -0.5, "both"), (0.5, 2, "both")))
  let vertical = _highlight-band((0, 0), (0, 3), 0.7, exclusions: exclusions)
  for i in range(61) {
    assert(contains(vertical, (0, i/20)), message: "The background must reach the symbol center without a gap")
  }
  assert(not contains(vertical, (-0.6, 0.5)) and not contains(vertical, (0.6, 0.5)),
    message: "A wide highlight must not cover the auxiliary label columns")
  assert(contains(vertical, (0.6, 1.2)), message: "Clipping must be limited to the label's height")
  let diagonal = _highlight-band((-3, -1.5), (0, 0), 0.2,
    exclusions: _highlight-label-exclusions((-0.5, -1, 2, 1), (0, 0), 1, 0.04, ((0.5, 2, "both"),)))
  for i in range(61) {
    assert(contains(diagonal, (-3 + i/20, -1.5 + i/40)),
      message: "A diagonal highlight must connect through the foreground label margin")
  }
  let horizontal = _highlight-band((0, 0), (3, 0), 0.2,
    exclusions: _highlight-label-exclusions((-0.5, -1, 2, 1), (0, 0), 1, 0.04, ((0.5, 2, "lower"),)))
  for i in range(61) {
    assert(contains(horizontal, (i/20, 0.1)), message: "A mapping subscript must not sever a horizontal highlight")
  }
  assert(not contains(horizontal, (1, -0.1)), message: "The mapping subscript must remain unhighlighted")
  let fixed = coordinate-plan(commands, 3em, _ => [Cl], config: (
    layout: "coordinates", components: "preserve", atom-positions: (a0: (-2, 1), a1: (2, 1)),
  ))
  assert(fixed.atoms.at(0).pos == (-2, 1) and fixed.atoms.at(1).pos == (2, 1))
  let callouts = place-callouts((callout-annotation(bond-anchor(0), [bond]),), fixed, 3em, (_, body) => body)
  assert(callouts.failures.len() == 0, message: "A bond callout must be able to reach its own target bond")
}
#render-smiles("CC(=O)O", skeletal: true, show-indices: "all", config: (layout: "avoid"), annotations: (
  callout-annotation(bond-anchor(1), [carbonyl bond]),
))
