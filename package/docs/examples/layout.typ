#import "../../lib.typ": *

#let molecule = "[13CH3:7]C(=O)O"
#context {
  let first = render-smiles(molecule, abbreviate: true,
    config: (layout: "coordinates", atom-sep: 1.2em))
  let second = render-smiles(molecule, abbreviate: true,
    config: (layout: "avoid", atom-sep: 1.2em))
  let first-width = measure(first).width
  let second-width = measure(second).width
  let width = calc.max(first-width, second-width)
  box(width: 2 * width + 10mm, {
    grid(columns: (width, width), column-gutter: 10mm, align: center,
      [*Coordinates*], [*Avoid*])
    v(3mm)
    h((width - first-width) / 2)
    first
    h(width - (first-width + second-width) / 2 + 10mm)
    second
    h((width - second-width) / 2)
  })
}
