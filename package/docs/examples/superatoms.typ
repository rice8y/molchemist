#import "../../lib.typ": *

#let data = read("../assets/Sgroups_Abbreviations.mol")
#context {
  let first = render-mol(data, skeletal: true, fidelity: "strict",
    config: (baseline-atom: "a0"))
  let second = render-mol(data, skeletal: true, fidelity: "strict",
    config: (sgroups: "expanded", baseline-atom: "a0"))
  let first-width = measure(first).width
  let second-width = measure(second).width
  let width = calc.max(first-width, second-width)
  box(width: 2 * width + 10mm, {
    grid(columns: (width, width), column-gutter: 10mm, align: center,
      [*Contracted*], [*Expanded*])
    v(3mm)
    h((width - first-width) / 2)
    first
    h(width - (first-width + second-width) / 2 + 10mm)
    second
    h((width - second-width) / 2)
  })
}
