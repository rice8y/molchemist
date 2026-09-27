#import "../../../../package/lib.typ": *
#set page(width: auto, height: auto, margin: 5mm)
#set text(font: "Libertinus Serif", size: 10pt)

#render-smiles("CC(=O)O", config: (layout: "avoid", collision-policy: "error"), annotations: (
  callout-annotation(atom-anchor(2), [Carbonyl]),
  callout-annotation(atom-anchor(3), [Hydroxyl]),
))
#pagebreak()
#render-smiles("[13CH3:7][C@@H](O)C(=O)O.[Na+]", skeletal: true, config: (layout: "avoid", collision-policy: "error"))
#pagebreak()
#render-smiles("[Pt@SP1]1(F)(Cl)CCC1", skeletal: true, config: (layout: "avoid"))
#pagebreak()
#render-mol2(read("extended/water.mol2"), config: (layout: "avoid", collision-policy: "error"))
#pagebreak()
#render-mol(read("extended/tetra-3d.mol"), skeletal: true, config: (layout: "reflow", infer-stereo: true))
#pagebreak()
#render-reaction("[CH3:1][CH2:2][CH2:3]Br>>[CH3:1][CH2:2][CH2:3]O", conditions: [substitution])
#pagebreak()
#render-reaction(read("extended/oxidation.rxn"), format: "rxn", conditions: [oxidation])
#pagebreak()
#render-rgroup(read("extended/alternatives.mol"))
#pagebreak()
#render-mol(read("rdkit/Sgroups_Abbreviations.mol"), skeletal: true, fidelity: "strict", config: (layout: "avoid", sgroups: "expanded"))
#pagebreak()
#grid(columns: 2, gutter: 8mm, align: center + horizon,
  ..((a0: (0, 0), a1: (0, 1)), (a0: (0, 0), a1: (-1, 1))).map(positions => {
    render-smiles("[13CH3:7]O", skeletal: true, config: (
      atom-positions: positions, highlight-atoms: ("a0", "a1"),
      highlight-bonds: (("a0", "a1"),), collision-policy: "error",
    ))
  }))
