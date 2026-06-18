#import "@preview/mantys:1.0.2": *
#import "@preview/codly:1.3.0"
#import "../lib.typ" as molchemist
#let cetz = molchemist.cetz

#let infos = toml("../typst.toml")
#let docs-cid-241-sdf = read("assets/Structure2D_COMPOUND_CID_241.sdf")
#let docs-cid-93406-sdf = read("assets/Structure2D_COMPOUND_CID_93406.sdf")
#let docs-sid-93298-sdf = read("assets/DepositedStructure_SUBSTANCE_SID_93298_Version_3.sdf")

#show: mantys(
  ..infos,
  title: [#infos.package.name],
  subtitle: [Molecule rendering from Molfile, SDF, and SMILES],
  date: datetime.today(),
  abstract: [
    `molchemist` renders chemical structures from Molfile / SDF data. It can also parse SMILES strings by generating a 2D layout first. The package turns molecular graphs into an `alchemist` drawing program and renders the final figure with CeTZ.

    The package is aimed at Typst documents that need compact molecule figures, publication-oriented skeletal formulae, and light annotation without leaving Typst.
  ],
  wrap-snippets: true,
  examples-scope: (
    scope: (
      molchemist: molchemist,
      cetz: cetz,
      docs-cid-241-sdf: docs-cid-241-sdf,
      docs-cid-93406-sdf: docs-cid-93406-sdf,
      docs-sid-93298-sdf: docs-sid-93298-sdf,
    ),
    imports: (
      molchemist: "*",
    ),
  ),
  theme: create-theme(
    fonts: (
      serif: ("Times New Roman", "Georgia"),
      sans: ("Helvetica Neue", "Arial"),
      mono: ("Menlo", "Courier New"),
    ),
    text: (
      size: 11pt,
      font: ("Times New Roman", "Georgia"),
      fill: rgb(35, 31, 32),
    ),
    heading: (
      font: ("Helvetica Neue", "Arial"),
      fill: rgb(35, 31, 32),
    ),
    emph: (
      link: rgb("#1f4f73"),
    ),
    code: (
      size: 9pt,
      font: ("Menlo", "Courier New"),
      fill: rgb("#555555"),
    ),
  ),
)

#let example = example.with(side-by-side: false, breakable: true)
#let doc-code(..args, body) = frame(
  breakable: true,
  codly.local(number-format: none, breakable: true, ..args, body),
)

#import molchemist: *

= Getting Started

Import `molchemist` and choose the renderer that matches your input: @cmd:render-mol[-] for Molfile or SDF text, and @cmd:render-smiles[-] for inline SMILES text.

#example[
  ```typ
  #import "@preview/molchemist:0.1.2": *

  #let mol-data = read("Structure2D_COMPOUND_CID_93406.sdf")
  #render-mol(mol-data, abbreviate: true)
  ```
][
  #render-mol(docs-cid-93406-sdf, abbreviate: true)
]

The examples below assume this import unless they need an additional package such as CeTZ.

On Typst 0.15.0 and later, @cmd:render-mol[-] can also receive `path("molecule.sdf")` directly. This manual keeps using `read(...)` in examples for compatibility with older Typst versions.

SMILES is useful for compact inline examples or generated documents. Because SMILES stores connectivity rather than drawing coordinates, `molchemist` first computes a 2D layout and then renders the structure.

#example(```typ
#render-smiles(
  "CCC1=C(N=C2C=CC=C(N2C1=O)C)C",
  abbreviate: true,
)
```)

== Choosing an Input Format

- `Molfile / SDF`: coordinate-bearing structure text. Use it for database exports or drawing-tool output when preserving the supplied layout matters.
- `SMILES`: compact inline source text. Use it for examples, generated documents, and quick sketches. `molchemist` computes 2D coordinates before rendering.

= Example Data

The SDF examples in this manual use real PubChem records included in the repository test data.

- #link("https://pubchem.ncbi.nlm.nih.gov/compound/241")[PubChem CID 241]: benzene.
- #link("https://pubchem.ncbi.nlm.nih.gov/compound/93406")[PubChem CID 93406]: `3-ethyl-2,6-dimethylpyrido[1,2-a]pyrimidin-4-one`.
- #link("https://pubchem.ncbi.nlm.nih.gov/substance/93298")[PubChem SID 93298]: a deposited DTP/NCI substance record associated with CID 235403.

#example[
  ```typ
  #let benzene = read("Structure2D_COMPOUND_CID_241.sdf")
  #let fused = read("Structure2D_COMPOUND_CID_93406.sdf")
  #let deposited = read("DepositedStructure_SUBSTANCE_SID_93298_Version_3.sdf")

  #grid(
    columns: 3,
    gutter: 7mm,
    align: horizon + center,
    render-mol(benzene, skeletal: true, config: (atom-sep: 1.55em)),
    render-mol(fused, skeletal: true, config: (atom-sep: 1.55em)),
    render-mol(deposited, skeletal: true, config: (atom-sep: 1.55em)),
  )
  ```
][
  #grid(
    columns: 3,
    gutter: 7mm,
    align: horizon + center,
    render-mol(docs-cid-241-sdf, skeletal: true, config: (atom-sep: 1.55em)),
    render-mol(docs-cid-93406-sdf, skeletal: true, config: (atom-sep: 1.55em)),
    render-mol(docs-sid-93298-sdf, skeletal: true, config: (atom-sep: 1.55em)),
  )
]

= Rendering Modes

Both renderers support the same three modes. Full mode is useful for small molecules and debugging. Abbreviated mode folds common hydrogens and terminal carbons into labels. Skeletal mode is usually the best default for publication figures.

#example[
  ```typ
  #let mol-data = read("Structure2D_COMPOUND_CID_241.sdf")

  #grid(
    columns: 3,
    gutter: 8mm,
    align: horizon + center,
    render-mol(mol-data),
    render-mol(mol-data, abbreviate: true),
    render-mol(mol-data, skeletal: true),
  )
  ```
][
  #grid(
    columns: 3,
    gutter: 8mm,
    align: horizon + center,
    render-mol(docs-cid-241-sdf),
    render-mol(docs-cid-241-sdf, abbreviate: true),
    render-mol(docs-cid-241-sdf, skeletal: true),
  )
]

= Styling

Pass #arg[config] for visual adjustments that should reach `alchemist`, such as atom spacing, fragment color, and bond styles.

#example[
  ```typ
  #let mol-data = read("Structure2D_COMPOUND_CID_93406.sdf")

  #render-mol(
    mol-data,
    skeletal: true,
    config: (
      atom-sep: 2.6em,
      fragment-color: rgb("#1b4d5a"),
      single: (stroke: 0.75pt + rgb("#1b4d5a")),
      double: (stroke: 0.75pt + rgb("#b14b2f")),
    ),
  )
  ```
][
  #render-mol(
    docs-cid-93406-sdf,
    skeletal: true,
    config: (
      atom-sep: 2.6em,
      fragment-color: rgb("#1b4d5a"),
      single: (stroke: 0.75pt + rgb("#1b4d5a")),
      double: (stroke: 0.75pt + rgb("#b14b2f")),
    ),
  )
]

#info-alert[
  Molfile and SDF input already contains coordinates. Routing-oriented `alchemist` settings do not reshape that graph; prefer `atom-sep`, font, and stroke settings for visual tuning.
]

= Annotations

Annotations are drawn after the molecule. They are meant for sparse publication callouts, not for replacing a full figure editor.

== Finding Atom and Bond Indices

Enable #arg[show-indices] while authoring. The visible labels are the indices used by @cmd:atom-anchor[-] and @cmd:bond-anchor[-].

#example[
  ```typ
  #let mol-data = read("Structure2D_COMPOUND_CID_93406.sdf")
  #render-mol(mol-data, abbreviate: true, show-indices: true)
  ```
][
  #render-mol(docs-cid-93406-sdf, abbreviate: true, show-indices: true)
]

== Callouts

Use @cmd:callout-annotation[-] for external labels. Defaults are intentionally quiet: no arrowhead, a thin leader line, and an unboxed label.

#example[
  ```typ
  #let mol-data = read("Structure2D_COMPOUND_CID_93406.sdf")

  #render-mol(
    mol-data,
    abbreviate: true,
    annotations: (
      callout-annotation(
        bond-anchor(0),
        [carbonyl],
        label-at: (to: molecule-anchor(anchor: "south"), rel: (0.82, -0.48)),
        leader-start: (to: molecule-anchor(anchor: "south"), rel: (0.7, -0.48)),
        leader-end: (to: bond-anchor(0), rel: (0.2, -0.24)),
        leader: "curve",
        stroke: luma(58%) + 0.28pt,
        label-size: 0.76em,
        label-anchor: "west",
      ),
      callout-annotation(
        atom-anchor(12),
        [methyl substituent],
        label-at: (to: molecule-anchor(anchor: "east"), rel: (0.9, 1.04)),
        leader-start: (to: molecule-anchor(anchor: "east"), rel: (0.78, 1.04)),
        leader-end: (to: atom-anchor(12), rel: (0.28, -0.38)),
        leader: "curve",
        stroke: luma(58%) + 0.28pt,
        label-size: 0.76em,
        label-anchor: "west",
      ),
    ),
  )
  ```
][
  #render-mol(
    docs-cid-93406-sdf,
    abbreviate: true,
    annotations: (
      callout-annotation(
        bond-anchor(0),
        [carbonyl],
        label-at: (to: molecule-anchor(anchor: "south"), rel: (0.82, -0.48)),
        leader-start: (to: molecule-anchor(anchor: "south"), rel: (0.7, -0.48)),
        leader-end: (to: bond-anchor(0), rel: (0.2, -0.24)),
        leader: "curve",
        stroke: luma(58%) + 0.28pt,
        label-size: 0.76em,
        label-anchor: "west",
      ),
      callout-annotation(
        atom-anchor(12),
        [methyl substituent],
        label-at: (to: molecule-anchor(anchor: "east"), rel: (0.9, 1.04)),
        leader-start: (to: molecule-anchor(anchor: "east"), rel: (0.78, 1.04)),
        leader-end: (to: atom-anchor(12), rel: (0.28, -0.38)),
        leader: "curve",
        stroke: luma(58%) + 0.28pt,
        label-size: 0.76em,
        label-anchor: "west",
      ),
    ),
  )
]

== Arrows

Use @cmd:arrow-annotation[-] for molecule-level process arrows or simple directional marks.

#example[
  ```typ
  #let mol-data = read("Structure2D_COMPOUND_CID_93406.sdf")

  #render-mol(
    mol-data,
    skeletal: true,
    annotations: arrow-annotation(
      (to: molecule-anchor(anchor: "east"), rel: (0.48, 0)),
      (to: molecule-anchor(anchor: "east"), rel: (2.65, 0)),
      label: [derivatization],
      label-offset: (0, -0.46),
      label-anchor: "north",
    ),
  )
  ```
][
  #render-mol(
    docs-cid-93406-sdf,
    skeletal: true,
    annotations: arrow-annotation(
      (to: molecule-anchor(anchor: "east"), rel: (0.48, 0)),
      (to: molecule-anchor(anchor: "east"), rel: (2.65, 0)),
      label: [derivatization],
      label-offset: (0, -0.46),
      label-anchor: "north",
    ),
  )
]

== Reaction Schemes

For reactions, compose separate molecule renderings with normal Typst/CeTZ layout. This keeps reaction arrows independent from molecule annotations.

#example(```typ
#let reaction-arrow(above, below: none) = cetz.canvas({
  import cetz.draw: *
  line(
    (0.1, 0),
    (2.95, 0),
    stroke: 0.65pt + black,
    mark: (end: ">>", scale: 0.72, fill: black),
  )
  content((1.52, 0.42), text(size: 0.78em)[#above], anchor: "south")
  if below != none {
    content((1.52, -0.38), text(size: 0.72em)[#below], anchor: "north")
  }
})

#grid(
  columns: (auto, 28mm, auto),
  column-gutter: 4mm,
  align: horizon + center,
  render-smiles("O=CC1=CC=CC=C1", abbreviate: true),
  reaction-arrow([NaBH#sub[4]], below: [MeOH]),
  render-smiles("OCC1=CC=CC=C1", abbreviate: true),
)
```)

== Mechanism Example: von Richter Reaction

Long mechanisms can combine SMILES rendering, @cmd:cetz-annotation[-] for electron movement, and ordinary CeTZ layout. This example shows the classical conversion of p-bromonitrobenzene to m-bromobenzoic acid.

#example[
  ```typ
  // Electron-flow helpers.
  #let electron-arrows(body) = cetz-annotation(mol => {
    import cetz.draw: *
    body(mol)
  })

  #let electron-arrow(
    mol,
    from,
    to,
    from-offset: (0, 0),
    to-offset: (0, 0),
    controls: ((0.35, 0.6), (0.35, 0.6)),
  ) = {
    import cetz.draw: *
    let start = (to: (name: mol, anchor: from), rel: from-offset)
    let end = (to: (name: mol, anchor: to), rel: to-offset)
    if controls.len() == 1 {
      bezier(
        start,
        end,
        (to: start, rel: controls.at(0)),
        stroke: 0.48pt + black,
        mark: (end: ">>", scale: 0.62, fill: black),
      )
    } else {
      bezier(
        start,
        end,
        (to: start, rel: controls.at(0)),
        (to: end, rel: controls.at(1)),
        stroke: 0.48pt + black,
        mark: (end: ">>", scale: 0.62, fill: black),
      )
    }
  }

  #let reagent(mol, at, label, offset: (0, 0)) = {
    import cetz.draw: *
    content(
      (to: (name: mol, anchor: at), rel: offset),
      text(size: 8pt)[#label],
      anchor: "center",
    )
  }

  // Reaction-arrow and row-layout helpers.
  #let reaction-arrow(above: none, below: none) = cetz.canvas({
    import cetz.draw: *
    line(
      (0.08, 0),
      (1.08, 0),
      stroke: 0.6pt + black,
      mark: (end: ">>", scale: 0.68, fill: black),
    )
    if above != none {
      content((0.58, 0.3), text(size: 7.2pt)[#above], anchor: "south")
    }
    if below != none {
      content((0.58, -0.27), text(size: 6.8pt)[#below], anchor: "north")
    }
  })

  #let molecule(smiles, annotations: none) = box(
    width: 32mm,
    align(center)[
      #set text(size: 8pt)
      #render-smiles(
        smiles,
        abbreviate: true,
        config: (atom-sep: 1.7em),
        annotations: annotations,
      )
    ],
  )

  #let mechanism-row(..items) = {
    let cells = items.pos()
    grid(
      columns: cells.map(
        item => if item.at(0) == "molecule" { 32mm } else { 12mm },
      ),
      column-gutter: 0.6mm,
      align: horizon + center,
      ..cells.map(item => item.at(1)),
    )
  }

  // Substrate and early intermediates.
  #let substrate = molecule(
    "O=[N+]([O-])c1ccc(Br)cc1",
    annotations: electron-arrows(mol => {
      import cetz.draw: *
      reagent(mol, "east", [CN#super[-]], offset: (0.52, 0.22))
      electron-arrow(
        mol,
        "east",
        "b4.50%",
        from-offset: (0.4, 0.42),
        to-offset: (0.1, 0.16),
        controls: ((0, 0.48), (0.18, 0.68)),
      )
    }),
  )

  #let sigma-adduct = molecule(
    "O=[N+]([O-])C1=CC(Br)=CC=C1C#N",
    annotations: electron-arrows(mol => {
      electron-arrow(
        mol,
        "a2.south",
        "b10.20%",
        from-offset: (-0.05, -0.08),
        to-offset: (0.18, 0),
        controls: ((0, -0.4),),
      )
    }),
  )

  #let cyclic-imidate = molecule("N=C1ON(=O)c2ccc(Br)cc21")

  // Later intermediates and products.
  #let nitroso-amide = molecule("O=Nc1c(C(=O)N)cc(Br)cc1")
  #let hydroxy-azo = molecule("O=C1NN(O)C2=C1C=CC(Br)=C2")

  #let azoketone = molecule(
    "O=C1N=NC2=C1C=CC(Br)=C2",
    annotations: electron-arrows(mol => {
      import cetz.draw: *
      reagent(mol, "east", [OH#super[-]], offset: (0.72, 0.18))
      electron-arrow(
        mol,
        "east",
        "b0.50%",
        from-offset: (0.48, 0.34),
        to-offset: (0.28, 0.38),
        controls: ((0, 0.36), (0.16, 0.5)),
      )
    }),
  )

  #let carboxylate = molecule("O=C([O-])c1cc(Br)ccc1")
  #let product = molecule("O=C(O)c1cc(Br)ccc1")

  // Compose the mechanism.
  #block(breakable: false, width: 100%)[
    #stack(
      dir: ttb,
      spacing: 5mm,
      mechanism-row(
        ("molecule", substrate),
        ("arrow", reaction-arrow(above: [KCN], below: [EtOH / H#sub[2]O])),
        ("molecule", sigma-adduct),
        ("arrow", move(dy: -2.1mm, reaction-arrow(above: [cyclization]))),
        ("molecule", cyclic-imidate),
      ),
      mechanism-row(
        ("arrow", reaction-arrow(above: [ring opening])),
        ("molecule", nitroso-amide),
        ("arrow", reaction-arrow(above: [cyclization])),
        ("molecule", hydroxy-azo),
        ("arrow", reaction-arrow(above: [-H#sub[2]O])),
        ("molecule", azoketone),
      ),
      align(center)[
        #mechanism-row(
          ("arrow", reaction-arrow(above: [OH#super[-]], below: [-N#sub[2]])),
          ("molecule", carboxylate),
          ("arrow", move(dy: -2mm, reaction-arrow(above: [H#super[+]]))),
          ("molecule", product),
        )
      ],
    )
  ]
  ```
]

#let von-richter-mechanism() = {
  let electron-arrows(body) = cetz-annotation(mol => {
    import cetz.draw: *
    body(mol)
  })

  let electron-arrow(
    mol,
    from,
    to,
    from-offset: (0, 0),
    to-offset: (0, 0),
    controls: ((0.35, 0.6), (0.35, 0.6)),
  ) = {
    import cetz.draw: *
    let start = (to: (name: mol, anchor: from), rel: from-offset)
    let end = (to: (name: mol, anchor: to), rel: to-offset)
    if controls.len() == 1 {
      bezier(
        start,
        end,
        (to: start, rel: controls.at(0)),
        stroke: 0.48pt + black,
        mark: (end: ">>", scale: 0.62, fill: black),
      )
    } else {
      bezier(
        start,
        end,
        (to: start, rel: controls.at(0)),
        (to: end, rel: controls.at(1)),
        stroke: 0.48pt + black,
        mark: (end: ">>", scale: 0.62, fill: black),
      )
    }
  }

  let reagent(mol, at, label, offset: (0, 0)) = {
    import cetz.draw: *
    content(
      (to: (name: mol, anchor: at), rel: offset),
      text(size: 8pt)[#label],
      anchor: "center",
    )
  }

  let reaction-arrow(above: none, below: none) = cetz.canvas({
    import cetz.draw: *
    line(
      (0.08, 0),
      (1.08, 0),
      stroke: 0.6pt + black,
      mark: (end: ">>", scale: 0.68, fill: black),
    )
    if above != none {
      content((0.58, 0.3), text(size: 7.2pt)[#above], anchor: "south")
    }
    if below != none {
      content((0.58, -0.27), text(size: 6.8pt)[#below], anchor: "north")
    }
  })

  let molecule(smiles, annotations: none) = box(
    width: 32mm,
    align(center)[
      #set text(size: 8pt)
      #render-smiles(
        smiles,
        abbreviate: true,
        config: (atom-sep: 1.7em),
        annotations: annotations,
      )
    ],
  )

  let mechanism-row(..items) = {
    let cells = items.pos()
    grid(
      columns: cells.map(
        item => if item.at(0) == "molecule" { 32mm } else { 12mm },
      ),
      column-gutter: 0.6mm,
      align: horizon + center,
      ..cells.map(item => item.at(1)),
    )
  }

  let substrate = molecule(
    "O=[N+]([O-])c1ccc(Br)cc1",
    annotations: electron-arrows(mol => {
      import cetz.draw: *
      reagent(mol, "east", [CN#super[-]], offset: (0.52, 0.22))
      electron-arrow(
        mol,
        "east",
        "b4.50%",
        from-offset: (0.4, 0.42),
        to-offset: (0.1, 0.16),
        controls: ((0, 0.48), (0.18, 0.68)),
      )
    }),
  )

  let sigma-adduct = molecule(
    "O=[N+]([O-])C1=CC(Br)=CC=C1C#N",
    annotations: electron-arrows(mol => {
      electron-arrow(
        mol,
        "a2.south",
        "b10.20%",
        from-offset: (-0.05, -0.08),
        to-offset: (0.18, 0),
        controls: ((0, -0.4),),
      )
    }),
  )

  let cyclic-imidate = molecule("N=C1ON(=O)c2ccc(Br)cc21")

  let nitroso-amide = molecule("O=Nc1c(C(=O)N)cc(Br)cc1")

  let hydroxy-azo = molecule("O=C1NN(O)C2=C1C=CC(Br)=C2")

  let azoketone = molecule(
    "O=C1N=NC2=C1C=CC(Br)=C2",
    annotations: electron-arrows(mol => {
      import cetz.draw: *
      reagent(mol, "east", [OH#super[-]], offset: (0.72, 0.18))
      electron-arrow(
        mol,
        "east",
        "b0.50%",
        from-offset: (0.48, 0.34),
        to-offset: (0.28, 0.38),
        controls: ((0, 0.36), (0.16, 0.5)),
      )
    }),
  )

  let carboxylate = molecule("O=C([O-])c1cc(Br)ccc1")
  let product = molecule("O=C(O)c1cc(Br)ccc1")

  block(breakable: false, width: 100%)[
    #stack(
      dir: ttb,
      spacing: 5mm,
      mechanism-row(
        ("molecule", substrate),
        ("arrow", reaction-arrow(above: [KCN], below: [EtOH / H#sub[2]O])),
        ("molecule", sigma-adduct),
        ("arrow", move(dy: -2.1mm, reaction-arrow(above: [cyclization]))),
        ("molecule", cyclic-imidate),
      ),
      mechanism-row(
        ("arrow", reaction-arrow(above: [ring opening])),
        ("molecule", nitroso-amide),
        ("arrow", reaction-arrow(above: [cyclization])),
        ("molecule", hydroxy-azo),
        ("arrow", reaction-arrow(above: [-H#sub[2]O])),
        ("molecule", azoketone),
      ),
      align(center)[
        #mechanism-row(
          ("arrow", reaction-arrow(above: [OH#super[-]], below: [-N#sub[2]])),
          ("molecule", carboxylate),
          ("arrow", move(dy: -2mm, reaction-arrow(above: [H#super[+]]))),
          ("molecule", product),
        )
      ],
    )
  ]
}

The mechanistic sequence follows M. Rosenblum, _The Mechanism of the von Richter Reaction_, J. Am. Chem. Soc. 82 (1960), 3796-3798, #link("https://doi.org/10.1021/ja01499a090")[doi:10.1021/ja01499a090]. Molecule-specific anchors and curved-arrow routing can be adjusted directly in the figure source.

== Low-Level Labels

Use @cmd:label-annotation[-] when no leader line is needed.

#example[
  ```typ
  #let mol-data = read("DepositedStructure_SUBSTANCE_SID_93298_Version_3.sdf")

  #render-mol(
    mol-data,
    skeletal: true,
    annotations: label-annotation(
      molecule-anchor(anchor: "south"),
      [PubChem SID 93298],
      offset: (0, -0.65),
      label-anchor: "north",
    ),
  )
  ```
][
  #render-mol(
    docs-sid-93298-sdf,
    skeletal: true,
    annotations: label-annotation(
      molecule-anchor(anchor: "south"),
      [PubChem SID 93298],
      offset: (0, -0.65),
      label-anchor: "north",
    ),
  )
]

== Custom CeTZ Overlays

For final figure polishing, @cmd:cetz-annotation[-] exposes the generated molecule name, so you can draw directly against CeTZ anchors.

#example[
  ```typ
  #let mol-data = read("DepositedStructure_SUBSTANCE_SID_93298_Version_3.sdf")

  #render-mol(
    mol-data,
    skeletal: true,
    annotations: cetz-annotation(mol => {
      import cetz.draw: *
      content(
        (to: (name: mol, anchor: "north"), rel: (0, 0.54)),
        text(size: 0.82em)[database structure],
        anchor: "south",
      )
    }),
  )
  ```
][
  #render-mol(
    docs-sid-93298-sdf,
    skeletal: true,
    annotations: cetz-annotation(mol => {
      import cetz.draw: *
      content(
        (to: (name: mol, anchor: "north"), rel: (0, 0.54)),
        text(size: 0.82em)[database structure],
        anchor: "south",
      )
    }),
  )
]

= Dump Mode

With #arg[dump], `molchemist` returns generated `alchemist` source instead of a drawing. Use this when a figure needs manual surgery beyond the annotation API.

#example[
  ```typ
  #let mol-data = read("Structure2D_COMPOUND_CID_241.sdf")
  #render-mol(mol-data, skeletal: true, dump: true)
  ```
][
  #render-mol(docs-cid-241-sdf, skeletal: true, dump: true)
]

= Publication Guidance

For paper figures, start with #arg[skeletal] for hydrocarbon-heavy structures and #arg[abbreviate] when heteroatom hydrogens or terminal groups should remain explicit. Keep full mode for small structures and debugging.

Keep annotations sparse. In most cases, a thin unboxed @cmd:callout-annotation[-] is better than a boxed label or an arrowhead. If a leader line looks like a chemical bond, move the label with #arg[label-at] or stop the leader outside the structure with #arg[leader-end].

#warning-alert[
  Dense structures can still overlap in full mode, especially after SMILES implicit hydrogens are expanded. Prefer abbreviated or skeletal mode when a molecule is meant for a publication figure.
]

= SMILES Notes

SMILES support is a parse-and-layout pipeline. The parser accepts common SMILES notation, aromatic rings, charges, tetrahedral `@` / `@@` centers, and `/` / `\` double-bond geometry.

#example(```typ
#grid(
  columns: 2,
  gutter: 10mm,
  align: horizon + center,
  render-smiles("O=[N+]([O-])c1ccccc1", abbreviate: true),
  render-smiles("N[C@@H](C)C(=O)O", abbreviate: true),
)
```)

Extended OpenSMILES chirality classes such as `@AL`, `@SP`, `@TB`, and `@OH` are accepted, but they are currently preserved as textual stereo annotations rather than native geometric depictions.

= API Reference

#custom-type("anchor", color: aqua)
#custom-type("annotation", color: orange)

== Rendering Functions

#command(
  "render-mol",
  arg("data"),
  arg(abbreviate: false),
  arg(skeletal: false),
  arg(dump: false),
  arg(config: (:)),
  arg(annotations: none),
  arg(show-indices: false),
  ret: content,
)[
  Render a molecule from raw Molfile or SDF text.

  #argument("data", types: (str, bytes, "path"))[
    Raw `.mol` or `.sdf` data. Typst 0.15.0 and later may pass `path(...)` directly; older versions should pass `read(...)` output.
  ]

  #argument("abbreviate", types: bool, default: false)[
    Enables abbreviated rendering.
  ]

  #argument("skeletal", types: bool, default: false)[
    Enables skeletal rendering. This overrides #arg[abbreviate].
  ]

  #argument("dump", types: bool, default: false)[
    Returns generated `alchemist` source code instead of rendering the molecule.
  ]

  #argument("config", types: dictionary, default: (:))[
    Visual configuration passed to `alchemist`.
  ]

  #argument("annotations", types: ("annotation", array, none), default: none)[
    Optional overlay annotation or array of annotations.
  ]

  #argument("show-indices", types: (bool, str), default: false)[
    Debug overlay for annotation authoring. Use `true`, `"all"`, `"atoms"`, or `"bonds"`.
  ]
]

#command(
  "render-smiles",
  arg("smiles"),
  arg(abbreviate: false),
  arg(skeletal: false),
  arg(dump: false),
  arg(config: (:)),
  arg(annotations: none),
  arg(show-indices: false),
  ret: content,
)[
  Parse a SMILES string, generate 2D coordinates, and render the resulting molecule.

  #argument("smiles", types: str)[
    A SMILES string.
  ]

  #argument("abbreviate", types: bool, default: false)[
    Enables abbreviated rendering after layout generation.
  ]

  #argument("skeletal", types: bool, default: false)[
    Enables skeletal rendering. This overrides #arg[abbreviate].
  ]

  #argument("dump", types: bool, default: false)[
    Returns generated `alchemist` source code instead of rendering the molecule.
  ]

  #argument("config", types: dictionary, default: (:))[
    Visual configuration passed to `alchemist`.
  ]

  #argument("annotations", types: ("annotation", array, none), default: none)[
    Optional overlay annotation or array of annotations.
  ]

  #argument("show-indices", types: (bool, str), default: false)[
    Debug overlay for annotation authoring. Use `true`, `"all"`, `"atoms"`, or `"bonds"`.
  ]
]

== Anchor Helpers

#command("atom-anchor", arg("index"), arg(anchor: "mid"), ret: "anchor")[
  Select an atom anchor for annotation placement.

  #argument("index", types: int)[
    Atom index shown by #arg[show-indices].
  ]

  #argument("anchor", types: str, default: "mid")[
    CeTZ anchor on the atom object, such as `"north"`, `"east"`, or `"mid"`.
  ]
]

#command("bond-anchor", arg("index"), arg(anchor: "50%"), ret: "anchor")[
  Select a bond anchor for annotation placement.

  #argument("index", types: int)[
    Bond index shown by #arg[show-indices].
  ]

  #argument("anchor", types: str, default: "50%")[
    CeTZ anchor on the bond object.
  ]
]

#command("molecule-anchor", arg(anchor: "center"), ret: "anchor")[
  Select an anchor on the whole rendered molecule.

  #argument("anchor", types: str, default: "center")[
    CeTZ group anchor such as `"center"`, `"north"`, `"east"`, `"south"`, or `"west"`.
  ]
]

#command("atom-ref", arg("index"), ret: str)[
  Return the internal atom anchor name as a string.
]

#command("bond-ref", arg("index"), ret: str)[
  Return the internal bond anchor name as a string.
]

== Annotation Builders

#command(
  "callout-annotation",
  arg("at"),
  arg("label"),
  arg(anchor: "mid"),
  arg(side: "north-east"),
  arg(label-at: auto),
  arg(leader: "curve"),
  arg(mark: none),
  ret: "annotation",
)[
  Build an external label with a leader line.

  #argument("at", types: "anchor")[
    Target anchor.
  ]

  #argument("label", types: content)[
    Label content.
  ]

  #argument("side", types: str, default: "north-east")[
    Preset label side. Supported values are `"east"`, `"west"`, `"north"`, `"south"`, and diagonal combinations.
  ]

  #argument("label-at", types: (auto, dictionary), default: auto)[
    Explicit CeTZ coordinate for the label.
  ]

  #argument("leader", types: str, default: "curve")[
    Leader style: `"curve"`, `"straight"`, or `"elbow"`.
  ]

  #argument("leader-end", types: (auto, dictionary), default: auto)[
    Manual CeTZ coordinate where the leader should stop.
  ]

  #argument("target-gap", types: float, default: 0.2)[
    Clearance from the target when #arg[leader-end] is automatic.
  ]
]

#command(
  "arrow-annotation",
  arg("from"),
  arg("to"),
  arg(label: none),
  arg(mark: (end: ">>", fill: black)),
  ret: "annotation",
)[
  Build a free arrow overlay.
]

#command(
  "label-annotation",
  arg("at"),
  arg("label"),
  arg(offset: (0, 0.45)),
  ret: "annotation",
)[
  Build a free text label overlay without a leader line.
]

#command("cetz-annotation", arg("body"), ret: "annotation")[
  Run custom CeTZ code after the molecule is drawn.

  #argument("body", types: function)[
    Function receiving the generated molecule name.
  ]
]

= Limitations and Roadmap

- Full mode can become crowded for large molecules because explicit hydrogens and atom labels occupy real page space.
- SMILES layout is generated internally and may not match the exact layout from an external chemical drawing program.
- Extended OpenSMILES chirality classes are currently preserved as textual stereo annotations rather than native geometric depictions.
- Annotation helpers are intentionally minimal. For complex figure composition, use @cmd:cetz-annotation[-] or dump the generated `alchemist` code.

= License and Dependencies

`molchemist` is distributed under the MIT license. Molfile / SDF parsing is powered by `sdfrust`; SMILES parsing is based on `opensmiles`; SMILES 2D coordinate generation uses `CoordgenLibs`; rendering is handled by `alchemist` and CeTZ. See the third-party notices distributed with the package for full license details.

The example SDF files and rendered example images are attributed separately from the package code. In particular, this manual includes PubChem-derived example structures such as #link("https://pubchem.ncbi.nlm.nih.gov/compound/241")[CID 241]. See `THIRD_PARTY_NOTICES.md` and `docs/assets/README.md` for source URLs and the relevant NCBI data-usage policy.
