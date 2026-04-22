# molchemist

**molchemist** is a Typst package for rendering chemical structures from Molfile / SDF data and from SMILES strings.

It uses a Rust/WASM core to parse molecular graphs and generate `alchemist` ASTs, together with a companion WASM layout plugin for SMILES 2D coordinate generation. Molfile / SDF parsing is powered by [`sdfrust`](https://github.com/hfooladi/sdfrust), SMILES parsing is based on [`opensmiles`](https://crates.io/crates/opensmiles), SMILES 2D coordinate generation uses [`CoordgenLibs`](https://github.com/schrodinger/coordgenlibs), and the final rendering is handled by the declarative drawing engine of [`alchemist`](https://github.com/Typsium/alchemist).

Third-party license notices for redistributed dependencies are collected in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

## Usage

Import `render-mol` for Molfile/SDF inputs, or `render-smiles` for SMILES inputs.

```typ
#import "@preview/molchemist:0.1.1": render-mol, render-smiles

// Read your molecule data
// Example: https://pubchem.ncbi.nlm.nih.gov/compound/93406
#let mol-data = read("Structure2D_COMPOUND_CID_93406.sdf")
```

For SMILES, `molchemist` generates a 2D layout internally before sending the structure to `alchemist`.

```typ
// Example: https://pubchem.ncbi.nlm.nih.gov/compound/896
#render-smiles("CC(=O)NCCC1=CNC2=C1C=C(C=C2)OC", abbreviate: true)
```

![SMILES Example](images/ex06.png)

### Rendering Modes

`molchemist` supports three distinct rendering styles to suit your document's needs:

#### 1. Full Mode (Default)

Draws every single atom and bond explicitly exactly as defined in the source file, including all carbons and hydrogens.

*Note: For complex molecules, text overlapping may occur. See [Known Limitations](#known-limitations) for workarounds.*

```typ
#render-mol(mol-data)
```

![Full Mode](images/ex01.png)

#### 2. Abbreviated Mode

A standard chemical representation. It hides the carbon backbone, wraps explicit hydrogens into their parent heteroatoms (e.g., `O` + `H` becomes `OH`), and neatly formats terminal carbon groups (e.g., `CH3`).

```typ
#render-mol(mol-data, abbreviate: true)
```

![Abbreviated Mode](images/ex02.png)

#### 3. Skeletal Mode

A pure skeletal formula. All backbone carbons and their attached hydrogens are completely hidden, leaving only the zigzag lines and heteroatoms.

```typ
#render-mol(mol-data, skeletal: true)
```

![Skeletal Mode](images/ex03.png)

### Customizing Appearance

Under the hood, `molchemist` parses the graph and generates native `alchemist` elements. You can customize the look of your molecules by passing styling arguments via the `config` dictionary, which are passed directly to `alchemist`'s `skeletize` function.

```typ
#render-mol(
  mol-data, 
  skeletal: true,
  config: (
    atom-sep: 2em,
    fragment-margin: 2pt,
    fragment-color: blue,
    fragment-font: "New Computer Modern",
    single: (stroke: 1pt + black),
    double: (gap: 0.3em, stroke: 1pt + red)
  )
)
```

![Custom Appearance](images/ex04.png)

**Important Note on Configuration:**

- **Routing overrides:** Because `molchemist` maps the exact 2D absolute coordinates from the source `.sdf`/`.mol` file, `alchemist`'s automatic routing configs (like `angle-increment`, `base-angle`) are bypassed and have no effect.
- **Lewis Structures:** `molchemist` does not automatically infer or generate Lewis structures from SDF files, so `lewis-*` configs are not applicable out of the box.

### Advanced: Ejecting to Alchemist Code (Dump Mode)

If you need to manually fine-tune a molecule, add a specific Lewis structure, or integrate the structure into a larger custom `alchemist` drawing, you can use the `dump` parameter.

When `dump: true` is passed, `molchemist` will not render the molecule. Instead, it will output the generated native `alchemist` code block into your document. You can then copy, paste, and modify this code directly.

```typ
#render-mol(mol-data, dump: true)
```

![Dump Mode](images/ex05.png)

## Known Limitations

When rendering highly complex or dense molecules (e.g., polycyclic compounds, dense substituents) in the default **Full Mode**, you may encounter overlapping atoms or intersecting bonds. This occurs because the 2D absolute coordinates provided in the source `.sdf`/`.mol` files might not allocate enough physical space on the canvas to draw every explicit text label without collisions.

**Recommended Workarounds:**

1. **Use Abbreviated or Skeletal Mode:** For complex organic structures, it is highly recommended to set `abbreviate: true` or `skeletal: true`. This hides redundant atoms, dramatically improving readability and preventing overlaps, which aligns with standard chemical drawing practices.
2. **Increase Bond Length:** If you strictly require Full Mode, you can increase the distance between atoms to create more physical space for the text labels by adjusting the `atom-sep` property in the `config` argument:
    ```typ
    // The default atom-sep is 3em
    #render-mol(mol-data, config: (atom-sep: 4.5em))
    ```

For SMILES input, the default `render-smiles(...)` mode expands implicit hydrogens into explicit `H` atoms so that `full` mode stays closer to the behavior of `render-mol(...)`. Highly complex or dense molecules can still become visually busy in `full` mode, so `abbreviate: true` or `skeletal: true` will often produce a clearer result. The current implementation also supports tetrahedral `@` / `@@` centers and `/` / `\` double-bond geometry as stereochemical depictions. Extended OpenSMILES chirality classes such as `@AL`, `@SP`, `@TB`, and `@OH` are accepted as well; because the current `alchemist`-based renderer does not have native glyphs for those geometries, they are preserved as stereo annotations below the rendered structure instead of wedge/dash depictions.

## API Reference

### `render-mol(data, abbreviate: false, skeletal: false, dump: false, config: (:))`

- **`data`** (`str`): The raw string content of a `.mol` or `.sdf` file.
- **`abbreviate`** (`bool`): If `true`, applies standard chemical abbreviations (e.g., folding H into heteroatoms, labeling terminal CH3). Default is `false`.
- **`skeletal`** (`bool`): If `true`, renders a pure skeletal structure, overriding `abbreviate`. Hides all backbone C and H atoms. Default is `false`.
- **`dump`** (`bool`): If `true`, outputs the generated `alchemist` source code as a formatted Typst code block instead of rendering the molecule graphic. Useful for manual tweaking. Default is `false`.
- **`config`** (`dictionary`): A dictionary of visual styling options passed directly to the `alchemist` package.

### `render-smiles(smiles, abbreviate: false, skeletal: false, dump: false, config: (:))`

- **`smiles`** (`str`): A SMILES string to be parsed, laid out in 2D, and rendered.
- **`abbreviate`** (`bool`): If `true`, applies standard chemical abbreviations after SMILES parsing and layout generation. Default is `false`.
- **`skeletal`** (`bool`): If `true`, renders a pure skeletal structure, overriding `abbreviate`. Default is `false`.
- **`dump`** (`bool`): If `true`, outputs the generated `alchemist` source code as a formatted Typst code block instead of rendering the molecule graphic. Default is `false`.
- **`config`** (`dictionary`): A dictionary of visual styling options passed directly to the `alchemist` package.

## License

This project is distributed under the MIT License. See [LICENSE](LICENSE) for details.
