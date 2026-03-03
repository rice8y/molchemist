# molchemist

**molchemist** is a Typst package for rendering chemical structures directly from Molfile (`.mol`) and Structure-Data File (`.sdf`) formats. 

It leverages a blazing-fast Rust/WASM core (powered by [`sdfrust`](https://github.com/hfooladi/sdfrust)) to parse molecular graphs and detect cycles, and seamlessly renders them using the declarative drawing engine of [`alchemist`](https://typst.app/universe/package/alchemist).

## Usage

Import the `render-mol` function from the package and pass the raw string data of your `.mol` or `.sdf` file.

```typ
#import "molchemist.typ": render-mol

// Read your molecule data
#let mol-data = read("molecule.sdf")
```

### Rendering Modes

`molchemist` supports three distinct rendering styles to suit your document's needs:

#### 1. Full Mode (Default)

Draws every single atom and bond explicitly exactly as defined in the source file, including all carbons and hydrogens.

```typ
#render-mol(mol-data)
```

#### 2. Abbreviated Mode

A standard chemical representation. It hides the carbon backbone, wraps explicit hydrogens into their parent heteroatoms (e.g., `O` + `H` becomes `OH`), and neatly formats terminal carbon groups (e.g., `CH3`).

```typ
#render-mol(mol-data, abbreviate: true)
```

#### 3. Skeletal Mode

A pure skeletal formula. All backbone carbons and their attached hydrogens are completely hidden, leaving only the zigzag lines and heteroatoms.

```typ
#render-mol(mol-data, skeletal: true)
```

### Customizing Appearance

Under the hood, `molchemist` passes styling arguments directly to `alchemist`'s `skeletize` function. You can customize the look of your molecules using the `config` dictionary.

```typ
#render-mol(
  mol-data, 
  skeletal: true,
  config: (
    atom-sep: 2em,
    fragment-color: blue,
    angle-increment: 30deg
  )
)
```

## API Reference

* **`data`** (`str`): The raw string content of a `.mol` or `.sdf` file.
* **`abbreviate`** (`bool`): If `true`, applies standard chemical abbreviations (e.g., folding H into heteroatoms, labeling terminal CH3). Default is `false`.
* **`skeletal`** (`bool`): If `true`, renders a pure skeletal structure, overriding `abbreviate`. Hides all backbone C and H atoms. Default is `false`.
* **`config`** (`dictionary`): A dictionary of styling options passed directly to the `alchemist` package.

## Under the Hood

1. **Parse:** The WASM plugin (`sdfrust`) parses the 3D/2D coordinates and bonding topology from the Molfile.
2. **Traverse & Calculate:** A Depth-First Search (DFS) algorithm traverses the graph, calculating absolute angles (`dx, dy -> atan2`) for every bond to ensure accurate placement regardless of automatic layout engines.
3. **Serialize:** The resulting Abstract Syntax Tree (AST) is serialized into a lightweight CBOR binary.
4. **Render:** Typst decodes the CBOR natively and maps the instructions directly to `alchemist` functions (`fragment`, `hook`, `branch`, `single`, `double`, `triple`).

## License

This project is distributed under the MIT License. See [LICENSE](LICENSE) for details.