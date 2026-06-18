# Alchemist 0.2.0 Release Switch Notes

`package/lib.typ` currently imports the local development copy:

```typ
#import "../artifacts/alchemist/lib.typ": *
```

This is intentional while `alchemist:0.2.0` is not yet available in Typst Universe. Before releasing `molchemist`, switch it to:

```typ
#import "@preview/alchemist:0.2.0": *
```

The implementation depends on the alchemist 0.2.0 anchor model for named molecules, atoms, and links. In particular, `bond-anchor(...)` expects named links such as `b0`, `b1`, ... to be exposed through the named molecule group.
