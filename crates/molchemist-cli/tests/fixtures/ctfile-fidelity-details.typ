#import "../../../../package/lib.typ": render-mol

#set page(width: auto, height: auto, margin: 3mm)

#render-mol(
  read("ctfile-fidelity.sdf"),
  fidelity: "strict",
  config: (
    ctfile: (
      highlight-paint: rgb("#74c0fc").transparentize(45%),
      query-details: true,
    ),
  ),
)
