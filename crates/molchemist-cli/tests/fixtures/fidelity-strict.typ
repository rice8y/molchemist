#import "../../../../package/lib.typ": render-mol

#let sdf = "strict fidelity\n  molchemist\n\n  1  0  0  0  0  0  0  0  0  0  0 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  APO  1   1   1\nM  END\n"

#render-mol(sdf, fidelity: "strict")
