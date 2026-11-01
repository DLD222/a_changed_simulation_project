#!/usr/bin/env python3
"""T7 (Meep side): 2D ring resonator resonance frequency vs resolution.

Geometry and source follow the standard Meep ring example: a high-index ring
(n=3.4) side-coupled to a straight waveguide, excited by a broadband point source
in the waveguide, with the resonance extracted by harminv.

Usage (needs the Meep environment):
  ~/miniconda3/envs/meep313/bin/python meep_ring.py --res 10 20 40
"""
import argparse
import json

import meep as mp


def run(res):
    n = 3.4
    w = 1.0
    r = 1.0
    pad = 4.0
    dpml = 2.0
    sx = 2 * (r + pad + dpml)
    cell = mp.Vector3(sx, sx)

    geometry = [
        mp.Block(mp.Vector3(mp.inf, w, mp.inf),
                 material=mp.Medium(epsilon=n ** 2),
                 center=mp.Vector3(0, -r - w / 2)),
        mp.Cylinder(r, material=mp.Medium(epsilon=n ** 2)),
    ]
    fcen, df = 0.15, 0.10
    sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df),
                         component=mp.Ez,
                         center=mp.Vector3(0, -r - w / 2))]
    sim = mp.Simulation(cell_size=cell, geometry=geometry, sources=sources,
                        boundary_layers=[mp.PML(dpml)], resolution=res)
    h = mp.Harminv(mp.Ez, mp.Vector3(0, r, 0), fcen, df)
    sim.run(mp.after_sources(h), until_after_sources=2000)

    modes = sorted([(m.freq, m.decay) for m in h.modes], key=lambda t: -abs(t[1]))
    return {"resolution": res, "sx": sx,
            "modes": [{"freq": f, "decay": d, "Q": (3.141592653589793 * f / d) if d else None}
                      for f, d in modes]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--res", type=int, nargs="+", default=[10, 20, 40])
    ap.add_argument("--out", default="/tmp/t7_meep.json")
    args = ap.parse_args()

    out = []
    for res in args.res:
        print(f"--- resolution {res} ---", flush=True)
        res_d = run(res)
        out.append(res_d)
        for m in res_d["modes"][:3]:
            print(f"   f={m['freq']:.6f}  decay={m['decay']:.3e}  Q={m['Q']:.1f}")
    with open(args.out, "w") as f:
        json.dump(out, f, indent=1)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
