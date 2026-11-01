#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""T7 (Meep side): dielectric disk resonances vs resolution.

Geometry: n=3.4 disk of radius r=1 in vacuum, PML, Ez point source at (0.5, 0)
inside the disk. harminv extracts the resonances over 0.10 < f < 0.60.

Note: the search window matters. An earlier ring attempt returned zero modes because
harminv was asked for 0.05 < f < 0.25 while the disk modes sit at 0.23-0.35.

Usage:
  ~/miniconda3/envs/meep313/bin/python meep_disk.py --res 10 20 40 --out /tmp/t7_disk.json
"""
import argparse
import json

import meep as mp


def run(res, n=3.4, r=1.0, dpml=1.0):
    cell = mp.Vector3(2 * (r + dpml + 1), 2 * (r + dpml + 1), 0)
    si = mp.Medium(index=n)
    geometry = [mp.Cylinder(radius=r, height=mp.inf, material=si)]
    pt = mp.Vector3(0.5, 0.0)
    src = mp.Source(mp.GaussianSource(0.3, fwidth=0.25), component=mp.Ez, center=pt)
    sim = mp.Simulation(cell_size=cell, resolution=res, geometry=geometry,
                        sources=[src], boundary_layers=[mp.PML(dpml)])
    h = mp.Harminv(mp.Ez, pt, 0.10, 0.60)
    sim.run(mp.after_sources(h),
            until_after_sources=mp.stop_when_fields_decayed(20, mp.Ez, pt, 1e-7))
    modes = sorted([{"freq": m.freq, "decay": m.decay, "Q": m.Q} for m in h.modes],
                   key=lambda x: x["freq"])
    return {"resolution": res, "modes": modes}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--res", type=int, nargs="+", default=[10, 20, 40])
    ap.add_argument("--out", default="/tmp/t7_disk.json")
    args = ap.parse_args()

    out = []
    for res in args.res:
        print(f"--- resolution {res} ---", flush=True)
        d = run(res)
        out.append(d)
        for m in d["modes"]:
            print(f"   f={m['freq']:.5f}  Q={m['Q']:.1f}", flush=True)
    with open(args.out, "w") as f:
        json.dump(out, f, indent=1)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
