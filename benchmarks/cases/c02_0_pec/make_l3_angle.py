#!/usr/bin/env python3
"""Generate an L3 input for an arbitrary incidence angle from inputs_2d_l3.

The probe positions are placed on the specularly reflected beam
k_r_hat = (-sin 2a, cos 2a), at arc length s = 3..8 m from the specular point:

    P(s) = s * k_r_hat

Usage:  python3 make_l3_angle.py <angle_deg> [out_file]
Writes inputs_2d_l3_<angle>deg by default.
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TEMPLATE = os.path.join(HERE, "inputs_2d_l3")
S_LIST = (3, 4, 5, 6, 7, 8)


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    angle = float(sys.argv[1])
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
        HERE, f"inputs_2d_l3_{int(round(angle))}deg")

    alpha = math.radians(angle)
    text = open(TEMPLATE).read()
    head = text.split("warpx.reduced_diags_names")[0]
    head = head.replace("my_constants.alpha = 0.5235987756     # 30 degrees",
                        f"my_constants.alpha = {alpha!r}     # {angle:g} degrees")

    names = ["fp_inc", "fp_shadow"] + [f"fp_ref{s}" for s in S_LIST]
    lines = [f"# Auto-generated from inputs_2d_l3 for alpha = {angle:g} deg "
             f"(see make_l3_angle.py).",
             f"warpx.reduced_diags_names = {' '.join(names)}",
             "fp_inc.type = FieldProbe", "fp_inc.intervals = 2",
             "fp_inc.probe_geometry = Point", "fp_inc.x_probe = 0.0",
             "fp_inc.z_probe = -5.0", ""]
    for s in S_LIST:
        x = -s * math.sin(2 * alpha)
        z = s * math.cos(2 * alpha)
        lines += [f"# Reflected beam, s = {s} m: k_r_hat = "
                  f"({-math.sin(2 * alpha):+.6f}, {math.cos(2 * alpha):+.6f})",
                  f"fp_ref{s}.type = FieldProbe", f"fp_ref{s}.intervals = 2",
                  f"fp_ref{s}.probe_geometry = Point",
                  f"fp_ref{s}.x_probe = {x!r}", f"fp_ref{s}.z_probe = {z!r}", ""]
    lines += ["# Inside the conductor: must stay zero.",
              "fp_shadow.type = FieldProbe", "fp_shadow.intervals = 2",
              "fp_shadow.probe_geometry = Point", "fp_shadow.x_probe = 0.0",
              "fp_shadow.z_probe = 4.0", ""]

    open(out, "w").write(head + "\n".join(lines))
    print(f"wrote {out} (alpha = {angle:g} deg)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
