#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""Compute WarpX-style field checksums (sum |field|) from an AMReX plotfile."""
import sys, json
import numpy as np
import yt
yt.funcs.mylog.setLevel(50)

plotfile = sys.argv[1]
ds = yt.load(plotfile)
grid = ds.covering_grid(level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
out = {}
for f in ["Ex", "Ey", "Ez", "Bx", "By", "Bz"]:
    if f in ds.field_list or ("boxlib", f) in ds.field_list:
        v = np.asarray(grid[f].v)
        out[f] = float(np.sum(np.abs(v)))
print(json.dumps({"lev=0": out}, indent=2))
print("time =", float(ds.current_time.to_value()))
