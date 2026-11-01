#!/usr/bin/env python3
"""Minimal AMReX plotfile checksum: parse Header + Level_0/Cell_D_00000.
Computes sum|field| per component, comparable to WarpX Regression benchmark JSON."""
import numpy as np, sys, re, json

plotdir = sys.argv[1]
hdr = open(plotdir + "/Header").read().split("\n")
nvars = int(hdr[1])
varnames = [hdr[2 + k] for k in range(nvars)]
dim = int(hdr[2 + nvars])
# level 0 box appears later; find the line with '((0,0' pattern
boxline = [l for l in hdr if l.startswith("((")][0]
# e.g. ((0,0) (31,31) (0,0))
nums = re.findall(r"\((\d+),(\d+)\)", boxline)
lo = (int(nums[0][0]), int(nums[0][1]))
hi = (int(nums[1][0]), int(nums[1][1]))
nx = hi[0] - lo[0] + 1
ny = hi[1] - lo[1] + 1
ncells = nx * ny
print(f"plotfile={plotdir}")
print(f"vars={varnames} dim={dim} box={(lo,hi)} nx*ny={nx}x{ny}")

raw = open(plotdir + "/Level_0/Cell_D_00000", "rb").read()
# header is a text FAB descriptor ending with the ncomp integer + '\n'
# find the newline that terminates the header: data size must equal ncells*nvars*8
need = ncells * nvars * 8
off = len(raw) - need
assert off > 0, f"unexpected file size {len(raw)} vs need {need}"
print(f"data offset = {off} bytes (file {len(raw)}, header {off})")
arr = np.frombuffer(raw[off:off+need], dtype="<f8")
arr = arr.reshape(nvars, ncells)   # component-major
sums = np.sum(np.abs(arr), axis=1)
res = {}
for k, name in enumerate(varnames):
    res[name] = float(sums[k])
    print(f"  {name}: sum|.| = {sums[k]:.12e}")
print(json.dumps(res, indent=2))
