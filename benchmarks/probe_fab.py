#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
import numpy as np, struct, sys

path = sys.argv[1]
raw = open(path, 'rb').read()
print("file size", len(raw))

# known first-cell values from Cell_H ASCII preview
target = -2.62785367770417452e+02

# search every byte offset for a matching float64
found = []
for off in range(0, 200):
    if off + 8 <= len(raw):
        v = struct.unpack('<d', raw[off:off+8])[0]
        if abs(v - target) < 1e-6 * abs(target):
            found.append(off)
print("offset(s) where first value found:", found)

# try interpreting the rest as doubles
for off in found:
    arr = np.frombuffer(raw[off:], dtype='<f8')
    print(f"offset {off}: n doubles = {len(arr)}")
    print("  first 6:", arr[:6])
    print("  target first 6: -262.785, 0, -4.2756, 0, -8.2435e-07, 0")
