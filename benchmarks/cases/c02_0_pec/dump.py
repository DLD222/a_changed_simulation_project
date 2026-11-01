#!/usr/bin/env python3
import numpy as np
d = np.loadtxt("diags/reducedfiles/fp_inc.txt", skiprows=1)
t, Ey = d[:,1], d[:,6]
# print extrema of the signal
print("t: min=%.3e max=%.3e ns" % (t.min()*1e9, t.max()*1e9))
print("Ey: max=%.4e min=%.4e" % (Ey.max(), Ey.min()))
# locate the main positive and negative peaks
ip = np.argmax(Ey); ineg = np.argmin(Ey)
print(f"main positive peak: t={t[ip]*1e9:.3f} ns, Ey={Ey[ip]:.4e}")
print(f"main negative peak: t={t[ineg]*1e9:.3f} ns, Ey={Ey[ineg]:.4e}")
if ip != ineg:
    print(f"ratio (neg/pos) = {Ey[ineg]/Ey[ip]:.6f}  (expect -1)")
# print a coarse trace
print("\nEy(t) trace (every 10th):")
for k in range(0, len(t), 10):
    print(f"  t={t[k]*1e9:6.2f} ns  Ey={Ey[k]:+.4e}")
