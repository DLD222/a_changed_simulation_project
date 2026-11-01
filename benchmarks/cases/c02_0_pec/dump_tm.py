#!/usr/bin/env python3
"""TM-mode C2-0: use Ex (col 5) for incident/reflected pulses."""
import numpy as np
d = np.loadtxt("diags/reducedfiles/fp_inc.txt", skiprows=1)
t = d[:,1]; Ex = d[:,5]
print("t range: %.2f .. %.2f ns" % (t.min()*1e9, t.max()*1e9))
print("Ex max=%.4e  min=%.4e" % (Ex.max(), Ex.min()))
ip = np.argmax(Ex); ineg = np.argmin(Ex)
print(f"positive peak: t={t[ip]*1e9:.2f} ns  Ex={Ex[ip]:.4e}")
print(f"negative peak: t={t[ineg]*1e9:.2f} ns  Ex={Ex[ineg]:.4e}")
if ip != ineg and Ex[ip] != 0:
    print(f"ratio (neg/pos) = {Ex[ineg]/Ex[ip]:.4f}   (PEC expects -1)")
# coarse trace
print("\nEx(t):")
for k in range(0, len(t), 40):
    print(f"  t={t[k]*1e9:6.2f} ns  Ex={Ex[k]:+.4e}")
