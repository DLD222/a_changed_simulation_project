# SPDX-License-Identifier: BSD-3-Clause
"""1D toy check of the time-averaged conduction current.

Yee 1D in z: Ez at cell centres (j), Dz = eps*Ez, rho and div D at nodes.
Without curl, the exponential sigma update is exact:

    Ez^{n+1} = e^{-a dt} Ez^n ,   a = sigma/eps

Charge continuity with the step-averaged current:

    rho^{n+1} = rho^n - dt * sigma * c * (div Ez)^n ,
    c = (1 - e^{-a dt}) / (a dt)

Then div D - rho should stay at its initial value (zero) for ever.
Also compares against using the instantaneous current (c=1) to show the O(dt) drift.
"""
import numpy as np

def run(c_mode, n=200, dz=0.05, eps=4 * 8.8541878188e-12, sig=0.01, dt=1e-13, nstep=200):
    a = sig / eps
    # initial Ez: a smooth bump, so div Ez is non-zero
    z = (np.arange(n) + 0.5) * dz
    Ez = np.exp(-((z - n * dz / 2) / (0.15 * n * dz)) ** 2).astype(float)

    def div(f):                      # node -> (f[j]-f[j-1])/dz
        return np.diff(f) / dz

    rho = eps * np.concatenate(([0.0], div(Ez)))   # rho0 = div D0 (node j=1..n-1)

    drift = []
    for _ in range(nstep):
        if c_mode == "avg":
            c = (1 - np.exp(-a * dt)) / (a * dt) if a * dt > 1e-300 else 1.0
        elif c_mode == "inst_n":
            c = 1.0
        else:  # instantaneous at n+1
            c = np.exp(-a * dt)
        rho[1:] -= dt * sig * c * div(Ez)
        Ez *= np.exp(-a * dt)
        Drho = eps * np.concatenate(([0.0], div(Ez)))
        drift.append(np.abs(Drho - rho).max() / (eps * np.abs(Ez).max() / dz))
    return np.array(drift)

for mode in ("avg", "inst_n", "inst_np1"):
    d = run(mode)
    print(f"{mode:9s}: drift/norm  first={d[0]:.3e}  last={d[-1]:.3e}  max={d.max():.3e}  "
          f"{'GROWS' if d[-1] > 10*d[0] else 'flat'}")

# analytic factor check
eps = 4 * 8.8541878188e-12; sig = 0.01; dt = 1e-13
a = sig / eps
c = (1 - np.exp(-a * dt)) / (a * dt)
print(f"\na*dt = {a*dt:.3e}   c = {c:.15f}   (1-c)/a*dt = {(1-c)/(a*dt):.3e}")
print(f"identity check: exp(-a dt) vs 1 - a dt c : {np.exp(-a*dt):.17g} vs {1-a*dt*c:.17g}")
