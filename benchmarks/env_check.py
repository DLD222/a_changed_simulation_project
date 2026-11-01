# SPDX-License-Identifier: BSD-3-Clause
import importlib
for m in ["numpy", "yt", "openpmd_viewer", "scipy", "h5py", "adios2"]:
    try:
        mod = importlib.import_module(m)
        print(f"{m}: OK {getattr(mod, '__version__', '')}")
    except Exception as e:
        print(f"{m}: MISSING ({type(e).__name__})")
