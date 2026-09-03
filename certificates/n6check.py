import os as _os; _ROOT = _os.path.normpath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # release: paths relative to repo root
import sys
sys.path.insert(0, _ROOT + "/certificates")
from phase2 import exact_Smax
import numpy as np
TGT = 2*np.sqrt(2.0)
# two-point check around (sqrt2-1)/4 = 0.10355339059327376
for F0 in [0.10300, (np.sqrt(2)-1)/4]:
    S = exact_Smax(6, F0)
    print(f"N=6  F={F0:.8f}  Smax={S:.8f}  target={TGT:.8f}", flush=True)
