import os as _os; _ROOT = _os.path.normpath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # release: paths relative to repo root
"""Gap-fill checks for the no-signalling variant:
  (1) N=2 NS: S_max over a fine F grid (confirm S=2 for all F, no Bell violation)
  (2) N=5 NS: S_max at F=(sqrt2-1)/4 (confirm = 2*sqrt(2), closing N=3,4,5 -> all N>=3)
"""
import sys
import numpy as np
sys.path.insert(0, _ROOT + "/certificates")
from phase3 import exact_Smax_ns

TGT = 2 * np.sqrt(2.0)

print("=== NS N=2 grid ===", flush=True)
maxS = -np.inf
for F0 in np.linspace(0.0, 0.4, 9):
    S = exact_Smax_ns(2, float(F0))
    maxS = max(maxS, S)
    print(f"NS N=2  F={F0:.4f}  Smax={S:.8f}", flush=True)
print(f"NS N=2 max over grid: {maxS:.8f}  (claim: exactly 2 for all F)", flush=True)

print("\n=== NS N=5 at F*=(sqrt2-1)/4 ===", flush=True)
Fstar4 = (np.sqrt(2) - 1) / 4.0
S5 = exact_Smax_ns(5, Fstar4)
print(f"NS N=5  F={Fstar4:.8f}  Smax={S5:.8f}  target={TGT:.8f}", flush=True)
print(f"NS N=5  F=(sqrt2-1)/5={ (np.sqrt(2)-1)/5:.8f}  Smax={exact_Smax_ns(5,(np.sqrt(2)-1)/5):.8f}  (below target expected)", flush=True)
