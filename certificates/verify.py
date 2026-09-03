import os as _os; _ROOT = _os.path.normpath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # release: paths relative to repo root
"""Self-contained verification of the main result.

Run (fast, default):   python3 verify.py
Run (full incl. N=5,6): python3 verify.py --heavy

Default checks (exact, single-process, ~a few minutes):
  (A) N=4 construction: explicit (p,u,v) achieves S = 2*sqrt(2) at F = (sqrt2-1)/4,
      verified by DIRECT arithmetic (not LP), rows sum to 1, TV recomputed independently.
  (B) Bell bound: with F=0, the maximum achievable S is exactly 2 (N=2,3,4, exact enumeration).
  (C) Upper bound S_max(N,F) <= 2 + 8F and tightness S_max(N,F) = 2 + 2*min(N,4)*F
      at several F for N=2,3,4 (exact enumeration over all 2^{4N} vertices + LP over p).
  (D) Headline: F*(2*sqrt2) = (sqrt2-1)/min(N,4) for N=2,3,4 (exact bisection).

--heavy additionally (parallel, ~30-60 min):
  (E) S_max(5,F), S_max(6,F) == 2+8F at F in {0.05, 0.10}  and F*(2*sqrt2) for N=5,6.

Every number is recomputed from raw data at run time.
"""
import sys
import numpy as np
from scipy.optimize import linprog

SIGMA = np.array([[1., 1.], [1., -1.]])
AB = [(a, b) for a in range(2) for b in range(2)]
R2 = np.sqrt(2.0)
TGT = 2.0 * R2

def lp_max_S(N, u, v, F0):
    """Exact max S over p (row-stoch, per-row TV(p(.|ab)||1/N) <= F0) for fixed u,v."""
    rho = np.full(N, 1.0 / N)
    c = np.zeros(4 * N)
    for ab in range(4):
        a, b = AB[ab]
        c[ab * N:(ab + 1) * N] = SIGMA[a, b] * u[a] * v[b]
    n = 12 * N
    c_full = np.zeros(n); c_full[:4 * N] = c
    A_eq = np.zeros((4 + 4 * N, n)); b_eq = np.zeros(4 + 4 * N)
    A_ub = np.zeros((4, n)); b_ub = np.full(4, 2.0 * F0)
    for ab in range(4):
        A_eq[ab, ab * N:(ab + 1) * N] = 1.0; b_eq[ab] = 1.0
        base = 4 + ab * N
        for lam in range(N):
            A_eq[base + lam, ab * N + lam] = 1.0
            A_eq[base + lam, 4 * N + ab * N + lam] = -1.0
            A_eq[base + lam, 8 * N + ab * N + lam] = 1.0
            b_eq[base + lam] = rho[lam]
        A_ub[ab, 4 * N + ab * N:4 * N + (ab + 1) * N] = 1.0
        A_ub[ab, 8 * N + ab * N:8 * N + (ab + 1) * N] = 1.0
    res = linprog(-c_full, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                  bounds=[(0, 1)] * (4 * N) + [(0, None)] * (8 * N), method='highs')
    return -res.fun if res.success else -np.inf

def exact_Smax_small(N, F0):
    """Exhaustive over all 2^{4N} bang-bang (u,v); LP over p. Exact for N<=4 here."""
    best = -np.inf
    for m in range(2 ** (4 * N)):
        bits = ((m >> np.arange(4 * N)) & 1)
        u = np.where(bits[:2 * N].reshape(2, N) == 0, -1.0, 1.0)
        v = np.where(bits[2 * N:].reshape(2, N) == 0, -1.0, 1.0)
        best = max(best, lp_max_S(N, u, v, F0))
    return best

def S_and_TV(N, p, u, v):
    """Independent recompute of S and F from raw (p,u,v) — no LP involved."""
    S = 0.0; tv = 0.0; rows_ok = True
    for ab in range(4):
        a, b = AB[ab]
        S += SIGMA[a, b] * float(p[ab] @ (u[a] * v[b]))
        tv = max(tv, 0.5 * float(np.abs(p[ab] - 1.0 / N).sum()))
        rows_ok = rows_ok and abs(float(p[ab].sum()) - 1.0) < 1e-12
    return S, tv, rows_ok

ok = True
def check(name, cond, detail=""):
    global ok
    ok = ok and bool(cond)
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}  {detail}", flush=True)

print("=== (A) N=4 construction, DIRECT arithmetic ===")
Fstar4 = (R2 - 1.0) / 4.0
a0 = 0.25 - Fstar4
b0 = 0.25 + Fstar4
p4 = np.array([
    [0.25, a0, b0, 0.25],
    [a0, 0.25, b0, 0.25],
    [b0, 0.25, a0, 0.25],
    [b0, 0.25, 0.25, a0],
])
u4 = np.array([[1., -1., -1., -1.],
               [1.,  1.,  1., -1.]])
v4 = np.array([[1.,  1., -1., -1.],
               [-1., -1., -1., -1.]])
S4, F4, rows_ok = S_and_TV(4, p4, u4, v4)
check("rows of p sum to 1", rows_ok)
check("S = 2*sqrt(2)", abs(S4 - TGT) < 1e-12, f"S={S4:.12f}  target={TGT:.12f}")
check("F = (sqrt2-1)/4", abs(F4 - Fstar4) < 1e-12, f"F={F4:.12f}  target={Fstar4:.12f}")

print("\n=== (B) Bell bound: F=0 => S_max = 2 (exact enumeration) ===")
for N in [2, 3, 4]:
    S0 = exact_Smax_small(N, 0.0)
    check(f"N={N}  S_max(F=0)=2", abs(S0 - 2.0) < 1e-6, f"S={S0:.6f}")

print("\n=== (C) S_max(N,F) = 2 + 2*min(N,4)*F (exact enumeration, N=2,3,4) ===")
for N in [2, 3, 4]:
    for F in [0.02, 0.05, 0.08, 0.10]:
        S = exact_Smax_small(N, F)
        pred = 2.0 + 2.0 * min(N, 4) * F
        check(f"N={N} F={F:.2f}  S=2+2*min(N,4)*F", abs(S - pred) < 1e-5,
              f"S={S:.6f}  pred={pred:.6f}  (2+8F={2+8*F:.4f} upper bound)")

print("\n=== (D) headline: F*(2*sqrt2) = (sqrt2-1)/min(N,4), N=2,3,4 ===")
for N in [2, 3, 4]:
    lo, hi = 0.0, 0.30
    for _ in range(20):
        mid = 0.5 * (lo + hi)
        S = exact_Smax_small(N, mid)
        if S >= TGT - 1e-7:
            hi = mid
        else:
            lo = mid
    pred = (R2 - 1.0) / min(N, 4)
    check(f"N={N}  F*=(sqrt2-1)/min(N,4)", abs(hi - pred) < 1e-5,
          f"F*~={hi:.6f}  pred={pred:.6f}")

if "--heavy" in sys.argv:
    print("\n=== (E) heavy: N=5,6 confirm S_max=2+8F and F* (parallel exact) ===")
    sys.path.insert(0, _ROOT + "/certificates")
    from phase2 import exact_Smax
    for N in [5, 6]:
        for F in [0.05, 0.10]:
            S = exact_Smax(N, F)
            check(f"N={N} F={F:.2f}  S=2+8F", abs(S - (2 + 8 * F)) < 1e-5,
                  f"S={S:.6f}  2+8F={2+8*F:.6f}")
        lo, hi = 0.05, 0.30
        for _ in range(14):
            mid = 0.5 * (lo + hi)
            S = exact_Smax(N, mid)
            if S >= TGT - 1e-6:
                hi = mid
            else:
                lo = mid
        pred = (R2 - 1.0) / 4.0
        check(f"N={N}  F*=(sqrt2-1)/4", abs(hi - pred) < 5e-5,
              f"F*~={hi:.6f}  pred={pred:.6f}")

print("\n" + ("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED"))
raise SystemExit(0 if ok else 1)
