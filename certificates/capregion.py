import os as _os; _ROOT = _os.path.normpath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # release: paths relative to repo root
"""TASK 1: LP-vs-exact cross-check at bend points + cap points by explicit construction.

Writes capregion.out. Context (date-independent):
  - Signal-allowed problem, uniform source rho=1/N, per-row TV(p(.|ab)||rho) <= F.
  - Exact curve N=5 (PROOFS.md Theorem 3): 2+8F on [0,1/5]; 16/5+2F on [1/5,2/5]; then 4.
    => S_max(5,0.2)=3.6, S_max(5,0.25)=3.7, S_max(5,0.3)=3.8  (interior bend region).
  - Exact curve N=6: 2+8F on [0,1/6]; 8/3+4F on [1/6,1/3]; then 4.
    => S_max(6,0.2)=8/3+0.8=3.4666667 (OPTIONAL item).
  - Cap points (Theorem 2): S_max(N,F)=4 at F=ceil(N/4)/N: N=5 -> F=2/5=0.4, N=6 -> F=1/3.
    Proved by explicit construction (no LP): round-robin miss patterns Q_N; each row ab is
    supported uniformly on the class {lambda : q_ab(lambda) = sigma_ab}; then E[ab]=sigma_ab
    for all ab (S=4) and per-row TV = n_ab/N <= F. Trivial S<=4 gives equality.
  - LP route: phase2.exact_Smax(N,F) enumerates all 2^(4N) response vertices, one HiGHS LP each.
  - PASS tolerance for LP values: |LP - exact| <= 1e-5 (HiGHS tol ~1e-9 expected).
"""
import sys, time
import numpy as np

sys.path.insert(0, _ROOT + "/certificates")
from phase2 import exact_Smax

D = _ROOT + "/logs/"
OUT = open(D + "capregion.out", "w")
NPROC = 96  # <=120 processes total (shared machine)

def p(*a):
    print(*a, file=OUT, flush=True)

p("=== TASK 1: cap/bend region cross-check ===")
p(f"LP route: phase2.exact_Smax, Pool(nproc={NPROC}), HiGHS; PASS tol |diff|<=1e-5")
p("")

# ---------- Part A: LP vs exact at N=5 bend points ----------
p("--- Part A: N=5 LP evaluations at bend-region F (exact: 16/5+2F on [1/5,2/5]) ---")
t_n5 = []
ok_all = True
for F in [0.2, 0.25, 0.3]:
    t0 = time.time()
    S = exact_Smax(5, float(F), nproc=NPROC)
    dt = time.time() - t0
    t_n5.append(dt)
    exact = 16.0 / 5.0 + 2.0 * F
    d = abs(S - exact)
    ok = d <= 1e-5
    ok_all = ok_all and ok
    p(f"F={F:.2f}   LP Smax={S:.9f}   exact={exact:.9f}   |diff|={d:.3e}   "
      f"{'PASS' if ok else 'FAIL'}   (runtime {dt:.1f}s)")
p("")

# ---------- Part B: cap points by explicit construction (no LP) ----------
p("--- Part B: cap points by explicit construction (S=4, per-row TV = n_ab/N <= F) ---")
SIGMA = np.array([[1.0, 1.0], [1.0, -1.0]])
AB = [(a, b) for a in range(2) for b in range(2)]
# miss-type patterns q=(q00,q01,q10,q11), each with C(q)=+2 (PROOFS.md section 1)
M_A = np.array([-1., 1., 1., -1.])   # misses A: q00=+1
M_B = np.array([1., -1., 1., -1.])   # misses B: q01=+1
M_C = np.array([1., 1., -1., -1.])   # misses C: q10=+1
M_D = np.array([1., 1., 1., 1.])     # misses D: q11=-1
SEQ = [M_A, M_B, M_C, M_D]

def build_cap(N, F):
    """Q_N = first N of round-robin sequence; p[ab] uniform on correct-sign class."""
    Q = np.stack([SEQ[i % 4] for i in range(N)])          # (N,4) q patterns
    # realizability (PROOFS.md Theorem 1 (iii)): u0=+1, v0=q00, v1=q01, u1=q00*q10
    u = np.zeros((2, N)); v = np.zeros((2, N))
    u[0] = 1.0; v[0] = Q[:, 0]; v[1] = Q[:, 1]; u[1] = Q[:, 0] * Q[:, 2]
    pmat = np.zeros((4, N))
    n_wrong = np.zeros(4)
    for ab in range(4):
        a, b = AB[ab]
        correct = (Q[:, ab] == SIGMA[a, b])               # q_ab = sigma_ab
        n_wrong[ab] = int((~correct).sum())
        pmat[ab, correct] = 1.0 / correct.sum()           # uniform on correct class
    return u, v, pmat, n_wrong

for N, F in [(5, 0.4), (6, 1.0 / 3.0)]:
    u, v, pmat, n_wrong = build_cap(N, F)
    S = 0.0; tvs = []
    for ab in range(4):
        a, b = AB[ab]
        E = float(pmat[ab] @ (u[a] * v[b]))
        S += SIGMA[a, b] * E
        tvs.append(0.5 * float(np.abs(pmat[ab] - 1.0 / N).sum()))
    rowsum_ok = all(abs(float(pmat[ab].sum()) - 1.0) < 1e-12 for ab in range(4))
    tv_ok = all(t <= F + 1e-12 for t in tvs)
    s_ok = abs(S - 4.0) < 1e-12
    ok = rowsum_ok and tv_ok and s_ok
    ok_all = ok_all and ok
    p(f"N={N} F={F:.6f}: S={S:.12f}  per-row TV={[f'{t:.9f}' for t in tvs]}  "
      f"n_ab/N={[f'{w}/{N}' for w in n_wrong]}  rows-sum-1={rowsum_ok}")
    p(f"   construction gives S=4 (=> S_max>=4) and S<=4 trivial => S_max({N},{F:.6f})=4 exactly: "
      f"{'PASS' if ok else 'FAIL'}")
p("")

# ---------- Part C (optional): N=6 at F=0.2 ----------
med5 = float(np.median(t_n5))
est6 = med5 * 16.0 * 1.2   # 16x vertices, ~20% larger LPs -> safety factor 1.2
p("--- Part C (optional): N=6 at F=0.2, exact = 8/3+4*0.2 = 3.4666667 ---")
p(f"median N=5 eval runtime = {med5:.1f}s -> estimated N=6 runtime ~ {est6/60.0:.1f} min "
  f"(threshold 30 min)")
if est6 <= 30 * 60:
    t0 = time.time()
    S6 = exact_Smax(6, 0.2, nproc=NPROC)
    dt6 = time.time() - t0
    exact6 = 8.0 / 3.0 + 4.0 * 0.2
    d6 = abs(S6 - exact6)
    ok6 = d6 <= 1e-5
    ok_all = ok_all and ok6
    p(f"F=0.2   LP Smax={S6:.9f}   exact={exact6:.9f}   |diff|={d6:.3e}   "
      f"{'PASS' if ok6 else 'FAIL'}   (runtime {dt6:.1f}s)")
else:
    p(f"skipped: estimated runtime {est6/60.0:.1f} min > 30 min")
p("")
p(f"OVERALL TASK 1: {'PASS' if ok_all else 'FAIL'}")
OUT.close()
print("capregion.out written; overall:", "PASS" if ok_all else "FAIL", flush=True)
