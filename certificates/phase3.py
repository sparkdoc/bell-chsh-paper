import os as _os; _ROOT = _os.path.normpath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # release: paths relative to repo root
"""Phase 3: (A) no-signalling variant: add NS constraints to the LP.
       (B) N=5 bisection to pin F* for the signalling-allowed problem.
       (C) fine-grained N=5 curve points near the predicted break.

NS condition (Hall sense): p(b|a) independent of a  <=>
  for all a:  sum_lam p[ab,lam]*u[a,lam] is the same for b=0,1
(and symmetric for v).  Note: both sides must be imposed EXPLICITLY — one side does NOT imply the
other (under A-side-only NS even N=2 reaches S = 2+4F; see PROOFS.md §9 side finding). All four NS
rows are added in lp_max_S_ns.
"""
import numpy as np
from scipy.optimize import linprog
from multiprocessing import Pool
import sys, time

SIGMA = np.array([[1.0, 1.0], [1.0, -1.0]])
AB = [(a, b) for a in range(2) for b in range(2)]

def lp_max_S_ns(N, u, v, F0):
    rho = np.full(N, 1.0 / N)
    c = np.zeros(4 * N)
    for ab in range(4):
        a, b = AB[ab]
        c[ab * N:(ab + 1) * N] = SIGMA[a, b] * u[a] * v[b]
    n = 12 * N
    c_full = np.zeros(n); c_full[:4 * N] = c
    nrows = 4 + 4 * N + 4        # row-sums, p-s+t=rho, four NS rows
    A_eq = np.zeros((nrows, n)); b_eq = np.zeros(nrows)
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
    r = 4 + 4 * N
    # AB blocks: 0=(0,0) 1=(0,1) 2=(1,0) 3=(1,1)
    # A-side NS: p(A|a,b) indep of b  ->  sum_lam p[ab,lam]u[a,lam] equal in b, for each a
    # B-side NS: p(B|a,b) indep of a  ->  sum_lam p[ab,lam]v[b,lam] equal in a, for each b
    ns_rows = [
        (u[0], 0, 1),   # a=0: sum p[00]u[0] = sum p[01]u[0]
        (u[1], 2, 3),   # a=1: sum p[10]u[1] = sum p[11]u[1]
        (v[0], 0, 2),   # b=0: sum p[00]v[0] = sum p[10]v[0]
        (v[1], 1, 3),   # b=1: sum p[01]v[1] = sum p[11]v[1]
    ]
    for (w, abL, abR) in ns_rows:
        A_eq[r, abL * N:(abL + 1) * N] = w
        A_eq[r, abR * N:(abR + 1) * N] = -w
        b_eq[r] = 0.0
        r += 1
    res = linprog(-c_full, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                  bounds=[(0, 1)] * (4 * N) + [(0, None)] * (8 * N), method='highs')
    return -res.fun if res.success else -np.inf

def vertex_sv(m, N):
    bits = ((m >> np.arange(4 * N)) & 1)
    u = np.where(bits[:2 * N].reshape(2, N) == 0, -1.0, 1.0)
    v = np.where(bits[2 * N:].reshape(2, N) == 0, -1.0, 1.0)
    return u, v

def _work(args):
    N, F0, m0, m1 = args
    best = -np.inf
    for m in range(m0, m1):
        u, v = vertex_sv(m, N)
        s = lp_max_S_ns(N, u, v, F0)
        if s > best:
            best = s
    return best

def exact_Smax_ns(N, F0, nproc=120):
    n = 2 ** (4 * N)
    chunks = np.linspace(0, n, nproc + 1, dtype=np.int64)
    tasks = [(N, F0, int(chunks[i]), int(chunks[i + 1])) for i in range(nproc) if chunks[i] < chunks[i + 1]]
    with Pool(nproc) as pool:
        return max(pool.map(_work, tasks))

def bisect_Fstar_ns(N, TGT, lo=0.0, hi=0.5, iters=13):
    for it in range(iters):
        mid = 0.5 * (lo + hi)
        S = exact_Smax_ns(N, mid)
        print(f"  [NS N={N}] it{it} F={mid:.8f} Smax={S:.8f}  {'<-' :}", flush=True)
        if S >= TGT - 1e-6:
            hi = mid
        else:
            lo = mid
    return hi

if __name__ == "__main__":
    mode = sys.argv[1]
    TGT = 2 * np.sqrt(2.0)
    if mode == "ns_curve":
        # sample the exact NS curve for N=2..4 over a fine grid
        for N in [2, 3, 4]:
            Fs = np.linspace(0.10, 0.17, 8)
            print(f"=== NS N={N} curve ===", flush=True)
            for F0 in Fs:
                S = exact_Smax_ns(N, float(F0))
                print(f"NS N={N}  F={F0:.6f}  Smax={S:.8f}  2+6F={2+6*F0:.8f}", flush=True)
    elif mode == "ns_bisect":
        N = int(sys.argv[2])
        fstar = bisect_Fstar_ns(N, TGT)
        ref = (np.sqrt(2) - 1) / min(N, 4)
        print(f"F*_NS(N={N}) ~= {fstar:.8f}   (sqrt2-1)/min(N,4)={ref:.8f}  2sqrt2={TGT:.8f}", flush=True)
    elif mode == "n5_bisect":
        # signalling-allowed exact F* for N=5: 16 bisection evals of 2^20 vertices
        sys.path.insert(0, _ROOT + "/certificates")
        from phase2 import exact_Smax
        N = 5
        lo, hi = 0.0845, 0.140833  # hi is the F where Smax=2.8449957>target from curve
        for it in range(16):
            mid = 0.5 * (lo + hi)
            S = exact_Smax(N, mid)
            print(f"  [N=5] it{it} F={mid:.8f} Smax={S:.8f}", flush=True)
            if S >= TGT - 1e-6:
                hi = mid
            else:
                lo = mid
        print(f"F*_N5 ~= {hi:.8f}  (sqrt2-1)/5={((np.sqrt(2)-1)/5):.8f}  (3-2sqrt2)/3={((3-2*np.sqrt(2))/3):.8f}", flush=True)
    elif mode == "n5_fine":
        # exact N=5 single evaluations at fine F points near the break
        sys.path.insert(0, _ROOT + "/certificates")
        from phase2 import exact_Smax
        for F0 in [0.0855, 0.0860, 0.0865, 0.0870]:
            S = exact_Smax(5, F0)
            print(f"NSIG N=5  F={F0:.6f}  Smax={S:.8f}  target={TGT:.8f}", flush=True)
