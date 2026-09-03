import os as _os; _ROOT = _os.path.normpath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # release: paths relative to repo root
"""Phase 2: sample S_max(F; N) curves for N=2..6 (exact) + alternating method for N>=7.

Exact: enumerate all 2^{4N} bang-bang vertices (u,v) in parallel, LP over p each.
Alternating: for large N, coordinate ascent (LP over p given u,v; per-coordinate
sign flips for u,v given p) to find good constructions (lower bounds on F*).
"""
import numpy as np
from scipy.optimize import linprog
from multiprocessing import Pool
import sys, time

SIGMA = np.array([[1.0, 1.0], [1.0, -1.0]])
AB = [(a, b) for a in range(2) for b in range(2)]

def lp_max_S(N, u, v, F0):
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
        s = lp_max_S(N, u, v, F0)
        if s > best:
            best = s
    return best

def exact_Smax(N, F0, nproc=120, quiet=True):
    n = 2 ** (4 * N)
    chunks = np.linspace(0, n, nproc + 1, dtype=np.int64)
    tasks = [(N, F0, int(chunks[i]), int(chunks[i + 1])) for i in range(nproc) if chunks[i] < chunks[i + 1]]
    t0 = time.time()
    with Pool(nproc) as pool:
        best = max(pool.map(_work, tasks))
    if not quiet:
        print(f"   N={N} F={F0:.5f}: Smax={best:.7f}  ({time.time()-t0:.0f}s)", flush=True)
    return best

def alternating(N, F0, iters=200, restarts=8, seed=0):
    """Coordinate ascent on (u,v) bang-bang, LP over p each step. Returns best S."""
    rng = np.random.default_rng(seed)
    best = -np.inf
    for r in range(restarts):
        u = rng.integers(0, 2, (2, N)).astype(float) * 2 - 1
        v = rng.integers(0, 2, (2, N)).astype(float) * 2 - 1
        p = None
        for it in range(iters):
            S, _ = lp_full(N, u, v, F0)
            if S is None:
                break
            p = _
            # optimal u_{a, lam} = sign( sum_b SIGMA[a,b] p[ab,lam] v[b,lam] )
            for a in range(2):
                grad = np.zeros(N)
                for b in range(2):
                    grad += SIGMA[a, b] * p[a * 2 + b] * v[b]
                u[a] = np.where(grad >= 0, 1.0, -1.0)
            for b in range(2):
                grad = np.zeros(N)
                for a in range(2):
                    grad += SIGMA[a, b] * p[a * 2 + b] * u[a]
                v[b] = np.where(grad >= 0, 1.0, -1.0)
            best = max(best, S)
    return best

def lp_full(N, u, v, F0):
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
    if not res.success:
        return None, None
    return -res.fun, res.x[:4 * N].reshape(4, N)

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "curve"
    TGT = 2 * np.sqrt(2.0)
    if mode == "curve":
        for N in [2, 3, 4, 5, 6]:
            Fs = np.linspace(0.0, (2 * np.sqrt(2) - 2) / (2.0 * N) * 1.02, 11)
            print(f"=== N={N} curve ===", flush=True)
            for F0 in Fs:
                S = exact_Smax(N, float(F0), quiet=True)
                print(f"N={N}  F={F0:.6f}  Smax={S:.7f}  2+2NF={2+2*N*F0:.7f}  2+8F={2+8*F0:.7f}", flush=True)
    elif mode == "alt":
        N = int(sys.argv[2])
        for k in range(0, 13):
            F0 = (2 * np.sqrt(2) - 2) / (2.0 * N) * (1 - 1e-3) * (1 + 0.02 * k)
            S = alternating(N, F0, restarts=6)
            print(f"N={N}  F={F0:.6f}  S_alt={S:.7f}  (target {TGT:.7f})", flush=True)
    elif mode == "Fstar":
        # bisection on F for exact F* for given N using exact_Smax
        N = int(sys.argv[2])
        lo, hi = 0.0, 1.0
        for _ in range(16):
            mid = 0.5 * (lo + hi)
            S = exact_Smax(N, mid)
            if S >= TGT - 1e-6:
                hi = mid
            else:
                lo = mid
        ref = (np.sqrt(2) - 1) / min(N, 4)
        print(f"F*_N{N} ~= {hi:.8f}  ((sqrt2-1)/min(N,4)={ref:.8f})", flush=True)
