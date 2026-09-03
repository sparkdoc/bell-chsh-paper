import os as _os; _ROOT = _os.path.normpath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # release: paths relative to repo root
"""Phase 1 (fixed): exact minimum-fine-tuning CHSH for small N.

Bell scenario. Settings a,b in {0,1}, outcomes +-1, N hidden states.
  p[ab, lam] >= 0, sum_lam p[ab,lam] = 1        joint weights p(lam|a,b)
  u[a,lam], v[b,lam] in [-1,1]                  stochastic responses
  E[ab] = sum_lam p[ab,lam] u[a,lam] v[b,lam]
  S = E00 + E01 + E10 - E11                     CHSH  (classic |S|<=2, QM 2*sqrt(2))
  F = max_ab TV( p(.|a,b) || rho ),  rho uniform,  TV = 1/2 sum |p-rho|

Exactness: for fixed (u,v), S linear in p over a polytope -> LP.
At fixed p, S is multi-affine in the box variables u,v -> max at bang-bang
vertices {+/-1}^{2N x 2}. Enumerating all 2^{4N} vertices + LP is EXACT.

Also computes Hall-style metric: max_ab fraction of lam where p(.|a,b) != rho.
"""
import numpy as np
from scipy.optimize import linprog

SIGMA = np.array([[1, 1], [1, -1]])          # S = sum_ab SIGMA[a,b] E[ab]
AB = [(a, b) for a in range(2) for b in range(2)]

def lp_max_S(N, u, v, F0):
    """Exact max S over p at fine-tuning <= F0 (uniform source rho)."""
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
        return -np.inf, None
    return -res.fun, res.x[:4 * N].reshape(4, N)

def exact_optimum(N, F0, verbose=False):
    best, best_uv, best_p = -np.inf, None, None
    for m in range(2 ** (4 * N)):
        bits = ((m >> np.arange(4 * N)) & 1)
        u = np.where(bits[:2 * N].reshape(2, N) == 0, -1.0, 1.0)
        v = np.where(bits[2 * N:].reshape(2, N) == 0, -1.0, 1.0)
        S, p = lp_max_S(N, u, v, F0)
        if S > best:
            best, best_uv, best_p = S, (u.copy(), v.copy()), p
        if verbose and m % 4096 == 0:
            print(f"  N={N} F0={F0:.5f} m={m}/{2**(4*N)} best={best:.6f}", flush=True)
    return best, best_uv, best_p

def frac_finetuning(N, p):
    f = 0.0
    for ab in range(4):
        f = max(f, np.mean(np.abs(p[ab] - 1.0 / N) > 1e-9))
    return f

def min_F0_for_target(N, target, lo, hi, tol=1e-5, verbose=False):
    Slo, _, _ = exact_optimum(N, lo)
    Shigh, _, _ = exact_optimum(N, hi)
    assert Slo <= 2.0 + 1e-6, f"Bell bound violated at F=0: {Slo}"
    assert Shigh <= 4.0 + 1e-6, f"trivial bound violated: {Shigh}"
    print(f"N={N}: S_max(F=0)={Slo:.6f}  S_max(F=1)={Shigh:.6f}", flush=True)
    if Shigh < target - 1e-9:
        return None
    if Slo >= target - 1e-9:
        return 0.0
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        S, _, _ = exact_optimum(N, mid, verbose=verbose)
        if S >= target - 1e-9:
            hi = mid
        else:
            lo = mid
    return hi

if __name__ == "__main__":
    import sys
    Nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    TGT = 2 * np.sqrt(2)
    for N in range(2, Nmax + 1):
        print(f"=== N={N} ===", flush=True)
        fstar = min_F0_for_target(N, TGT, 0.0, 1.0)
        if fstar is None:
            print(f"N={N}: unreachable even at F=1", flush=True); continue
        S, (u, v), p = exact_optimum(N, fstar)
        # independent re-verification from the returned data
        Scheck = 0.0
        tv = 0.0
        for ab in range(4):
            a, b = AB[ab]
            E = float(p[ab] @ (u[a] * v[b]))
            Scheck += SIGMA[a, b] * E
            tv = max(tv, 0.5 * float(np.abs(p[ab] - 1.0 / N).sum()))
        assert abs(Scheck - S) < 1e-7 and abs(Scheck - TGT) < 1e-3, (S, Scheck)
        print(f"N={N}: F*~={fstar:.6f}  S={S:.8f} (recheck {Scheck:.8f})  "
              f"actual-TV={tv:.6f}  target={TGT:.8f}  "
              f"frac-lam-tuned={frac_finetuning(N, p):.5f}", flush=True)
        d = _ROOT + "/witnesses/"
        np.save(d + f"solN{N}_u.npy", u); np.save(d + f"solN{N}_v.npy", v)
        np.save(d + f"solN{N}_p.npy", p)
        print("u=\n" + np.array2string(u, precision=2) +
              "\nv=\n" + np.array2string(v, precision=2) +
              "\np[ab,lam]=\n" + np.array2string(p, precision=5), flush=True)
