import os as _os; _ROOT = _os.path.normpath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # release: paths relative to repo root
"""TASK 3 (re-verification half): re-verify ns_solN{3,4}_{p,u,v}.npy from the saved files ALONE.

Runs as a FRESH process (no shared state with ns_witness.py, no LP involved). Appends its
PASS/FAIL section to ns_witness.out. Checks, per N in {3,4}, at F = (sqrt(2)-1)/N:
  - E[ab] = sum_lam p[ab,lam]*u[a,lam]*v[b,lam]; S = E00+E01+E10-E11  (~= 2*sqrt(2))
  - per-row TV = 0.5*sum|p[ab]-1/N| <= F + 1e-9
  - row sums = 1 +/- 1e-9
  - all four NS equality residuals (phase3 ns_rows) < 1e-8:
      sum p[00]u[0] - sum p[01]u[0];  sum p[10]u[1] - sum p[11]u[1];
      sum p[00]v[0] - sum p[10]v[0];  sum p[01]v[1] - sum p[11]v[1]
"""
import sys
import numpy as np

WIT = _ROOT + "/witnesses/"
LOG = _ROOT + "/logs/"
SIGMA = np.array([[1.0, 1.0], [1.0, -1.0]])
AB = [(a, b) for a in range(2) for b in range(2)]
TGT = 2.0 * np.sqrt(2.0)

out = open(LOG + "ns_witness.out", "a")
def p(*a):
    print(*a, file=out, flush=True)

all_ok = True
for N in [3, 4]:
    F = (np.sqrt(2.0) - 1.0) / float(N)
    u = np.load(WIT + f"ns_solN{N}_u.npy")
    v = np.load(WIT + f"ns_solN{N}_v.npy")
    pmat = np.load(WIT + f"ns_solN{N}_p.npy")

    p(f"=== TASK 3 (re-verification, fresh process): N={N}, F=(sqrt2-1)/{N}={F:.15f} ===")
    p(f"loaded ns_solN{N}_{{u,v,p}}.npy: shapes u{u.shape} v{v.shape} p{pmat.shape}, "
      f"dtypes {u.dtype}/{v.dtype}/{pmat.dtype}")

    shape_ok = (u.shape == (2, N) and v.shape == (2, N) and pmat.shape == (4, N)
                and set(np.unique(u)) <= {-1.0, 1.0} and set(np.unique(v)) <= {-1.0, 1.0})
    p(f"[{'PASS' if shape_ok else 'FAIL'}] shapes/conventions: u(2,N), v(2,N), p(4,N), "
      f"entries +/-1")
    all_ok = all_ok and shape_ok

    S = 0.0
    for ab in range(4):
        a, b = AB[ab]
        E = float(pmat[ab] @ (u[a] * v[b]))
        S += SIGMA[a, b] * E
        tv = 0.5 * float(np.abs(pmat[ab] - 1.0 / N).sum())
        rsum = float(pmat[ab].sum())
        ok_tv = tv <= F + 1e-9
        ok_rs = abs(rsum - 1.0) <= 1e-9
        all_ok = all_ok and ok_tv and ok_rs
        p(f"[{'PASS' if ok_tv else 'FAIL'}] row {ab} TV={tv:.12f} <= F+1e-9={F + 1e-9:.12f}")
        p(f"[{'PASS' if ok_rs else 'FAIL'}] row {ab} sum={rsum:.15f} = 1 +/- 1e-9")

    nsr = [
        float(pmat[0] @ u[0]) - float(pmat[1] @ u[0]),
        float(pmat[2] @ u[1]) - float(pmat[3] @ u[1]),
        float(pmat[0] @ v[0]) - float(pmat[2] @ v[0]),
        float(pmat[1] @ v[1]) - float(pmat[3] @ v[1]),
    ]
    for i, q in enumerate(nsr):
        ok = abs(q) < 1e-8
        all_ok = all_ok and ok
        p(f"[{'PASS' if ok else 'FAIL'}] NS residual {i} = {q:+.3e} (|.| < 1e-8)")

    ok_S = abs(S - TGT) < 1e-5
    all_ok = all_ok and ok_S
    p(f"[{'PASS' if ok_S else 'FAIL'}] S = {S:.12f} ~= 2*sqrt(2) = {TGT:.12f} "
      f"(diff {S - TGT:+.3e}, tol 1e-5)")
    p(f"--- N={N}: {'PASS' if all_ok else 'FAIL'} ---")
    p("")

p(f"OVERALL TASK 3 RE-VERIFICATION: {'PASS' if all_ok else 'FAIL'}")
out.close()
print("ns_verify done; overall:", "PASS" if all_ok else "FAIL", flush=True)
