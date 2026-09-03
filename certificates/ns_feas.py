import os as _os; _ROOT = _os.path.normpath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # release: paths relative to repo root
"""General-N NS saturation feasibility (v2, corrected gauge).

Goal: for every N >= 4, an explicit NS model with S = 2*sqrt(2) and per-row TV =
F* = (sqrt2-1)/4, built on the round-robin pattern Q_N (state i has miss-type
M_{i mod 4}; all +2 patterns).

Any model saturating S = 2 + 8F* must (Lemma 1 / PROOFS.md) have S0 = 2 and each
row ab bang-bang at TV = F*: p[ab,i] = 1/N + x[ab,i], x <= 0 on the wrong-sign
class r(ab), x >= 0 elsewhere, sum_W x = -F*, sum_C x = +F*.

CAVEAT handled here: no-signalling is NOT gauge-invariant.  The products
q_ab = u[a]v[b] fix S but not the single-party marginals that NS constrains.
We therefore use the GAUGE of the exact N=4 NS witness (ns_solN4_{u,v}.npy),
repeated per class:
    U0 = (+,+,-,-)   U1 = (-,+,+,-)
    V0 = (-,+,-,-)   V1 = (+,-,-,-)
With the symmetric ansatz x[ab,i] depending only on i mod 4 (16 unknowns), the
row constraints and all four NS equations are LINEAR -> feasibility LP.

NS equations (block index = concatenated ab):
    a=0: sum p[00] U0 = sum p[01] U0   (blocks 0,1)
    a=1: sum p[10] U1 = sum p[11] U1   (blocks 2,3)
    b=0: sum p[00] V0 = sum p[10] V0   (blocks 0,2)
    b=1: sum p[01] V1 = sum p[11] V1   (blocks 1,3)
(the uniform 1/N parts cancel in each equation).
"""
import numpy as np
from scipy.optimize import linprog

Fstar = (np.sqrt(2.0) - 1.0) / 4.0
U0 = np.array([ 1.0,  1.0, -1.0, -1.0])   # u[0,j] per class
U1 = np.array([-1.0,  1.0,  1.0, -1.0])   # u[1,j]
V0 = np.array([-1.0,  1.0, -1.0, -1.0])   # v[0,j]
V1 = np.array([ 1.0, -1.0, -1.0, -1.0])   # v[1,j]
MISS = {0: 0, 1: 1, 2: 2, 3: 3}           # row ab (00,01,10,11) -> wrong-sign class

def nclasses(N):
    return np.array([(N - j + 3) // 4 for j in range(4)], dtype=float)

def idx(ab, j):
    return ab * 4 + j

allok = True
for N in range(4, 25):
    n = nclasses(N)
    d = Fstar
    c = np.zeros(16)
    A_eq, b_eq = [], []
    bounds = []
    for ab in range(4):
        r = MISS[ab]
        row = np.zeros(16); row[idx(ab, r)] = n[r]; A_eq.append(row); b_eq.append(-d)
        row = np.zeros(16)
        for j in range(4):
            if j != r:
                row[idx(ab, j)] = n[j]
        A_eq.append(row); b_eq.append(d)
        for j in range(4):
            bounds.append((-1.0 / N, 0.0) if j == r else (0.0, 1.0 - 1.0 / N))
    # four NS equations with the witness gauge
    for (w, abL, abR) in [(U0, 0, 1), (U1, 2, 3), (V0, 0, 2), (V1, 1, 3)]:
        row = np.zeros(16)
        for j in range(4):
            row[idx(abL, j)] = n[j] * w[j]
            row[idx(abR, j)] = -n[j] * w[j]
        A_eq.append(row); b_eq.append(0.0)
    res = linprog(c, A_eq=np.array(A_eq), b_eq=np.array(b_eq), bounds=bounds, method='highs')
    if not res.success:
        allok = False
        print(f"N={N:>2}: INFEASIBLE under symmetric ansatz + witness gauge")
        continue
    x = res.x.reshape(4, 4)
    # build and verify the full model explicitly (no LP trust)
    cls = np.arange(N) % 4
    u = np.vstack([U0[cls], U1[cls]])
    v = np.vstack([V0[cls], V1[cls]])
    P = np.zeros((4, N))
    for ab in range(4):
        P[ab] = 1.0 / N + x[ab, cls]
    E = np.array([(P[ab] * u[a] * v[b]).sum() for ab, (a, b) in enumerate([(0, 0), (0, 1), (1, 0), (1, 1)])])
    S = E[0] + E[1] + E[2] - E[3]
    tv = np.array([0.5 * np.abs(P[ab] - 1.0 / N).sum() for ab in range(4)])
    ns_res = [
        abs((P[0] * u[0]).sum() - (P[1] * u[0]).sum()),
        abs((P[2] * u[1]).sum() - (P[3] * u[1]).sum()),
        abs((P[0] * v[0]).sum() - (P[2] * v[0]).sum()),
        abs((P[1] * v[1]).sum() - (P[3] * v[1]).sum()),
    ]
    ok = (abs(S - 2 * np.sqrt(2)) < 1e-9 and tv.max() <= Fstar + 1e-9
          and max(ns_res) < 1e-9 and P.min() > -1e-12)
    allok &= ok
    print(f"N={N:>2}: {'PASS' if ok else 'CHECK'}  S={S:.12f}  tv_max={tv.max():.12f} "
          f"(F*={Fstar:.12f})  ns_max={max(ns_res):.2e}  p_min={P.min():.3e}")

print("\nALL PASS" if allok else "\nSOME FAILED")
