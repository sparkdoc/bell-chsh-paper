import os as _os; _ROOT = _os.path.normpath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # release: paths relative to repo root
"""Exact certificates for the no-signalling headline (Theorem 5, PROOFS.md).

Part A: general-N >= 4 closed-form NS construction, verified EXACTLY.
  Patterns: round-robin Q_N — state i has miss-type M_{i mod 4}, where
    M_0 = (-,+,+,-)  M_1 = (+,-,+,-)  M_2 = (+,+,-,-)  M_3 = (+,+,+,+)
  (all +2 patterns; class j = i mod 4, count n_j).
  Gauge (per class): U0=(+,+,-,-)  U1=(-,+,+,-)  V0=(-,+,-,-)  V1=(+,-,-,-)
  Rows: p[ab,i] = 1/N + x[ab, i mod 4], d = F* = (sqrt2-1)/4, with
    row 00: x[0,0] = -d/n_0        x[0,1] = +d/n_1
    row 01: x[1,0] = +d/n_0        x[1,1] = -d/n_1
    row 10: x[2,1] = +d/n_1        x[2,2] = -d/n_2
    row 11: x[3,0] = +d/n_0        x[3,3] = -d/n_3
  (all other x zero).  Verified exactly (Fraction, d formal):
    - row constraints  n_r x[ab,r] = -d,  sum_{j!=r} n_j x[ab,j] = +d
    - four NS identities (coefficients of d vanish identically in n_0..n_3)
    - bound  min_j n_j >= F* N   (case analysis on N mod 4; checked to N=999)
  Hence S = S0 + 8d = 2 + 8F* = 2 sqrt(2), TV_ab = d for all ab, NS holds.

Part B: N=3 closed-form witness in Q(sqrt2), d3 = (sqrt2-1)/3:
    row 00: [1/3-d3,   1/3+d3/2, 1/3+d3/2]
    row 01: [1/3+d3/2, 1/3-d3,   1/3+d3/2]
    row 10: [1/3-d3/2, 1/3+d3/2, 1/3     ]
    row 11: [1/3+d3/2, 1/3+d3/2, 1/3-d3  ]
  with the response/gauge of ns_solN3_{u,v}.npy; verified S = 2 sqrt(2) exactly
  in Q(sqrt2), per-row TV <= d3, NS residuals 0.

Part C: float cross-checks against the saved .npy witnesses (ns_solN3/N4).
"""
import math
from fractions import Fraction as F

# ---------------- Part A: exact, d formal ----------------
U0 = [1, 1, -1, -1]; U1 = [-1, 1, 1, -1]
V0 = [-1, 1, -1, -1]; V1 = [1, -1, -1, -1]

def coeffs(n):
    n0, n1, n2, n3 = (F(k) for k in n)
    c = [[F(0)] * 4 for _ in range(4)]
    c[0][0], c[0][1] = -1 / n0, 1 / n1
    c[1][0], c[1][1] = 1 / n0, -1 / n1
    c[2][1], c[2][2] = 1 / n1, -1 / n2
    c[3][0], c[3][3] = 1 / n0, -1 / n3
    return c

bad = []
for N in range(4, 65):
    n = [(N - j + 3) // 4 for j in range(4)]
    c = coeffs(n)
    for ab in range(4):                      # row constraints (coefficients of d)
        r = ab
        if F(n[r]) * c[ab][r] != -1:
            bad.append((N, ab, "R2"))
        if sum(F(n[j]) * c[ab][j] for j in range(4) if j != r) != 1:
            bad.append((N, ab, "R3"))
    for (w, abL, abR) in [(U0, 0, 1), (U1, 2, 3), (V0, 0, 2), (V1, 1, 3)]:
        if sum(F(n[j]) * w[j] * (c[abL][j] - c[abR][j]) for j in range(4)) != 0:
            bad.append((N, "NS", tuple(w)))
print(f"Part A: row constraints + 4 NS identities exact for N=4..64: "
      f"{'PASS (no violations)' if not bad else 'FAIL ' + str(bad[:5])}")
assert not bad, f"Part A violations: {bad[:5]}"

ok_bound = all(min((N - j + 3) // 4 for j in range(4)) * 4 >= (math.sqrt(2) - 1) * N - 1e-9
               for N in range(4, 1000))
print(f"Part A: min_j n_j >= F* N for all N=4..999: {ok_bound}")
assert ok_bound, "Part A: min_j n_j >= F* N failed for some N=4..999"

# ---------------- Part B: N=3 in Q(sqrt2) ----------------
def add(a, b): return (a[0] + b[0], a[1] + b[1])
def sub(a, b): return (a[0] - b[0], a[1] - b[1])
# (x + y*sqrt2)(u + v*sqrt2) = (xu + 2*yv) + (xv + yu)*sqrt2
def mul(a, b): return (a[0] * b[0] + 2 * a[1] * b[1], a[0] * b[1] + a[1] * b[0])
assert mul((F(0), F(1)), (F(0), F(1))) == (F(2), F(0)), "mul: sqrt2*sqrt2 != 2"
Z = (F(0), F(0)); O = (F(1), F(0))
d3 = (F(-1, 3), F(1, 3))
one3 = (F(1, 3), F(0))
hd3 = (d3[0] / 2, d3[1] / 2)
p3 = [
    [sub(one3, d3), add(one3, hd3), add(one3, hd3)],
    [add(one3, hd3), sub(one3, d3), add(one3, hd3)],
    [sub(one3, hd3), add(one3, hd3), one3],
    [add(one3, hd3), add(one3, hd3), sub(one3, d3)],
]

import numpy as np
u3 = np.load(_ROOT + "/witnesses/ns_solN3_u.npy"); v3 = np.load(_ROOT + "/witnesses/ns_solN3_v.npy")
AB = [(0, 0), (0, 1), (1, 0), (1, 1)]
Sx = Z
for ab, (a, b) in enumerate(AB):
    E = Z
    for i in range(3):
        E = add(E, mul(p3[ab][i], (F(u3[a, i] * v3[b, i]), F(0))))
    Sx = add(Sx, mul((F([1, 1, 1, -1][ab]), F(0)), E))
print(f"Part B: N=3 witness S == 2*sqrt(2) exactly in Q(sqrt2): {Sx == (F(0), F(2))}")
assert Sx == (F(0), F(2)), f"Part B: S != 2*sqrt2 exactly: {Sx}"

def sgn(a):
    x, y = a
    if x == 0 and y == 0: return 0
    if y == 0: return 1 if x > 0 else -1
    if x == 0: return 1 if y > 0 else -1
    lhs, rhs = x * x, 2 * y * y
    if lhs > rhs: return 1 if x > 0 else -1
    if lhs < rhs: return 1 if y > 0 else -1
    raise AssertionError("x^2 == 2y^2 with nonzero rationals — impossible")

def absz(a):
    s = sgn(a)
    return a if s >= 0 else (-a[0], -a[1])

# exact per-row TV in Q(sqrt2): TV_ab = (1/2) sum_i |p[ab,i] - 1/3|
tv3_exact = []
for ab in range(4):
    acc = Z
    for i in range(3):
        acc = add(acc, absz(sub(p3[ab][i], one3)))
    tv3_exact.append((acc[0] / 2, acc[1] / 2))
exp_tv = [d3, d3, hd3, d3]     # (d3, d3, d3/2, d3) — the paper's printed values
print(f"Part B: per-row TV exact in Q(sqrt2): {tv3_exact} vs expected {exp_tv}: "
      f"{tv3_exact == exp_tv}")
assert tv3_exact == exp_tv, f"Part B: per-row TV != (d3,d3,d3/2,d3) exactly: {tv3_exact}"

# exact NS residuals in Q(sqrt2): all four must vanish identically
def mrow(ab, w, i): return mul(p3[ab][i], (F(w[i]), F(0)))
nsr = [
    sub(add(add(mrow(0, u3[0], 0), mrow(0, u3[0], 1)), mrow(0, u3[0], 2)),
        add(add(mrow(1, u3[0], 0), mrow(1, u3[0], 1)), mrow(1, u3[0], 2))),
    sub(add(add(mrow(2, u3[1], 0), mrow(2, u3[1], 1)), mrow(2, u3[1], 2)),
        add(add(mrow(3, u3[1], 0), mrow(3, u3[1], 1)), mrow(3, u3[1], 2))),
    sub(add(add(mrow(0, v3[0], 0), mrow(0, v3[0], 1)), mrow(0, v3[0], 2)),
        add(add(mrow(2, v3[0], 0), mrow(2, v3[0], 1)), mrow(2, v3[0], 2))),
    sub(add(add(mrow(1, v3[1], 0), mrow(1, v3[1], 1)), mrow(1, v3[1], 2)),
        add(add(mrow(3, v3[1], 0), mrow(3, v3[1], 1)), mrow(3, v3[1], 2))),
]
print(f"Part B: all four NS residuals exactly 0 in Q(sqrt2): {nsr == [Z, Z, Z, Z]}")
assert nsr == [Z, Z, Z, Z], f"Part B: NS residual nonzero: {nsr}"

def fval(t): return float(t[0]) + float(t[1]) * math.sqrt(2)
tv3 = [max(abs(fval(p3[ab][i]) - 1 / 3) for i in range(3)) for ab in range(4)]
print(f"Part B: per-row TV <= d3 (float cross-check): {[f'{t:.15f}' for t in tv3]} vs d3={fval(d3):.15f}: "
      f"{all(t <= fval(d3) + 1e-15 for t in tv3)}")

# ---------------- Part C: float cross-check vs saved witnesses ----------------
pn3 = np.load(_ROOT + "/witnesses/ns_solN3_p.npy")
dev3 = max(abs(fval(p3[ab][i]) - pn3[ab, i]) for ab in range(4) for i in range(3))
print(f"Part C: closed-form N=3 vs saved ns_solN3_p.npy max deviation: {dev3:.3e}")
assert dev3 < 1e-4, f"Part C: closed form drifts from saved witness by {dev3}"

# general-N float spot check of the construction (independent re-evaluation)
d = (math.sqrt(2) - 1) / 4
for N in [5, 7, 10, 24]:
    n = [(N - j + 3) // 4 for j in range(4)]
    cls = np.arange(N) % 4
    u = np.vstack([np.array(U0)[cls], np.array(U1)[cls]])
    v = np.vstack([np.array(V0)[cls], np.array(V1)[cls]])
    P = np.zeros((4, N))
    for ab in range(4):
        P[ab] = 1.0 / N + d * np.array(coeffs(n)[ab])[cls]
    E = [sum(P[ab, i] * u[a, i] * v[b, i] for i in range(N)) for ab, (a, b) in enumerate(AB)]
    S = E[0] + E[1] + E[2] - E[3]
    tv = [0.5 * np.abs(P[ab] - 1.0 / N).sum() for ab in range(4)]
    nsr = [abs((P[0] * u[0]).sum() - (P[1] * u[0]).sum()),
           abs((P[2] * u[1]).sum() - (P[3] * u[1]).sum()),
           abs((P[0] * v[0]).sum() - (P[2] * v[0]).sum()),
           abs((P[1] * v[1]).sum() - (P[3] * v[1]).sum())]
    ok = abs(S - 2 * math.sqrt(2)) < 1e-12 and max(tv) <= d + 1e-12 and max(nsr) < 1e-12
    print(f"Part C: N={N:>2}: {'PASS' if ok else 'FAIL'}  S={S:.12f}  tv_max={max(tv):.12f}  ns_max={max(nsr):.2e}")
    assert ok, f"Part C: N={N} construction check failed (S={S}, tv={max(tv)}, ns={max(nsr)})"
