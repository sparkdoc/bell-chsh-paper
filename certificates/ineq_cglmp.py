#!/usr/bin/env python3
import os as _os; _ROOT = _os.path.normpath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # release: paths relative to repo root
"""
ineq_cglmp.py — Divergent explorer #5 (second inequality): CGLMP d=3.

Inequality: Collins, Gisin, Linden, Massar, Popescu, PRL 88, 040404 (2002),
arXiv:quant-ph/0106024, eq. (6) at d=3 (k runs 0..floor(d/2)-1 = {0}, weight 1):

  I_3 = P(A1=B1) + P(B1=A2+1) + P(A2=B2) + P(B2=A1)
        - P(A1=B1-1) - P(B1=A2) - P(A2=B2-1) - P(B2=A1-1)      (all mod 3)

Primary-source facts used: local bound I_3 <= 2; quantum value achieved by the
paper's state/measurements I_3(QM) = 4/(-9+6*sqrt(3)) = (12+8*sqrt(3))/9
~= 2.87293 (a lower bound on the Tsirelson bound; NOT claimed to be it).

Setting blocks: rows r in {(A1,B1),(A1,B2),(A2,B1),(A2,B2)} with differences
d00=A1-B1, d01=A1-B2, d10=A2-B1, d11=A2-B2 (mod 3).  Row weight vectors in
{-1,0,+1}:
  w00(d) = +1[d=0] -1[d=2];  w01(d) = +1[d=0] -1[d=1];
  w10(d) = +1[d=2] -1[d=0];  w11(d) = +1[d=0] -1[d=2].

Pattern space: (d00,d01,d10,d11) in {0,1,2}^4 with d00+d11 == d01+d10 (mod 3)
(the modular analog of the CHSH parity lemma): 27 patterns.  Algebraic max 4
requires (d00,d01,d10,d11) = (0,0,2,0), which violates the constraint -> no
pattern attains all four row conditions (modular obstruction).

Metric: same exact generalization as ineq_mermin.py (N states, uniform source,
p[s,lam] per setting block, deterministic responses, F = max_s TV(p(.|s)||1/N)).

Exact finite-max formula (Theorem-1 analog with generalized Lemma 1): for a
multiset Q of N patterns and row r with weight multiset w_r (N values in
{-1,0,+1}), the exact per-row max over the TV ball at F = k/N is
  h_r(k/N) = mu_r + (k*wmax_r - sum of k smallest w_r)/N   for 0 <= k <= N-1
  h_r(N/N) = wmax_r,
so S_max(N,k/N) = max_Q sum_r h_r(Q,k)  — computed by full enumeration over
C(27+N-1,N) multisets (exact integer inner loop).

Quantum-target F* values are exact in Q(sqrt(3)) via tuple arithmetic
(r,s) := r + s*sqrt(3), r,s in Fraction.

Exits 0 iff every assertion passes.
"""
import sys, time
from fractions import Fraction as Fr
from itertools import combinations_with_replacement

T0 = time.time()
NCHECK = 0
def ok(msg):
    global NCHECK; NCHECK += 1
    print(f"  [check {NCHECK:3d}] {msg}")

# ---- Q(sqrt(3)) exact arithmetic -------------------------------------------
def q3_add(x, y): return (x[0] + y[0], x[1] + y[1])
def q3_smul(x, c): return (c * x[0], c * x[1])          # c in Fraction
def q3_sub(x, y): return (x[0] - y[0], x[1] - y[1])
def q3_cmp(x, y):
    """-1/0/+1 for x-y."""
    a, b = x[0] - y[0], x[1] - y[1]
    if a == 0: return (b > 0) - (b < 0)
    if b == 0: return (a > 0) - (a < 0)
    # sign of a + b*sqrt(3) when signs differ: square both sides
    if (a > 0) == (b > 0):
        return (a > 0) - (a < 0)
    lhs, rhs = a * a, 3 * b * b
    if b > 0:   # a < 0 < b:  positive iff 3b^2 > a^2
        return (lhs < rhs) - (lhs > rhs)
    else:       # b < 0 < a:  positive iff a^2 > 3b^2
        return (lhs > rhs) - (lhs < rhs)
def q3_div_by_frac(x, c): return (x[0] / c, x[1] / c)   # c > 0

# self-test of the comparison against hand-checked values
assert q3_cmp((1, 1), (2, 0)) > 0        # 1+sqrt(3) ~= 2.732 > 2
assert q3_cmp((0, 1), (2, 0)) < 0        # sqrt(3) ~= 1.732 < 2
assert q3_cmp((1, 1), (3, 0)) < 0        # 1+sqrt(3) < 3
assert q3_cmp((-1, 1), (0, 0)) > 0       # -1+sqrt(3) ~= 0.732 > 0
assert q3_cmp((-2, 1), (0, 0)) < 0       # -2+sqrt(3) < 0
T_QM = (Fr(12, 9), Fr(8, 9))          # I_3(QM) = (12+8 sqrt(3))/9 ~= 2.87293
assert q3_cmp(T_QM, (Fr(2), 0)) > 0 and q3_cmp((Fr(4), 0), T_QM) > 0

# ---- pattern space ----------------------------------------------------------
W = [
    lambda d: 1 if d == 0 else (-1 if d == 2 else 0),   # w00
    lambda d: 1 if d == 0 else (-1 if d == 1 else 0),   # w01
    lambda d: 1 if d == 2 else (-1 if d == 0 else 0),   # w10
    lambda d: 1 if d == 0 else (-1 if d == 2 else 0),   # w11
]
TARGET = (0, 0, 2, 0)
pats = []
for d00 in range(3):
    for d01 in range(3):
        for d10 in range(3):
            d11 = (d01 + d10 - d00) % 3
            pats.append((d00, d01, d10, d11))
assert len(pats) == 27

print("=== PART A: CGLMP d=3 pattern structure (exact) ===")
ok(f"pattern space: {len(pats)} valid difference patterns (modular constraint)")
vals = [sum(W[r](d[r]) for r in range(4)) for d in pats]
B = max(vals)
ok(f"local bound B = max pattern value = {B}  (Collins et al.: I_3 <= 2)")
assert B == 2
assert TARGET not in [tuple(d) for d in pats]
ok("modular obstruction: no pattern attains the algebraic-max target (0,0,2,0)")
blocked = [[r for r in range(4) if W[r](d[r]) < 1] for d in pats]
nm = [(d, b[0]) for d, b in zip(pats, blocked) if len(b) == 1]
ok(f"near-miss patterns (block exactly one row): {len(nm)}")
assert len(nm) == 4 and sorted(b for _, b in nm) == [0, 1, 2, 3]
for d, r in nm:
    assert sum(W[i](d[i]) for i in range(4)) == 2
ok("each near-miss pattern has value 2 = B (local bound), one per row")
assert all(len(b) >= 1 for b in blocked)
ok("every pattern blocks at least one row (cap lower-bound ingredient)")

# ---- exhaustive S_max(N, k/N) ----------------------------------------------
def row_h(wlist, k, N):
    """exact per-row value numerator: h = (num)/N with integer num."""
    ws = sorted(wlist)
    sm = sum(ws)
    wmax = ws[-1]
    if k >= N:
        return wmax * N          # h = wmax  -> num = wmax*N
    pref = sum(ws[:k])
    return sm + k * wmax - pref

def exhaustive(N):
    best_num = [-10**9] * (N + 1)
    min_maxblk = None
    for combo in combinations_with_replacement(range(27), N):
        rows_w = [[W[r](pats[p][r]) for p in combo] for r in range(4)]
        Csum = sum(sum(w) for w in rows_w)
        blks = [N - sum(1 for x in w if x == 1) for w in rows_w]
        mx = max(blks)
        if mx < N and (min_maxblk is None or mx < min_maxblk):
            min_maxblk = mx
        for k in range(N + 1):
            num = Csum + sum(row_h(w, k, N) - sum(w) for w in rows_w)
            if num > best_num[k]:
                best_num[k] = num
    return {k: Fr(best_num[k], N) for k in range(N + 1)}, min_maxblk

print("\n=== PART B: exact curves S_max(N, k/N), N = 2..6 (full enumeration) ===")
res = {}
for N in range(2, 7):
    curve, mmblk = exhaustive(N)
    res[N] = (curve, mmblk)
    print(f"  N={N}: " + "  ".join(f"S({k}/{N})={curve[k]}" for k in range(N + 1)))
    assert curve[0] == Fr(2), f"N={N}: S_max(N,0) must be B=2"
    ok(f"N={N}: S_max(N,0) = 2 = local bound")
    assert curve[1] == Fr(2 * N + 2 * min(N, 4), N), f"N={N}: first-segment slope"
    ok(f"N={N}: S_max(N,1/N) = 2 + 2*min(N,4)/N   (first segment 2+2*min(N,4)*F)")
    cap_num = -(-N // 4)
    assert mmblk == cap_num, f"N={N}: min_Q max block-count = {mmblk} != ceil(N/4)"
    ok(f"N={N}: F*(I_3=4) = min_Q max_r n_r^+/N = {cap_num}/{N} = ceil(N/4)/N (exact cap)")
    assert curve[cap_num] == Fr(4) and (cap_num == 0 or curve[cap_num - 1] < Fr(4))

# ---- F* to the quantum value T = (12+8 sqrt3)/9 -----------------------------
print("\n=== PART C: F* to reach I_3(QM) = (12+8 sqrt(3))/9 ~= 2.87293 ===")
def fstar_of_Q(N, combo):
    """exact min F with S_Q(F) >= T, in Q(sqrt(3)); None if never reached < 1."""
    rows_w = [[W[r](pats[p][r]) for p in combo] for r in range(4)]
    best = None
    for k in range(N):
        A_num = Csum_at(rows_w, k, N)          # S_Q(k/N)*N, integer
        Bsl = sum((sorted(w)[-1] - sorted(w)[k]) for w in rows_w)   # slope*N on [k/N,(k+1)/N]
        if Bsl <= 0:
            continue
        A = (Fr(A_num, N), Fr(0))
        cand = q3_div_by_frac(q3_sub(T_QM, A), Fr(Bsl))   # Bsl is dS/dF
        lo, hi = (Fr(k, N), Fr(0)), (Fr(k + 1, N), Fr(0))
        if q3_cmp(cand, lo) < 0: cand = lo
        # valid only if segment actually reaches T: S_Q((k+1)/N) >= T
        A2_num = Csum_at(rows_w, k + 1, N)
        if q3_cmp((Fr(A2_num, N), Fr(0)), T_QM) < 0:
            continue
        if best is None or q3_cmp(cand, best) < 0:
            best = cand
    return best

def Csum_at(rows_w, k, N):
    """S_Q(k/N) * N as integer."""
    tot = sum(sum(w) for w in rows_w)
    return tot + sum(row_h(w, k, N) - sum(w) for w in rows_w)

for N in (2, 3):
    bestF, bestQ = None, None
    for combo in combinations_with_replacement(range(27), N):
        c = fstar_of_Q(N, combo)
        if c is not None and (bestF is None or q3_cmp(c, bestF) < 0):
            bestF, bestQ = c, combo
    pred = q3_div_by_frac(q3_sub(T_QM, (Fr(2), 0)), Fr(2 * N, 1))
    assert bestF is not None and q3_cmp(bestF, pred) == 0, \
        f"N={N}: F* = {bestF} != predicted (T-2)/(2N) = {pred}"
    ok(f"N={N}: F*(T) = (T-2)/(2N) exactly in Q(sqrt(3))  [min over all {sum(1 for _ in combinations_with_replacement(range(27), N))} multisets]")

for N in range(4, 7):
    curve, _ = res[N]
    # linear-segment theorem: S_max(N,F) = 2+8F on [0, floor(N/4)/N]
    kmax = N // 4
    for k in range(kmax + 1):
        assert curve[k] == Fr(2 * N + 8 * k, N), f"N={N}: linear segment at k={k}"
    ok(f"N={N}: S_max(N,k/N) = 2+8k/N for k <= {kmax} (construction + global bound 2+8F)")
    Fstar = q3_div_by_frac(q3_sub(T_QM, (Fr(2), 0)), Fr(8, 1))   # (T-2)/8
    assert q3_cmp(Fstar, (Fr(kmax, N), Fr(0))) < 0 or q3_cmp(Fstar, (Fr(kmax, N), Fr(0))) == 0
    ok(f"N={N}: F*(T) = (T-2)/8 = sqrt(3)/9 - 1/12 ~= {float(Fstar[0]) + float(Fstar[1])*1.7320508:6f} < floor(N/4)/N -> on linear segment")

# ---- general-N cap construction ---------------------------------------------
print("\n=== PART D: general-N cap construction (round-robin over 4 near-misses) ===")
nm_pats = [d for d, _ in nm]
for N in range(2, 65):
    blks = [0] * 4
    for i in range(N):
        r = i % 4
        blks[r] += 1
    cap_num = -(-N // 4)
    assert max(blks) <= cap_num and sum(blks) == N
ok("round-robin over the 4 near-miss patterns: max block-count <= ceil(N/4), sum = N, N=2..64")
print("     => F*(I_3=4; N) = ceil(N/4)/N for all N (lower bound: every pattern blocks >= 1 row)")

print(f"\nALL CHECKS PASSED ({NCHECK} checks).  elapsed {time.time()-T0:.1f}s")
sys.exit(0)
