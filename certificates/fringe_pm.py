#!/usr/bin/env python3
import os as _os; _ROOT = _os.path.normpath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # release: paths relative to repo root
"""fringe_pm.py — exact certificate for the Peres--Mermin square fine-tuning theory.

Companion to review/divergent/fringe.md (serendipity scout, 2026-08-30).

Scenario (the CONTEXTUALITY ANALOG of the CHSH fine-tuning question, task item (a)):
  Peres--Mermin square: 9 binary observables v[i,j] in a 3x3 grid. Six lines:
  rows R_0,R_1,R_2 (target +1), columns C_1,C_2 (target +1), column C_0 (target -1).
  A noncontextual (NC) hidden-variable model preassigns v in {+-1}^9 per state.
  Measurement dependence: the source may depend on which line c is measured,
  p[c, lam] >= 0 with row-sum 1; fine-tuning F = max_c TV(p(.|c) || 1/N),
  uniform baseline (project convention). Game: line c is won iff
  prod_{i in c} v_i = target(c); value G = (1/6) sum_c P(win c).
  NC bound (F=0): G <= 5/6. Quantum (the actual PM observables, ANY state): G = 1
  (products are +I on R0,R1,R2,C1,C2 and -I on C0 — asserted by PM-Qd), so the
  quantum value is the algebraic maximum here (state-independent contextuality).

Results certified below (exact Fraction arithmetic; exit 0 iff every assertion passes):
  PM-1  Parity lemma: every assignment violates an ODD number of lines (never 0);
        all 32 odd-weight miss-sets occur, in particular the six single-line
        "near-miss" types. Consequence S0(Q) = 1 - M/(6N) <= 5/6 for every Q.
  PM-2  Global bound G <= 5/6 + F (proof in fringe.md sec.3; ingredients here).
  PM-3  First segment: G_max(N,F) = 5/6 + F on [0, floor(N/6)/N] for N >= 6;
        for N = 2..5 a single segment G = 5/6 + (N/6)F on [0, 1/N] reaching 1.
  PM-4  Cap / headline: F*_PM(N) := min{F : G_max(N,F) = 1} = ceil(N/6)/N for all
        N >= 2 (parity lower bound + round-robin construction); saturates at 1/6.
  PM-5  Exact curves G_max(N,F), N = 2..8 (required) and 9,10 (best effort):
        full enumeration of achievable miss-count vectors w in N^6; the upper
        envelope is asserted to equal the round-robin curve at EVERY breakpoint
        k/N; the convexity argument of fringe.md sec.5 extends equality to each
        whole cell (same technique as verification/scripts/curves_highN.py).
  PM-6  Explicit witness models at the cap for N = 4..8: response tables + p rows,
        verified exactly (row sums 1, per-row TV = w_c/N, max TV = ceil(N/6)/N, G = 1).
  PM-7  BONUS (CHSH average-vs-worst-case, task item (d), exact Q(sqrt(2))): the
        N=4 witness of RESULTS.md sec.5 has ALL FOUR per-row TVs equal to
        d = (sqrt(2)-1)/4 and S = 2*sqrt(2) exactly, so the minimum AVERAGE-TV
        fine-tuning at N=4 is also (sqrt(2)-1)/4 — the headline is robust to
        max_ab -> (1/4) sum_ab aggregation. (General bound S <= 2 + 8*F_avg proved
        in fringe.md sec.6.)

Soft cross-checks (float, HiGHS via scipy; WARN only, never fail the run):
  LP-1  PM cap: with the round-robin response table at N=6, max G subject to
        per-line TV <= ceil(6/6)/6 = 1/6 should reach 1.
  LP-2  CHSH: at the known-optimal N=4 vertex, min average-TV subject to
        S >= 2*sqrt(2) should equal (sqrt(2)-1)/4.

No floating point in PM-1..PM-7. Re-run: /usr/bin/python3 verification/scripts/fringe_pm.py
"""
import sys
import time
from fractions import Fraction as Fr
from itertools import product

T0 = time.time()
LOG = []

def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.append(s)

CHECKS = [0]
FAILS = [0]

def check(name, cond, detail=""):
    CHECKS[0] += 1
    if cond:
        say(f"  [PASS] {name}")
    else:
        FAILS[0] += 1
        say(f"  [FAIL] {name} {detail}")

# ----------------------------------------------------------------------------
# Part 0/PM-1: the grid, the parity lemma, miss-set enumeration
# ----------------------------------------------------------------------------
say("=" * 78)
say("PART 0 (PM-1): Peres--Mermin square — parity lemma and miss-set census")
say("=" * 78)

# v60 fix (2026-09-01): the -I product is on C0 (XX*YY*ZZ = -I), not C2 — the old
# vector targeted C2, making the certified game have quantum G = 2/3 < 5/6. The cap
# ceil(N/6)/N is parity-invariant under this relabeling (product of targets stays -1).
TARGETS = (1, 1, 1, -1, 1, 1)   # lines: R0,R1,R2,C0,C1,C2 ; C0 target -1
LINENAMES = ("R0", "R1", "R2", "C0", "C1", "C2")

def line_cells(c):
    return [(c, j) for j in range(3)] if c < 3 else [(i, c - 3) for i in range(3)]

# sanity: each cell lies in exactly one row and one column (the parity identity)
cell_lines = {}
for c in range(6):
    for cell in line_cells(c):
        cell_lines.setdefault(cell, []).append(c)
check("each cell in exactly two lines (one row, one col)",
      all(len(v) == 2 and min(v) < 3 and max(v) >= 3 for v in cell_lines.values()))
check("product of the six line targets is -1",
      TARGETS[0] * TARGETS[1] * TARGETS[2] * TARGETS[3] * TARGETS[4] * TARGETS[5] == -1)

miss_masks = set()
near_miss_example = {}          # line c -> one assignment (bit tuple, 1 = value -1) violating ONLY line c
parity_ok = True
for bits in product((0, 1), repeat=9):
    m = 0
    for c in range(6):
        p = 1
        for (i, j) in line_cells(c):
            if bits[3*i + j]:
                p = -p
        if p != TARGETS[c]:
            m |= (1 << c)
    w = bin(m).count("1")
    if w % 2 == 0:
        parity_ok = False
    miss_masks.add(m)
    if w == 1:
        for c in range(6):
            if m == (1 << c) and c not in near_miss_example:
                near_miss_example[c] = bits

check("PM-1a: every one of the 512 assignments violates an ODD number of lines", parity_ok)
odd_masks = {m for m in range(64) if bin(m).count("1") % 2 == 1}
check("PM-1b: exactly the 32 odd-weight miss-sets occur (none missing, none extra)",
      miss_masks == odd_masks, f"got {len(miss_masks)}")
check("PM-1c: a near-miss type exists for every one of the 6 lines",
      set(near_miss_example) == set(range(6)))

say("\n  Near-miss example assignments (bit 1 = value -1; cells in row-major order):")
for c in range(6):
    bits = near_miss_example[c]
    grid = "\n".join("   " + " ".join("-" if bits[3*i+j] else "+" for j in range(3)) for i in range(3))
    say(f"  violates only {LINENAMES[c]}:\n{grid}")

# --- PM-Q: quantum realization — the Peres--Mermin observables, exact check ---
say("")
say("PART 0b (PM-Q): quantum value G = 1 — Peres--Mermin operator identities (exact)")

def _mm(A, B):
    n = len(A)
    C = [[(0, 0)] * n for _ in range(n)]
    for i in range(n):
        for k in range(n):
            if A[i][k] == (0, 0):
                continue
            for j in range(n):
                a, b = A[i][k]; c, d = B[k][j]
                r, im = C[i][j]
                C[i][j] = (r + a * c - b * d, im + a * d + b * c)
    return C

def _kr(A, B):
    n, m = len(A), len(B)
    C = [[(0, 0)] * (m * m) for _ in range(n * m)]
    for i in range(n):
        for k in range(n):
            if A[i][k] == (0, 0):
                continue
            for j in range(m):
                for l in range(m):
                    a, b = A[i][k]; c, d = B[j][l]
                    C[i * m + j][k * m + l] = (a * c - b * d, a * d + b * c)
    return C

_PX = [[(0, 0), (1, 0)], [(1, 0), (0, 0)]]
_PY = [[(0, 0), (0, -1)], [(0, 1), (0, 0)]]
_PZ = [[(1, 0), (0, 0)], [(0, 0), (-1, 0)]]
_PI = [[(1, 0), (0, 0)], [(0, 0), (1, 0)]]
_I4 = _kr(_PI, _PI)
_PN = {"X": _PX, "Y": _PY, "Z": _PZ, "I": _PI}

# one standard Peres--Mermin table (found by backtracking; any relabeling works):
PM_TABLE = [("XX", "XI", "IX"), ("YY", "XZ", "ZX"), ("ZZ", "IZ", "ZI")]
pm_mats = [[_kr(_PN[n[0]], _PN[n[1]]) for n in row] for row in PM_TABLE]

def _is_id(M, sign):
    return all(M[i][j] == ((sign if i == j else 0), 0) for i in range(4) for j in range(4))

_pm_lines = [pm_mats[r] for r in range(3)] + [[pm_mats[r][c] for r in range(3)] for c in range(3)]
_pm_ok_commute = all(_mm(A, B) == _mm(B, A) for line in _pm_lines for A, B in
                     [(line[0], line[1]), (line[0], line[2]), (line[1], line[2])])
_pm_prods = [_mm(_mm(L[0], L[1]), L[2]) for L in _pm_lines]
_pm_negs = [c for c, P in enumerate(_pm_prods) if _is_id(P, -1)]
check("PM-Qa: all 6 line triples are jointly measurable (pairwise commuting)", _pm_ok_commute)
check("PM-Qb: every line product is exactly +I or -I (operator identity)",
      all(_is_id(P, 1) or _is_id(P, -1) for P in _pm_prods))
check("PM-Qc: exactly one line has product -I (the parity obstruction)", len(_pm_negs) == 1,
      f"neg lines {_pm_negs}")
check("PM-Qd: TARGETS[c] == sign(line product c) for all c (game matches the operators)",
      all(TARGETS[c] == (-1 if _is_id(_pm_prods[c], -1) else 1) for c in range(6)),
      f"targets {TARGETS} vs products "
      + ",".join("-I" if _is_id(P, -1) else "+I" for P in _pm_prods))
say(f"  Table (rows R0-R2, cols C0-C2): {PM_TABLE}")
say(f"  Line products: " + ", ".join(f"{LINENAMES[c]}={'-I' if _is_id(P,-1) else '+I'}"
                                     for c, P in enumerate(_pm_prods)))
say("  => in ANY quantum state, measuring the three observables of line c gives outcomes")
say("     whose product equals the target t_c with probability 1: quantum G = 1 (algebraic max).")

MISSVECS = [tuple((m >> c) & 1 for c in range(6)) for m in sorted(miss_masks)]

# ----------------------------------------------------------------------------
# Parts PM-2..PM-5: exact curves by miss-count-vector enumeration
# ----------------------------------------------------------------------------
say("")
say("=" * 78)
say("PARTS PM-2..PM-5: exact G_max(N,F) — miss-count-vector DP + breakpoint envelope")
say("=" * 78)

def addw(w, v):
    return tuple(sorted(a + b for a, b in zip(w, v)))   # canonical: sorted coords

def G_of(w, N, k):
    """Exact value at F = k/N for miss-count vector w (sorted)."""
    M = sum(w)
    return Fr(1) - Fr(M, 6 * N) + Fr(1, 6) * sum(Fr(min(k, wc), N) for wc in w)

def round_robin_w(N):
    k, r = divmod(N, 6)
    return tuple(sorted((k + (1 if c < r else 0)) for c in range(6)))

REQUIRED_N = list(range(2, 9))     # N = 2..8: assertions required
OPTIONAL_N = [9, 10]               # best effort
DP_CAP = 400_000                   # safety cap on the DP set size

for N in REQUIRED_N + OPTIONAL_N:
    optional = N in OPTIONAL_N
    S = {tuple([0] * 6)}
    for step in range(N):
        S = {addw(w, v) for w in S for v in MISSVECS}
        if len(S) > DP_CAP:
            say(f"  [SKIP] N={N}: DP set exceeded {DP_CAP} at step {step+1}; stopping enumeration.")
            break
    else:
        ws = round_robin_w(N)
        check(f"N={N}: round-robin miss-count vector is achievable", ws in S, f"{ws}")
        W = [max(G_of(w, N, kk) for w in S) for kk in range(N + 1)]
        Cf = [G_of(ws, N, kk) for kk in range(N + 1)]
        check(f"N={N}: envelope == round-robin curve at all {N+1} breakpoints k/N", W == Cf,
              f"W={[str(x) for x in W]}")
        m = -(-N // 6)   # ceil(N/6)
        check(f"N={N}: NC bound exact, G_max(N,0) = 5/6", W[0] == Fr(5, 6))
        check(f"N={N}: cap reached, G_max(N, ceil(N/6)/N) = 1 (ceil(N/6)={m})", W[m] == 1)
        check(f"N={N}: strictly below 1 just before the cap, G_max(N,(ceil-1)/N) < 1",
              W[m - 1] < 1, f"got {W[m-1]}")
        # first-segment statement PM-3: for N >= 6, G = 5/6 + F on [0, floor(N/6)/N];
        # for N = 2..5 the round-robin has only N active lines and reaches 1 at F = 1/N.
        if N >= 6:
            kfloor = N // 6
            seg_ok = all(W[kk] == Fr(5, 6) + Fr(kk, N) for kk in range(kfloor + 1))
            check(f"N={N}: first segment G = 5/6 + F on [0, {kfloor}/{N}] (PM-3)", seg_ok)
        else:
            seg_ok = W[1] == Fr(5, 6) + Fr(N, 6) * Fr(1, N) and W[1] == 1
            check(f"N={N}: single segment G = 5/6 + (N/6)F on [0, 1/N], reaching 1 at the cap (PM-3)", seg_ok)
        # print the piecewise-linear curve (true slopes: N * per-breakpoint increment)
        pieces = []
        prev_slope = None
        for kk in range(N):                      # cell [kk/N, (kk+1)/N]
            slope = N * (W[kk + 1] - W[kk])
            if slope != prev_slope:
                pieces.append(f"[{kk}/{N},{kk+1}/{N}]: {W[kk]} + {slope}*F")
                prev_slope = slope
        say(f"  N={N}: w*={ws}, ceil(N/6)/N={Fr(m, N)}")
        say("    " + " ; ".join(pieces))

# ----------------------------------------------------------------------------
# Part PM-6: explicit witness models at the cap (N = 4..8)
# ----------------------------------------------------------------------------
say("")
say("=" * 78)
say("PART PM-6: explicit cap witnesses — response tables + p rows, exact verification")
say("=" * 78)

def build_witness(N):
    """Round-robin: state i carries the near-miss assignment of line (i mod 6).
    At F* = ceil(N/6)/N each row c is supported uniformly on its satisfying states."""
    states = [near_miss_example[i % 6] for i in range(N)]
    masks = []
    for bits in states:
        m = 0
        for c in range(6):
            p = 1
            for (i, j) in line_cells(c):
                if bits[3*i + j]:
                    p = -p
            if p != TARGETS[c]:
                m |= (1 << c)
        masks.append(m)
    rows = []
    for c in range(6):
        sup = [i for i in range(N) if not (masks[i] & (1 << c))]
        wc = N - len(sup)
        pvec = tuple(Fr(1, len(sup)) if i in sup else Fr(0) for i in range(N))
        rows.append((wc, pvec))
    return states, masks, rows

for N in (4, 5, 6, 7, 8):
    states, masks, rows = build_witness(N)
    mstar = -(-N // 6)
    ok_rowsum = all(sum(p) == 1 for _, p in rows)
    tvs = []
    ok_tv = True
    for c, (wc, pvec) in enumerate(rows):
        tv = Fr(1, 2) * sum(abs(x - Fr(1, N)) for x in pvec)
        tvs.append(tv)
        if tv != Fr(wc, N):
            ok_tv = False
    ok_maxtv = max(tvs) == Fr(mstar, N)
    # success probability of each line under its row (must be 1: support on satisfiers)
    ok_win = True
    for c in range(6):
        wc, pvec = rows[c]
        winmass = Fr(0)
        for i in range(N):
            if not (masks[i] & (1 << c)):
                winmass += pvec[i]
        if winmass != 1:
            ok_win = False
    check(f"N={N}: row sums 1", ok_rowsum)
    check(f"N={N}: per-row TV = w_c/N exactly, i.e. {tvs}", ok_tv)
    check(f"N={N}: max per-row TV = ceil(N/6)/N = {Fr(mstar, N)}", ok_maxtv)
    check(f"N={N}: every line won with probability 1 (G = 1)", ok_win)

say("\n  Full witness at N=4 (states = near-miss assignments of lines R0,R1,R2,C0;")
say("  p rows below give the exact cap model, F* = 1/4):")
states, masks, rows = build_witness(4)
for i, bits in enumerate(states):
    grid = " ".join("-" if bits[3*a+b] else "+" for a in range(3) for b in range(3))
    say(f"  state {i} (violates only {LINENAMES[i]}): v = [{grid}]")
for c in range(6):
    wc, pvec = rows[c]
    say(f"  row {LINENAMES[c]}: p = [{', '.join(str(x) for x in pvec)}]   (w={wc}, TV={Fr(wc,4)})")

# ----------------------------------------------------------------------------
# Part PM-7 (bonus): CHSH average-vs-worst-case at N=4, exact Q(sqrt(2))
# ----------------------------------------------------------------------------
say("")
say("=" * 78)
say("PART PM-7 (BONUS, task item d): CHSH N=4 witness — per-row TVs all equal, Q(sqrt(2))")
say("=" * 78)

sys.path.insert(0, _ROOT + "/certificates")
from qsqrt2 import zadd, zsub, zmul, zr, zsign   # noqa: E402

d = (Fr(-1, 4), Fr(1, 4))          # (sqrt2 - 1)/4
g = (Fr(1, 4), Fr(0))              # 1/4
a = zsub(g, d)                     # 1/4 - d
b = zadd(g, d)                     # 1/4 + d

prows = [
    (g, a, b, g),    # row 00
    (a, g, b, g),    # row 01
    (b, g, a, g),    # row 10
    (b, g, g, a),    # row 11
]
u = [(1, -1, -1, -1), (1, 1, 1, -1)]     # u[a, lam]
v = [(1, 1, -1, -1), (-1, -1, -1, -1)]   # v[b, lam]

def zabs(x):
    s = zsign(x)
    return x if s >= 0 else zsub((Fr(0), Fr(0)), x)

E = []
for ab in range(4):
    acc = (Fr(0), Fr(0))
    for lam in range(4):
        sign = u[ab // 2][lam] * v[ab % 2][lam]
        acc = zadd(acc, zmul(zr(sign), prows[ab][lam]))
    E.append(acc)
Sval = zadd(zadd(zadd(E[0], E[1]), E[2]), zsub((Fr(0), Fr(0)), E[3]))
check("PM-7a: S = 2*sqrt(2) exactly (Q(sqrt(2)) tuple (0, 2))", Sval == (Fr(0), Fr(2)), f"got {Sval}")

tvs_z = []
for ab in range(4):
    l1 = (Fr(0), Fr(0))
    for lam in range(4):
        l1 = zadd(l1, zabs(zsub(prows[ab][lam], g)))
    tvs_z.append(zsub(l1, zmul(l1, (Fr(1, 2), Fr(0)))))   # TV = L1/2
check("PM-7b: all four per-row TVs equal d = (sqrt2-1)/4 exactly",
      all(t == d for t in tvs_z), f"got {tvs_z}")
Favg = zadd(zadd(tvs_z[0], tvs_z[1]), zadd(tvs_z[2], tvs_z[3]))
Favg = zmul(Favg, (Fr(1, 4), Fr(0)))
check("PM-7c: average-TV cost of the witness = d (so F*_avg(N=4) = d = F*(N=4))",
      Favg == d, f"got {Favg}")

# ----------------------------------------------------------------------------
# Soft LP cross-checks (float; WARN only)
# ----------------------------------------------------------------------------
say("")
say("=" * 78)
say("SOFT LP CROSS-CHECKS (HiGHS via scipy; float; never fail the run)")
say("=" * 78)
try:
    from scipy.optimize import linprog

    # LP-1: PM cap feasibility at N=6 with the round-robin table.
    states, masks, rows = build_witness(6)
    win = [[1 if not (masks[i] & (1 << c)) else 0 for i in range(6)] for c in range(6)]
    # vars: x[c*6+i] = p[c, i] (36 of them); y[c*6+i] >= |p - 1/6| (36 of them)
    c_obj = [-win[c][i] for c in range(6) for i in range(6)] + [0.0] * 36   # max G = (1/6) sum winmass
    A_eq, b_eq = [], []
    A_ub, b_ub = [], []
    for c in range(6):
        A_eq.append([1.0 if (j // 6) == c else 0.0 for j in range(72)])
        b_eq.append(1.0)
    for c in range(6):
        for i in range(6):
            r1 = [0.0] * 72; r2 = [0.0] * 72
            r1[c * 6 + i] = 1.0;  r1[36 + c * 6 + i] = -1.0   # y >= p - 1/6
            r2[c * 6 + i] = -1.0; r2[36 + c * 6 + i] = -1.0   # y >= 1/6 - p
            A_ub.append(r1); b_ub.append(1.0 / 6.0)
            A_ub.append(r2); b_ub.append(-1.0 / 6.0)   # -p - y <= -1/6  i.e.  y >= 1/6 - p
        r3 = [0.0] * 72
        for i in range(6):
            r3[36 + c * 6 + i] = 1.0
        A_ub.append(r3); b_ub.append(2.0 * Fr(1, 6))           # sum y <= 2*TV <= 2F*
    res = linprog(c_obj, A_eq=A_eq, b_eq=b_eq, A_ub=A_ub, b_ub=b_ub,
                  bounds=[(0.0, None)] * 72, method="highs")
    if res.status == 0:
        Gval = -res.fun / 6.0
        check("LP-1 (soft): PM N=6 cap — max G at F*=1/6 reaches 1", Gval > 1.0 - 1e-7, f"G={Gval}")
    else:
        say(f"  [WARN] LP-1 did not converge: {res.message}")

    # LP-2: CHSH N=4 — min average TV subject to S >= 2*sqrt(2), fixed optimal vertex.
    import math
    SQ = 2.0 ** 0.5
    dfl = (SQ - 1) / 4.0
    pvals = [[float(x[0]) + float(x[1]) * SQ for x in row] for row in prows]
    # vars: x[ab*4+lam] = p ; y[ab*4+lam] >= |p - 1/4|
    c_obj2 = [0.0] * 32
    for idx in range(16):
        c_obj2[16 + idx] = 1.0 / 8.0     # (1/4) * (1/2) * sum |p - 1/4| over all ab,lam
    A_eq2, b_eq2 = [], []
    for ab in range(4):
        r = [0.0] * 32
        for lam in range(4):
            r[ab * 4 + lam] = 1.0
        A_eq2.append(r); b_eq2.append(1.0)
    A_ub2, b_ub2 = [], []
    for ab in range(4):
        for lam in range(4):
            r1 = [0.0] * 32; r2 = [0.0] * 32
            r1[ab * 4 + lam] = 1.0;  r1[16 + ab * 4 + lam] = -1.0
            r2[ab * 4 + lam] = -1.0; r2[16 + ab * 4 + lam] = -1.0
            A_ub2.append(r1); b_ub2.append(0.25)
            A_ub2.append(r2); b_ub2.append(-0.25)      # -p - y <= -1/4  i.e.  y >= 1/4 - p
    # S >= 2*sqrt(2):  S = sum_ab sigma_ab * sum_lam p[ab,lam] u[a,lam] v[b,lam]
    sig = [1, 1, 1, -1]
    rS = [0.0] * 32
    for ab in range(4):
        for lam in range(4):
            rS[ab * 4 + lam] = sig[ab] * u[ab // 2][lam] * v[ab % 2][lam]
    A_ub2.append([-x for x in rS]); b_ub2.append(-2.0 * SQ + 1e-9)
    res2 = linprog(c_obj2, A_eq=A_eq2, b_eq=b_eq2, A_ub=A_ub2, b_ub=b_ub2,
                   bounds=[(0.0, None)] * 32, method="highs")
    if res2.status == 0:
        check("LP-2 (soft): CHSH N=4 min average-TV at S>=2*sqrt(2) equals d",
              abs(res2.fun - dfl) < 1e-7, f"got {res2.fun}, want {dfl}")
    else:
        say(f"  [WARN] LP-2 did not converge: {res2.message}")
except Exception as e:   # noqa: BLE001 — soft section must never fail the run
    say(f"  [WARN] LP cross-checks skipped ({type(e).__name__}: {e})")

# ----------------------------------------------------------------------------
say("")
say("=" * 78)
dt = time.time() - T0
say(f"TOTAL: {CHECKS[0]} checks, {FAILS[0]} failures, wall time {dt:.1f} s")
if FAILS[0] == 0:
    say("ALL PASS")
sys.exit(1 if FAILS[0] else 0)
