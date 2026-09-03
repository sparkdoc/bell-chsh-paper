#!/usr/bin/env python3
import os as _os; _ROOT = _os.path.normpath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # release: paths relative to repo root
"""Step-1 math audit (publication review of the Bell/CHSH fine-tuning project).

Exact, from-first-principles re-derivation of S_max(N,F), INDEPENDENT of the LP
pipeline (phase1/2/3.py). All core arithmetic in fractions.Fraction (exact);
float used only for the final cross-checks at irrational F*=(sqrt(2)-1)/m.

Theorem chain (prose proofs in PROOFS.md; this script is the certificate):
  L1  Row TV optimum:  max_{TV(p||rho)<=F} p.f = mu + 2*min(F, m/N)  (k,m > 0)
  T1  S_max(N,F) = max over multisets Q of N even-parity product patterns of
        Smax_Q(Q,F) = S0(Q) + sum_ab g_ab(F),
        g_ab = 2*min(F, n_ab/N) if 0 < n_ab < N else 0,
        n_ab = #{lambda: q_ab != sigma_ab}.
  T2  Cap: F(S=4;N) = ceil(N/4)/N  (parity lemma + round-robin construction).
  T3  Exact curves N=2..8 via per-cell upper envelopes (Fraction-exact).
  C1  Headline F*=(sqrt(2)-1)/min(N,4) for N<=6 (float check w/ margin).
  W   Re-verification of solN{2,3,4} witnesses + exact Q(sqrt(2)) algebra.
  NS  NS sub-case N=2: S_NS(2,F) = 2 for ALL F in [0,1/2] (exact affine cert).

Output: audit1.out in the project folder. Scratch only; no project file touched.
"""
import itertools
import math
from fractions import Fraction as Fr

WIT = _ROOT + "/witnesses/"
LOG = _ROOT + "/logs/"
out = open(LOG + "audit1.out", "w")

def pr(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    out.write(s + "\n")

SIGMA = (1, 1, 1, -1)  # ab order 00, 01, 10, 11
PATTERNS = [q for q in itertools.product((1, -1), repeat=4)
            if q[0] * q[1] * q[2] * q[3] == 1]
assert len(PATTERNS) == 8

def C(q):
    """CHSH value of one deterministic vertex (pattern)."""
    return q[0] + q[1] + q[2] - q[3]

def n_wrong(Q):
    """n_ab = #{lambda : q_ab != sigma_ab} (wrong-sign / 'missing' count)."""
    return tuple(sum(1 for q in Q if q[ab] != SIGMA[ab]) for ab in range(4))

def S0_of(Q):
    N = len(Q)
    return Fr(sum(C(q) for q in Q), N)

def cell_line(Q, k, N):
    """Smax_Q on the cell F in [k/N, (k+1)/N] as (intercept, slope), exact.

    On this cell min(F, n/N) equals n/N (frozen) for 1<=n<=k and equals F
    (active) for k<n<N; rows with n in {0,N} contribute 0 always (Lemma 1).
    """
    n = n_wrong(Q)
    b = S0_of(Q)
    m = 0
    for ab in range(4):
        if 1 <= n[ab] <= k:
            b += Fr(2 * n[ab], N)
        elif k < n[ab] < N:
            m += 2
    return m, b   # (slope, intercept)

def Smax_Q_float(Q, Fv):
    """Direct float evaluation of Smax_Q (independent of cell decomposition)."""
    N = len(Q)
    n = n_wrong(Q)
    tot = float(S0_of(Q))
    for ab in range(4):
        if 0 < n[ab] < N:
            tot += 2.0 * min(Fv, n[ab] / N)
    return tot

def upper_hull(lines):
    """Upper envelope of lines y = m*F + b (m,b Fractions).

    Returns [(line, x_start)] sorted by slope; x_start is Fr or None (-inf).
    Standard monotone chain after dedup by slope (keep max intercept).
    """
    best = {}
    for m, b in lines:
        if m not in best or b > best[m]:
            best[m] = b
    ls = sorted((Fr(m), Fr(b)) for m, b in best.items())  # slope asc
    H = []
    for L in ls:
        if H and H[-1][0] == L[0]:
            continue
        while len(H) >= 2:
            m1, b1 = H[-2]
            m2, b2 = H[-1]
            m3, b3 = L
            x12 = (b1 - b2) / (m2 - m1)   # where H[-1] overtakes H[-2]
            x23 = (b2 - b3) / (m3 - m2)   # where L overtakes H[-1]
            if x12 >= x23:
                H.pop()                   # H[-1] never maximal
            else:
                break
        H.append(L)
    res = []
    for i, L in enumerate(H):
        x = None if i == 0 else (H[i - 1][1] - L[1]) / (L[0] - H[i - 1][0])
        res.append((L, x))
    return res

# ---------------------------------------------------------------------------
pr("=== PART A: exact curves S_max(N,F) for N=2..8 (Fraction arithmetic) ===")
SQ2 = math.sqrt(2)
Qs_by_N = {}
for N in range(2, 9):
    Qs = list(itertools.combinations_with_replacement(PATTERNS, N))
    Qs_by_N[N] = Qs
    pr(f"--- N={N}: {len(Qs)} pattern multisets ---")
    lines_by_cell = [[cell_line(Q, k, N) for Q in Qs] for k in range(N)]

    # Bell bound at F=0
    V0 = max(b for m, b in lines_by_cell[0])
    pr(f"  S_max(N,0) = {V0}   (expect exactly 2)")
    assert V0 == 2, "Bell bound failed"

    # exact piecewise-linear curve on [0,1]
    for k in range(N):
        A, B = Fr(k, N), Fr(k + 1, N)
        hull = upper_hull(lines_by_cell[k])
        pieces = []
        for i, (L, x) in enumerate(hull):
            x_next = hull[i + 1][1] if i + 1 < len(hull) else None
            lo = A if x is None else max(A, x)
            hi = B if x_next is None else min(B, x_next)
            if lo < hi:
                pieces.append((lo, hi, L))
        for lo, hi, (m, b) in pieces:
            pr(f"  F in [{lo}, {hi}]:   S_max = {b} + {m}*F")

    # breakpoint values V_k = S_max(N, k/N), exact
    Vks = []
    for k in range(N + 1):
        kk = min(k, N - 1)
        Vks.append(max(b + m * Fr(k, N) for m, b in lines_by_cell[kk]))
    pr("  S_max at F=k/N: " + ", ".join(str(v) for v in Vks))

    # pin the exact breakpoint values (Theorem 5 / PROOFS.md §8): 2+8F up to the
    # cap ceil(N/4)/N, then flat at 4
    EXPECT_VKS = {
        2: [2, 4, 4],
        3: [2, 4, 4, 4],
        4: [2, 4, 4, 4, 4],
        5: [2, Fr(18, 5)] + [4] * 4,
        6: [2, Fr(10, 3)] + [4] * 5,
        7: [2, Fr(22, 7)] + [4] * 6,
        8: [2, 3] + [4] * 7,
    }
    assert Vks == EXPECT_VKS[N], f"breakpoint values N={N}: got {Vks}"

    # cap value: F(S=4;N) = min_Q max_ab n_ab/N   (T2 necessity+sufficiency)
    cap = Fr(min(max(n_wrong(Q)) for Q in Qs), N)
    cap_formula = Fr(-(-N // 4), N)  # ceil(N/4)/N
    pr(f"  F(S=4;N) exact = {cap}     ceil(N/4)/N = {cap_formula}")
    assert cap == cap_formula, f"cap mismatch N={N}"

# closed-form curve checks for N=2,3,4 (T3): min(2 + 2*min(N,4)*F, 4)
pr("--- T3 closed-form checks ---")
for N in (2, 3, 4):
    mN = min(N, 4)
    ok = True
    for k in range(N + 1):
        expect = min(2 + 2 * mN * Fr(k, N), 4)
        kk = min(k, N - 1)
        vals = [b + m * Fr(k, N) for m, b in
                (cell_line(Q, kk, N) for Q in Qs_by_N[N])]
        if max(vals) != expect:
            ok = False
            pr(f"  MISMATCH N={N} k={k}: got {max(vals)} want {expect}")
    # also check the envelope line on cell 0 is exactly 2 + 2*mN*F
    hull0 = upper_hull([cell_line(Q, 0, N) for Q in Qs_by_N[N]])
    mid0 = Fr(1, 2 * N)  # interior of cell 0
    top0 = max(hull0, key=lambda t: t[0][1] + t[0][0] * mid0)
    pr(f"  N={N}: breakpoints match min(2+{2*mN}*F,4): {ok}; "
       f"cell-0 top line: S = {top0[0][1]} + {top0[0][0]}*F")

# ---------------------------------------------------------------------------
pr("=== PART B: headline F* at the exact values (float cross-check) ===")
for N in (2, 3, 4, 5, 6):
    mN = min(N, 4)
    Fst = (SQ2 - 1) / mN
    vals = sorted(((Smax_Q_float(Q, Fst), i) for i, Q in enumerate(Qs_by_N[N])),
                  reverse=True)
    top, i_top = vals[0]
    gap = top - vals[1][0]
    want = 2 + 2 * mN * Fst  # == 2*sqrt(2) exactly
    pr(f"  N={N}: S_max(F*) = {top:.12f}   (2+{2*mN}*F* = {want:.12f}, "
       f"2*sqrt2 = {2*SQ2:.12f}; margin over 2nd-best Q: {gap:.3e})")
    assert abs(top - 2 * SQ2) < 1e-9, f"Part B: S_max(F*) != 2*sqrt2 at N={N}: {top}"

# ---------------------------------------------------------------------------
pr("=== PART C: re-verification of solN{2,3,4} witnesses (float, raw .npy) ===")
import numpy as np
for N in (2, 3, 4):
    u = np.load(WIT + f"solN{N}_u.npy")   # shape (2,N): u[a,lambda]
    v = np.load(WIT + f"solN{N}_v.npy")   # shape (2,N): v[b,lambda]
    p = np.load(WIT + f"solN{N}_p.npy")   # shape (4,N): p[ab,lambda]
    f = np.stack([u[a] * v[b] for a in range(2) for b in range(2)])  # (4,N)
    E = (p * f).sum(axis=1)   # E[ab] = sum_lambda p[ab,lambda] f[ab,lambda]
    S = float(sum(SIGMA[ab] * E[ab] for ab in range(4)))
    TV = [0.5 * float(np.abs(p[ab] - 1.0 / N).sum()) for ab in range(4)]
    Fv = max(TV)
    Q = sorted(tuple(int(f[ab, l]) for ab in range(4)) for l in range(N))
    SmaxQ = Smax_Q_float(list(Q), Fv)
    mN = min(N, 4)
    Fstar = (SQ2 - 1) / mN
    pr(f"  N={N}: S = {S:.12f}   F = {Fv:.9f}   (F* = {Fstar:.9f})")
    pr(f"        per-row TV: " + ", ".join(f"{t:.6f}" for t in TV))
    pr(f"        product patterns per lambda: {Q}")
    pr(f"        Smax_Q(Q, F) = {SmaxQ:.12f}  -> witness p optimal for its "
       f"vertex: {abs(S - SmaxQ) < 1e-9};  |S - 2sqrt2| = {abs(S-2*SQ2):.2e}")
    # saved arrays come from the bisection search pipeline (F-tolerance 1e-5,
    # see supp Part C). F is asserted at that tolerance directly; S is a
    # downstream quantity with local slope 8 (S = 2+8F on the first segment),
    # so its deviation is amplified by up to ~8x: assert at 1e-4. Observed
    # max over N=2,3,4 is 3.4e-5.
    assert abs(S - 2 * SQ2) < 1e-4, f"Part C: N={N} witness S != 2*sqrt2: {S}"
    assert abs(Fv - Fstar) < 1e-5, f"Part C: N={N} witness F != F*: {Fv}"

# ---------------------------------------------------------------------------
pr("=== PART D: exact Q(sqrt(2)) algebra for the N=4 witness table ===")
def zadd(z, w): return (z[0] + w[0], z[1] + w[1])
def zneg(z):    return (-z[0], -z[1])
# elements are (r, s) = r + s*sqrt(2);  mul not needed for these checks.
q_ = Fr  # shorthand
Fstar_z = (q_(-1, 4), q_(1, 4))          # (sqrt2-1)/4
a_z = (q_(1, 2), q_(-1, 4))              # 1/4 - F* = (2-sqrt2)/4
b_z = (0, q_(1, 4))                      # 1/4 + F* = sqrt2/4
half = (q_(1, 2), 0)
# rows of p (RESULTS.md section 5) and product-sign vectors f per row:
rows_p = [
    ((q_(1, 4), 0), a_z, b_z, (q_(1, 4), 0)),          # 00
    (a_z, (q_(1, 4), 0), b_z, (q_(1, 4), 0)),          # 01
    (b_z, (q_(1, 4), 0), a_z, (q_(1, 4), 0)),          # 10
    (b_z, (q_(1, 4), 0), (q_(1, 4), 0), a_z),          # 11
]
rows_f = [
    (+1, -1, +1, +1),   # q00 of (P2,P4,P3,P1)
    (-1, +1, +1, +1),   # q01
    (+1, +1, -1, +1),   # q10
    (-1, -1, -1, +1),   # q11
]
Ez = []
for pp, ff in zip(rows_p, rows_f):
    e = (0, 0)
    for x, s in zip(pp, ff):
        e = zadd(e, x if s == 1 else zneg(x))
    Ez.append(e)
S_z = zadd(zadd(zadd(Ez[0], Ez[1]), Ez[2]), zneg(Ez[3]))
pr(f"  E[ab] exact: " + ", ".join(f"({r}, {s})*sqrt2-part" for r, s in Ez))
pr(f"  S exact = ({S_z[0]}, {S_z[1]}*sqrt(2))   (expect (0, 2) i.e. 2*sqrt2)")
assert S_z == (0, 2), "N=4 witness table algebra FAILED"
# TV of each row computed EXACTLY from rows_p: TV_ab = (1/2) sum_i |p[ab,i] - 1/4|
def zsign(z):
    x, y = z
    if x == 0 and y == 0: return 0
    if y == 0: return 1 if x > 0 else -1
    if x == 0: return 1 if y > 0 else -1
    lhs, rhs = x * x, 2 * y * y
    if lhs > rhs: return 1 if x > 0 else -1
    if lhs < rhs: return 1 if y > 0 else -1
    raise AssertionError("x^2 == 2y^2 with nonzero rationals — impossible")

def zabs(z):
    s = zsign(z)
    return z if s >= 0 else zneg(z)

q4 = (Fr(1, 4), 0)
tvs_z = []
for pp in rows_p:
    acc = (Fr(0), 0)
    for x in pp:
        acc = zadd(acc, zabs((x[0] - q4[0], x[1] - q4[1])))
    tvs_z.append((acc[0] / 2, acc[1] / 2))
pr(f"  TV per row exact (computed from rows_p): {tvs_z}   (expect all == F* = {Fstar_z})")
assert all(t == Fstar_z for t in tvs_z), f"N=4 witness row TVs != F*: {tvs_z}"

# ---------------------------------------------------------------------------
pr("=== PART E: NS sub-case, N=2 — exact for ALL F in [0,1/2] ===")
# Variables t_ab = p[ab, lambda=0] in [1/2-F, 1/2+F] (TV constraint for N=2).
# S = const + sum_ab c_ab * t_ab (affine). FULL NS (all four equations, as in
# phase3.py lines 43-48): A-side: (u[a,l0]-u[a,l1])(t[a,0]-t[a,1])=0 for a=0,1;
# B-side: (v[b,l0]-v[b,l1])(t[0,b]-t[1,b])=0 for b=0,1.
# Vertex of {box + equalities}: every equality-free variable at a bound.
cands = []  # (A, B) with S(F) = A + B*F, exact Fractions
for u00, u01, u10, u11 in itertools.product((1, -1), repeat=4):
    for v00, v01, v10, v11 in itertools.product((1, -1), repeat=4):
        U = {(0, 0): u00, (0, 1): u01, (1, 0): u10, (1, 1): u11}
        V = {(0, 0): v00, (0, 1): v01, (1, 0): v10, (1, 1): v11}
        f = {}
        for a in range(2):
            for b in range(2):
                for l in range(2):
                    f[(a, b, l)] = U[(a, l)] * V[(b, l)]
        const = Fr(0)
        c = {}
        for a in range(2):
            for b in range(2):
                s = SIGMA[a * 2 + b]
                const += s * f[(a, b, 1)]
                c[(a, b)] = s * (f[(a, b, 0)] - f[(a, b, 1)])
        # var idx: 0=t00, 1=t01, 2=t10, 3=t11;  idx -> (a,b) = (idx//2, idx%2)
        parent = {x: x for x in range(4)}
        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        def union(x, y):
            parent[find(x)] = find(y)
        if u00 != u01: union(0, 1)   # A-side a=0
        if u10 != u11: union(2, 3)   # A-side a=1
        if v00 != v01: union(0, 2)   # B-side b=0
        if v10 != v11: union(1, 3)   # B-side b=1
        groups = {}
        for x in range(4):
            groups.setdefault(find(x), []).append(x)
        glist = list(groups.values())
        cg = {g[0]: sum(c[(x // 2, x % 2)] for x in g) for g in glist}
        for signs in itertools.product((1, -1), repeat=len(glist)):
            A = const
            B = Fr(0)
            for g, sgn in zip(glist, signs):
                cc = cg[g[0]]
                A += Fr(cc, 2)          # t = 1/2 + sgn*F  ->  cc/2 + sgn*cc*F
                B += sgn * cc
            cands.append((A, B))
pr(f"  candidate affine functions S(F)=A+B*F collected: {len(cands)}")
viol = [(A, B) for A, B in set(cands)
        if max(A, A + B * Fr(1, 2)) > 2]   # line <=2 on [0,1/2] iff endpoints
pr(f"  candidates exceeding 2 anywhere on [0,1/2]: {len(viol)}")
for A, B in viol[:5]:
    pr(f"    VIOLATOR: S(F) = {A} + {B}*F")
V_at = {Fv: max(A + B * Fv for A, B in cands) for Fv in (0, Fr(1, 4), Fr(1, 2))}
pr(f"  S_NS(2,F) exact at F=0: {V_at[0]};  F=1/4: {V_at[Fr(1,4)]};  "
   f"F=1/2: {V_at[Fr(1,2)]}")
assert not viol and V_at[0] == 2 and V_at[Fr(1, 4)] == 2 and V_at[Fr(1, 2)] == 2
pr("  CONCLUSION: S_NS(2,F) = 2 for ALL F in [0,1/2], exactly. "
   "(closes the [0.4,0.5] grid gap of ns_gap.out)")

# ---------------------------------------------------------------------------
pr("=== SUMMARY ===")
pr("All exact (Fraction) checks passed: Bell bound N=2..8; cap F(S=4;N)=ceil(N/4)/N")
pr("for N=2..8; T3 closed forms N=2,3,4; N=4 witness table in Q(sqrt2); NS N=2.")
out.close()
print("audit1.out written to", LOG + "audit1.out")
