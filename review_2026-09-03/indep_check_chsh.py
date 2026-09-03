#!/usr/bin/env python3
"""
Independent re-verification of the main paper (paper/main.tex), written from
scratch on 2026-09-03 without reusing any project code.  Exact arithmetic uses
Fraction and a tiny Q(sqrt2) class; LP cross-checks use scipy/HiGHS.

Run:  python3 indep_check_chsh.py
Exit code 0 iff every assertion passes.  Findings that are *not* assertion
failures are printed with the prefix  NOTE:
"""
import itertools, math, sys
from fractions import Fraction as Fr
import numpy as np
from scipy.optimize import linprog

FAILS = []
def check(name, cond, detail=""):
    tag = "PASS" if cond else "FAIL"
    print(f"[{tag}] {name} {detail}")
    if not cond:
        FAILS.append(name)

# ---------------------------------------------------------------- Q(sqrt2)
class Q2:
    """a + b*sqrt(2) with a,b Fractions"""
    __slots__ = ("a", "b")
    def __init__(self, a, b=0):
        self.a, self.b = Fr(a), Fr(b)
    def __add__(s, o): o = _q(o); return Q2(s.a+o.a, s.b+o.b)
    __radd__ = __add__
    def __sub__(s, o): o = _q(o); return Q2(s.a-o.a, s.b-o.b)
    def __rsub__(s, o): return _q(o) - s
    def __neg__(s): return Q2(-s.a, -s.b)
    def __mul__(s, o): o = _q(o); return Q2(s.a*o.a + 2*s.b*o.b, s.a*o.b + s.b*o.a)
    __rmul__ = __mul__
    def __truediv__(s, o):
        o = _q(o); n = o.a*o.a - 2*o.b*o.b
        return s * Q2(o.a/n, -o.b/n)
    def __eq__(s, o): o = _q(o); return s.a == o.a and s.b == o.b
    def __float__(s): return float(s.a) + float(s.b)*math.sqrt(2)
    def sign(s):
        # exact sign of a + b sqrt2
        if s.b == 0: return (s.a > 0) - (s.a < 0)
        if s.a == 0: return (s.b > 0) - (s.b < 0)
        if (s.a > 0) == (s.b > 0): return 1 if s.a > 0 else -1
        # opposite signs: compare a^2 vs 2 b^2
        if s.a*s.a > 2*s.b*s.b: return 1 if s.a > 0 else -1
        return 1 if s.b > 0 else -1
    def __abs__(s): return s if s.sign() >= 0 else -s
    def __lt__(s, o): return (s - _q(o)).sign() < 0
    def __le__(s, o): return (s - _q(o)).sign() <= 0
    def __repr__(s): return f"({s.a} + {s.b}*sqrt2)"
def _q(x): return x if isinstance(x, Q2) else Q2(x)
SQRT2 = Q2(0, 1)

# ---------------------------------------------------------------- patterns
SIGMA = (1, 1, 1, -1)          # CHSH sign vector, order 00,01,10,11
def pattern(u0, u1, v0, v1):
    return (u0*v0, u0*v1, u1*v0, u1*v1)
ALL = sorted({pattern(*t) for t in itertools.product((1, -1), repeat=4)})
check("exactly 8 realizable patterns", len(ALL) == 8, ALL)
check("all realizable patterns have even parity", all(math.prod(q) == 1 for q in ALL))
def C(q): return q[0]+q[1]+q[2]-q[3]
check("C(q) in {+2,-2}", all(C(q) in (2, -2) for q in ALL))
MISS = [q for q in ALL if C(q) == 2]
check("four miss-type patterns each differ from sigma at exactly one place",
      len(MISS) == 4 and all(sum(a != b for a, b in zip(q, SIGMA)) == 1 for q in MISS))
MA, MB, MC, MD = (-1, 1, 1, -1), (1, -1, 1, -1), (1, 1, -1, -1), (1, 1, 1, 1)
check("miss types match paper Eq.(miss)", set(MISS) == {MA, MB, MC, MD})

# ---------------------------------------------------------------- Lemma 1 by LP
def lp_row_max(f, N, F, extra_eq=None):
    """max sum p_i f_i  s.t. p in simplex, TV(p||uniform) <= F. Variables p (N), t (N) with
       t_i >= |p_i - 1/N|, sum t <= 2F."""
    c = np.concatenate([-np.array(f, float), np.zeros(N)])
    A_ub, b_ub = [], []
    for i in range(N):
        r = np.zeros(2*N); r[i] = 1; r[N+i] = -1; A_ub.append(r); b_ub.append(1/N)
        r = np.zeros(2*N); r[i] = -1; r[N+i] = -1; A_ub.append(r); b_ub.append(-1/N)
    r = np.zeros(2*N); r[N:] = 1; A_ub.append(r); b_ub.append(2*F)
    A_eq = [np.concatenate([np.ones(N), np.zeros(N)])]; b_eq = [1.0]
    res = linprog(c, A_ub=np.array(A_ub), b_ub=b_ub, A_eq=np.array(A_eq), b_eq=b_eq,
                  bounds=[(0, None)]*(2*N), method="highs")
    assert res.status == 0, res.message
    return -res.fun

def lemma1(f, N, F):
    k = sum(1 for x in f if x == 1); m = N - k
    mu = Fr(k-m, N)
    if k == 0 or m == 0: return mu
    return mu + 2*min(F, Fr(m, N))

rng = np.random.default_rng(1)
maxdev = 0.0; cnt = 0
for N in range(2, 8):
    for f in itertools.product((1, -1), repeat=N):
        for F in [Fr(j, 4*N) for j in range(0, 4*N+1)] + [Fr(1, 3), Fr(2, 7)]:
            if F > 1: continue
            v = lp_row_max(f, N, float(F)); e = float(lemma1(f, N, F))
            maxdev = max(maxdev, abs(v-e)); cnt += 1
check("Lemma 1 (per-row TV optimum) vs LP", maxdev < 1e-9, f"{cnt} cases, max dev {maxdev:.2e}")

# ---------------------------------------------------------------- Theorem 1 / curves
def Smax_formula(N, F):
    """max over multisets Q of S0(Q) + sum_ab g_ab(F); F Fraction; exact."""
    best = None
    for Q in itertools.combinations_with_replacement(range(8), N):
        S0 = Fr(sum(C(ALL[i]) for i in Q), N)
        tot = S0
        for ab in range(4):
            n = sum(1 for i in Q if ALL[i][ab] != SIGMA[ab])
            if 0 < n < N: tot += 2*min(F, Fr(n, N))
        if best is None or tot > best: best = tot
    return best

def paper_curve(N, F):
    F = Fr(F)
    if N == 2: return min(2+4*F, Fr(4))
    if N == 3: return min(2+6*F, Fr(4))
    if N == 4: return min(2+8*F, Fr(4))
    if N == 5: return min(2+8*F, Fr(16, 5)+2*F, Fr(4))
    if N == 6: return min(2+8*F, Fr(8, 3)+4*F, Fr(4))
    if N == 7: return min(2+8*F, Fr(16, 7)+6*F, Fr(4))
    if N == 8: return min(2+8*F, Fr(4))

# brute force over ALL response vertices with per-row LP for N=2,3 (independent of Lemma 1 and of the
# multiset reduction): S_max(N,F) = max over (u,v) of sum_ab max_p sigma_ab * sum p q_ab
def Smax_bruteforce_lp(N, F):
    best = -9
    for bits in itertools.product((1, -1), repeat=4*N):
        u0, u1, v0, v1 = bits[:N], bits[N:2*N], bits[2*N:3*N], bits[3*N:]
        S = 0.0
        for ab, (ua, vb) in enumerate(((u0, v0), (u0, v1), (u1, v0), (u1, v1))):
            f = [SIGMA[ab]*ua[i]*vb[i] for i in range(N)]
            S += lp_row_max(f, N, F)
        best = max(best, S)
    return best

for N in (2, 3):
    dev = 0
    for F in [0, 0.1, 0.2, 1/3, 0.25, 0.45, 0.5, 0.7]:
        dev = max(dev, abs(Smax_bruteforce_lp(N, F) - float(paper_curve(N, Fr(F).limit_denominator(10**6)))))
    check(f"N={N}: brute-force vertex+LP S_max equals paper curve", dev < 1e-8, f"max dev {dev:.1e}")

for N in range(2, 9):
    grid = sorted({Fr(k, 2*N) for k in range(0, 2*N+1)} | {Fr(1, 4), Fr(1, 3), Fr(2, 5), Fr(2, 7), Fr(1, 2)})
    bad = [(F, Smax_formula(N, F), paper_curve(N, F)) for F in grid if Smax_formula(N, F) != paper_curve(N, F)]
    check(f"N={N}: exhaustive-multiset S_max equals paper's piecewise-linear curve on grid+midpoints",
          not bad, str(bad[:3]))

# breakpoints quoted in Theorem (curves)
check("S_max(5,1/5)=18/5", Smax_formula(5, Fr(1, 5)) == Fr(18, 5))
check("S_max(6,1/6)=10/3", Smax_formula(6, Fr(1, 6)) == Fr(10, 3))
check("S_max(7,1/7)=22/7", Smax_formula(7, Fr(1, 7)) == Fr(22, 7))
check("S_max(5,1/4)=3.7 (cap remark)", Smax_formula(5, Fr(1, 4)) == Fr(37, 10))

# cap: min over Q of max_ab n_ab == ceil(N/4)
for N in range(2, 9):
    best = min(max(sum(1 for i in Q if ALL[i][ab] != SIGMA[ab]) for ab in range(4))
               for Q in itertools.combinations_with_replacement(range(8), N))
    check(f"N={N}: cap min_Q max_ab n_ab = ceil(N/4)", best == -(-N//4), f"got {best}")

# headline F*: smallest F with S_max >= 2 sqrt2 ; check S_max(N,F*) == 2sqrt2 exactly and S_max(N,F*-eps) < 2sqrt2
def Fstar(N): return (SQRT2 - 1) / min(N, 4)
for N in range(2, 9):
    Fs = Fstar(N)
    # evaluate the formula in Q(sqrt2): S0(Q) + sum 2 min(F, n/N) -- need exact comparison of F* with n/N
    def Smax_q2(F):
        best = None
        for Q in itertools.combinations_with_replacement(range(8), N):
            tot = Q2(Fr(sum(C(ALL[i]) for i in Q), N))
            for ab in range(4):
                n = sum(1 for i in Q if ALL[i][ab] != SIGMA[ab])
                if 0 < n < N:
                    tot = tot + 2*(F if F <= Q2(Fr(n, N)) else Q2(Fr(n, N)))
            if best is None or best < tot: best = tot
        return best
    eps = Q2(Fr(1, 10**6))
    check(f"N={N}: S_max(N,F*) == 2sqrt2 exactly", Smax_q2(Fs) == 2*SQRT2, repr(Smax_q2(Fs)))
    check(f"N={N}: S_max(N,F*-1e-6) < 2sqrt2", Smax_q2(Fs - eps) < 2*SQRT2)

# ---------------------------------------------------------------- explicit N=4 signalling witness (Sec. VI)
Fs4 = Fstar(4); a = Q2(Fr(1, 4)) - Fs4; b = Q2(Fr(1, 4)) + Fs4; q4 = Q2(Fr(1, 4))
P = [[q4, a, b, q4], [a, q4, b, q4], [b, q4, a, q4], [b, q4, q4, a]]
U = [(1, -1, -1, -1), (1, 1, 1, -1)]; V = [(1, 1, -1, -1), (-1, -1, -1, -1)]
def E_of(P, U, V, ab):
    a_, b_ = divmod(ab, 2)
    tot = Q2(0)
    for lam in range(len(P[ab])):
        tot = tot + P[ab][lam] * (U[a_][lam]*V[b_][lam])
    return tot
def S_of(P, U, V): return sum((SIGMA[ab]*E_of(P, U, V, ab) for ab in range(4)), Q2(0))
def TV_row(row, N):
    tot = Q2(0)
    for x in row: tot = tot + abs(x - Q2(Fr(1, N)))
    return tot / 2
check("Sec.VI witness: rows sum to 1", all(sum(r, Q2(0)) == 1 for r in P))
check("Sec.VI witness: rows nonnegative", all(x.sign() >= 0 for r in P for x in r))
check("Sec.VI witness: S = 2sqrt2 exactly", S_of(P, U, V) == 2*SQRT2, repr(S_of(P, U, V)))
check("Sec.VI witness: E[ab] = (s,s,s,-s), s=sqrt2/2",
      [E_of(P, U, V, ab) == SIGMA[ab]*SQRT2/2 for ab in range(4)] == [True]*4)
check("Sec.VI witness: every row TV = F*", all(TV_row(r, 4) == Fs4 for r in P))
def NS_residuals(P, U, V, N):
    res = []
    for a_ in range(2):   # Alice marginal: rows a0 vs a1
        m = [sum((P[2*a_+b_][l]*U[a_][l] for l in range(N)), Q2(0)) for b_ in range(2)]
        res.append(m[0]-m[1])
    for b_ in range(2):   # Bob marginal: rows 0b vs 1b
        m = [sum((P[2*a_+b_][l]*V[b_][l] for l in range(N)), Q2(0)) for a_ in range(2)]
        res.append(m[0]-m[1])
    return res
r = NS_residuals(P, U, V, 4)
print("NOTE: Sec.VI (signalling-allowed) witness NS residuals:", r, " (nonzero expected)")
# pairwise TVs (Sec. VII.5 claim: values in {delta, 2delta})
pw = sorted({(TV_row_pair := (sum((abs(x-y) for x, y in zip(P[i], P[j])), Q2(0))/2)).__repr__()
             for i in range(4) for j in range(i+1, 4)})
pwvals = [sum((abs(x-y) for x, y in zip(P[i], P[j])), Q2(0))/2 for i in range(4) for j in range(i+1, 4)]
check("Sec.VI witness: pairwise TVs in {F*, 2F*}", all(v == Fs4 or v == 2*Fs4 for v in pwvals))
check("Sec.VI witness: M_model = 2F*", max(pwvals, key=float) == 2*Fs4)

# ---------------------------------------------------------------- NS witnesses (Thm NS headline)
def ns_witness_general(N):
    d = Fstar(4)
    cls = [i % 4 for i in range(N)]
    n = [cls.count(j) for j in range(4)]
    U = [tuple((1, 1, -1, -1)[c] for c in cls), tuple((-1, 1, 1, -1)[c] for c in cls)]
    V = [tuple((-1, 1, -1, -1)[c] for c in cls), tuple((1, -1, -1, -1)[c] for c in cls)]
    X = {(0, 0): -1, (0, 1): +1, (1, 0): +1, (1, 1): -1, (2, 1): +1, (2, 2): -1, (3, 0): +1, (3, 3): -1}
    P = []
    for ab in range(4):
        row = []
        for i in range(N):
            x = X.get((ab, cls[i]))
            row.append(Q2(Fr(1, N)) + (d*x/n[cls[i]] if x else Q2(0)))
        P.append(row)
    return P, U, V, d
for N in range(4, 14):
    P, U, V, d = ns_witness_general(N)
    ok = (all(sum(r, Q2(0)) == 1 for r in P) and all(x.sign() >= 0 for r in P for x in r)
          and S_of(P, U, V) == 2*SQRT2 and all(TV_row(r, N) == d for r in P)
          and all(x == 0 for x in NS_residuals(P, U, V, N)))
    pats = {tuple(U[a_][i]*V[b_][i] for a_, b_ in ((0, 0), (0, 1), (1, 0), (1, 1))) for i in range(N)}
    check(f"NS witness N={N}: rows ok, S=2sqrt2, TV=F*, NS residuals 0", ok)
    if N == 4: check("NS witness N=4 patterns are exactly M_A..M_D", pats == {MA, MB, MC, MD})
# N=3 NS witness (Eq. ns3gauge / ns3rows)
d3 = (SQRT2 - 1)/3; t = Q2(Fr(1, 3))
coef = [[-1, Fr(1, 2), Fr(1, 2)], [Fr(1, 2), -1, Fr(1, 2)], [-Fr(1, 2), Fr(1, 2), 0], [Fr(1, 2), Fr(1, 2), -1]]
P3 = [[t + d3*c for c in row] for row in coef]
U3 = [(-1, -1, 1), (1, -1, 1)]; V3 = [(1, -1, 1), (-1, 1, 1)]
check("NS N=3 witness: rows sum 1 & nonneg", all(sum(r, Q2(0)) == 1 for r in P3) and all(x.sign() >= 0 for r in P3 for x in r))
check("NS N=3 witness: S = 2sqrt2", S_of(P3, U3, V3) == 2*SQRT2, repr(S_of(P3, U3, V3)))
check("NS N=3 witness: per-row TV = (d3,d3,d3/2,d3)", [TV_row(r, 3) for r in P3] == [d3, d3, d3/2, d3])
check("NS N=3 witness: NS residuals 0", all(x == 0 for x in NS_residuals(P3, U3, V3, 3)))
pats3 = {tuple(U3[a_][i]*V3[b_][i] for a_, b_ in ((0, 0), (0, 1), (1, 0), (1, 1))) for i in range(3)}
check("NS N=3 witness patterns = {M_A,M_B,M_D}", pats3 == {MA, MB, MD}, str(pats3))
pw3 = [sum((abs(x-y) for x, y in zip(P3[i], P3[j])), Q2(0))/2 for i in range(4) for j in range(i+1, 4)]
check("NS N=3 witness: M_model = 3/2 F*(3)", max(pw3, key=float) == d3*Fr(3, 2))

# ---------------------------------------------------------------- NS N=2 flatness by LP over all vertices
def ns_lp_max(N, F, U, V, onesided=False):
    """max S over rows p[ab] (4N vars) + TV aux (4N) with NS equalities."""
    nv = 8*N
    c = np.zeros(nv)
    for ab in range(4):
        a_, b_ = divmod(ab, 2)
        for l in range(N): c[ab*N+l] = -SIGMA[ab]*U[a_][l]*V[b_][l]
    A_ub, b_ub, A_eq, b_eq = [], [], [], []
    for ab in range(4):
        for l in range(N):
            r = np.zeros(nv); r[ab*N+l] = 1; r[4*N+ab*N+l] = -1; A_ub.append(r); b_ub.append(1/N)
            r = np.zeros(nv); r[ab*N+l] = -1; r[4*N+ab*N+l] = -1; A_ub.append(r); b_ub.append(-1/N)
        r = np.zeros(nv); r[4*N+ab*N:4*N+(ab+1)*N] = 1; A_ub.append(r); b_ub.append(2*F)
        r = np.zeros(nv); r[ab*N:(ab+1)*N] = 1; A_eq.append(r); b_eq.append(1.0)
    for a_ in range(2):   # Alice: rows a0, a1
        r = np.zeros(nv)
        for l in range(N): r[(2*a_)*N+l] += U[a_][l]; r[(2*a_+1)*N+l] -= U[a_][l]
        A_eq.append(r); b_eq.append(0.0)
    if not onesided:
        for b_ in range(2):   # Bob: rows 0b, 1b
            r = np.zeros(nv)
            for l in range(N): r[b_*N+l] += V[b_][l]; r[(2+b_)*N+l] -= V[b_][l]
            A_eq.append(r); b_eq.append(0.0)
    res = linprog(c, A_ub=np.array(A_ub), b_ub=b_ub, A_eq=np.array(A_eq), b_eq=b_eq,
                  bounds=[(0, None)]*nv, method="highs")
    if res.status != 0: return -np.inf
    return -res.fun
def ns_max_over_vertices(N, F, onesided=False):
    best = -9
    for bits in itertools.product((1, -1), repeat=4*N):
        U = [bits[:N], bits[N:2*N]]; V = [bits[2*N:3*N], bits[3*N:]]
        best = max(best, ns_lp_max(N, F, U, V, onesided))
    return best
for F in (0.0, 0.1, 0.25, 0.4, 0.5):
    v = ns_max_over_vertices(2, F)
    check(f"NS N=2, F={F}: max S = 2 (flat)", abs(v-2) < 1e-9, f"got {v:.12f}")
    v1 = ns_max_over_vertices(2, F, onesided=True)
    check(f"one-sided-NS N=2, F={F}: max S = 2+4F (side finding)", abs(v1-(2+4*F)) < 1e-9, f"got {v1:.12f}")
# NS N=3: check S_NS(3,F) = 2+6F at a few F <= F*(3)  (paper: "common curve is 2+6F" by scaling)
for F in (0.05, 0.1, float((SQRT2-1)/3)):
    v = ns_max_over_vertices(3, F)
    check(f"NS N=3, F={F:.5f}: max S = 2+6F", abs(v-(2+6*F)) < 1e-8, f"got {v:.12f}")

# ---------------------------------------------------------------- NS S=4 witnesses (Thm ns4cap) printed in paper
def check_S4(name, P, U, V, N, cap):
    Pq = [[Q2(x) for x in r] for r in P]
    ok = (all(sum(r, Q2(0)) == 1 for r in Pq) and S_of(Pq, U, V) == 4
          and all(x == 0 for x in NS_residuals(Pq, U, V, N)))
    tvs = [TV_row(r, N) for r in Pq]
    check(f"{name}: S=4, NS residuals 0, rows sum 1", ok)
    check(f"{name}: max row TV = {cap}", max(tvs, key=float) == Q2(cap), str(tvs))
# N=3: classes A,B,C (sizes 1,1,1,0), flip C.  Base gauge per type from realisation rule of Thm1 (iii)
def base_uv(q):  # u0=+1, v0=q00, v1=q01, u1=q00*q10
    return (1, q[0]*q[2]), (q[0], q[1])
def build_uv(pats, flipclass):
    U0, U1, V0, V1 = [], [], [], []
    for i, q in enumerate(pats):
        (u0, u1), (v0, v1) = base_uv(q)
        s = -1 if q == flipclass else 1
        U0.append(s*u0); U1.append(s*u1); V0.append(s*v0); V1.append(s*v1)
    return [tuple(U0), tuple(U1)], [tuple(V0), tuple(V1)]
U, V = build_uv([MA, MB, MC], MC)
P = [[0, Fr(1, 2), Fr(1, 2)], [Fr(1, 2), 0, Fr(1, 2)], [Fr(1, 2), Fr(1, 2), 0], [Fr(1, 2), Fr(1, 2), 0]]
check_S4("NS S=4 witness N=3", P, U, V, 3, Fr(1, 3))
U, V = build_uv([MA, MB, MC, MD], MD)
P = [[0, Fr(1, 4), Fr(1, 4), Fr(1, 2)], [Fr(1, 4), 0, Fr(1, 4), Fr(1, 2)],
     [Fr(1, 4), Fr(1, 2), 0, Fr(1, 4)], [Fr(1, 4), Fr(1, 2), Fr(1, 4), 0]]
check_S4("NS S=4 witness N=4 (PR-box decomposition)", P, U, V, 4, Fr(1, 4))
U, V = build_uv([MA, MB, MC, MD, MA], MD)
cm = [[0, Fr(1, 5), Fr(3, 10), Fr(1, 2)], [Fr(1, 5), 0, Fr(3, 10), Fr(1, 2)],
      [Fr(2, 5), Fr(1, 2), 0, Fr(1, 10)], [Fr(2, 5), Fr(1, 2), Fr(1, 10), 0]]
sizes = [2, 1, 1, 1]; cls = [0, 1, 2, 3, 0]
P = [[cm[ab][cls[i]]/sizes[cls[i]] for i in range(5)] for ab in range(4)]
check_S4("NS S=4 witness N=5", P, U, V, 5, Fr(2, 5))
# does the N=4 S=4 witness reproduce the PR box?  P(x,y|ab)
Pq = [[Q2(x) for x in r] for r in P] if False else None

# ---------------------------------------------------------------- one-sided (Sec. VII.2) by brute force LP
def os_max(N, F):
    best = -9
    for bits in itertools.product((1, -1), repeat=4*N):
        u0, u1, v0, v1 = bits[:N], bits[N:2*N], bits[2*N:3*N], bits[3*N:]
        fA = [u0[i]*v0[i] + u0[i]*v1[i] for i in range(N)]
        fB = [u1[i]*v0[i] - u1[i]*v1[i] for i in range(N)]
        best = max(best, lp_row_max(fA, N, F) + lp_row_max(fB, N, F))
    return best
for N in (2, 3):
    for F in (0.0, 0.1, 0.2, 1/3, 0.4, 0.5):
        v = os_max(N, F)
        exp = min(2+4*F, 4) if N == 2 else (2+4*F if F <= 1/3 else min(10/3 + 2*(F-1/3), 4))
        check(f"one-sided N={N}, F={F:.4f}: S_os,max = {exp:.6f}", abs(v-exp) < 1e-8, f"got {v:.10f}")
# one-sided remark: cap ratio  ceil(N/2)/N  over  ceil(N/4)/N  for N = 1,2 mod 4
for N in range(2, 15):
    ratio = Fr(-(-N//2)) / Fr(-(-N//4))
    if N % 4 in (1, 2):
        print(f"NOTE: one-sided/full cap ratio N={N}: {ratio} (paper claims >= 3/2 for N = 1,2 mod 4)")

# ---------------------------------------------------------------- mutual information of witnesses (float)
def MI_bits(P):
    P = np.array([[float(x) for x in r] for r in P]); avg = P.mean(axis=0)
    H = lambda p: -sum(x*math.log2(x) for x in p if x > 0)
    return H(avg) - np.mean([H(r) for r in P])
P4, U4, V4, _ = ns_witness_general(4)
print(f"NOTE: MI bits: Sec.VI witness {MI_bits([[q4, a, b, q4], [a, q4, b, q4], [b, q4, a, q4], [b, q4, q4, a]]):.6f}, "
      f"NS N=4 witness {MI_bits(P4):.6f}, NS N=3 witness {MI_bits(P3):.6f}  (paper: 0.056036 / 0.056036 / 0.054810)")
# factorisation claim: p(x,y|lambda) = p(x|lambda)p(y|lambda) fails for deterministic-response models?
# Sec. VII.5 says H&B's causal structure imposes p(x,y|lambda)=p(x|lambda)p(y|lambda).  In H&B's notation x,y are
# SETTINGS; in main.tex x,y are OUTCOMES (Sec. II), for which the factorisation is trivially true (point masses).
# Check the meaningful version: does p(a,b|lambda) factorise for the witnesses?  (uniform settings prior)
def settings_factorise(P):
    P = np.array([[float(x) for x in r] for r in P]); N = P.shape[1]; worst = 0
    for l in range(N):
        col = P[:, l]/P[:, l].sum()          # p(ab|lambda) for ab = 00,01,10,11
        worst = max(worst, abs(col[0]*col[3]-col[1]*col[2]))
    return worst
print(f"NOTE: settings-factorisation defect max|p00 p11 - p01 p10| : Sec.VI witness {settings_factorise([[q4, a, b, q4], [a, q4, b, q4], [b, q4, a, q4], [b, q4, q4, a]]):.4f}, "
      f"NS N=4 {settings_factorise(P4):.4f}, NS N=3 {settings_factorise(P3):.4f} (nonzero => not H&B-causal; but note the x,y notation clash)")

print("\nSUMMARY:", "ALL PASS" if not FAILS else f"{len(FAILS)} FAIL(S): {FAILS}")
sys.exit(0 if not FAILS else 1)
