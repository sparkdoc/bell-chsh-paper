#!/usr/bin/env python3
import os as _os; _ROOT = _os.path.normpath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # release: paths relative to repo root
"""
metric_compare.py — Metric conversion: our TV-from-rho metric F  vs  Hall's M  vs  Hall-Branciard MI.

Purpose (divergent explorer #3): the paper currently calls our per-setting TV metric
"incommensurable" with Hall's M and with Hall & Branciard's mutual information. On
overlapping model classes (local deterministic CHSH models) there ARE explicit
inequalities relating them, and our explicit witnesses can be EVALUATED in the other
metrics. This script certifies those conversions and evaluations.

Conventions (all stated precisely in review/divergent/metric.md):
  * Our metric:      F       = max_ab TV(p(.|ab) || rho),  rho uniform,  TV = (1/2)||p-q||_1.
  * Hall's M (2010): M_full  = max_{ab,a'b'} ||p(.|ab) - p(.|a'b')||_1   (FULL L1, in [0,2]).
  * Task's analog:   M_model = max_{ab,a'b'} TV(p(.|ab) || p(.|a'b'))    (pairwise, NO baseline).
    Note M_full = 2 * M_model  (L1 = 2*TV for the same pair).
  * Hall-Branciard:  I(X,Y;Lambda) in bits, X,Y uniform independent settings.

Exactness:
  * Parts A-F use EXACT Q(sqrt(2)) arithmetic: a value x = a + b*sqrt(2) is stored as a
    pair of Fractions (a, b). No floating point anywhere in the exact parts. This is the
    project's standard certificate convention (see audit1.py); "Fraction arithmetic" here
    means exact rational / quadratic-field arithmetic.
  * Part G (mutual information) uses mpmath at 60 decimal digits (high-precision NUMERICS,
    stated as such -- MI involves log2 of Q(sqrt(2)) numbers with no closed form).

Exit code 0 iff every assertion passes.
"""
import sys
from fractions import Fraction as F

# ---------------------------------------------------------------------------
# PART A: exact Q(sqrt(2)) arithmetic.  x = a + b*sqrt(2),  a,b in Q.
# ---------------------------------------------------------------------------
def qadd(x, y):   return (x[0] + y[0], x[1] + y[1])
def qsub(x, y):   return (x[0] - y[0], x[1] - y[1])
def qneg(x):      return (-x[0], -x[1])
def qmul(x, y):   # (a+b r)(c+d r) = (ac+2bd) + (ad+bc) r ,  r=sqrt(2)
    a, b = x; c, d = y
    return (a * c + 2 * b * d, a * d + b * c)
def qdiv_int(x, n):  return (x[0] / F(n), x[1] / F(n))

def qsign(x):
    """Exact sign of a+b*sqrt(2): +1 / -1 / 0."""
    a, b = x
    if a == 0 and b == 0:
        return 0
    if b == 0:
        return 1 if a > 0 else -1
    if a == 0:
        return 1 if b > 0 else -1
    if (a > 0) == (b > 0):          # same sign -> sign of a
        return 1 if a > 0 else -1
    # opposite signs: |x| = ||a|-|b| r| ; compare a^2 with 2 b^2 exactly
    return (1 if a > 0 else -1) if a * a >= 2 * b * b else (1 if b > 0 else -1)

def qabs(x):
    s = qsign(x)
    return x if s >= 0 else qneg(x)

def qcmp(x, y):
    """sign of x-y."""
    return qsign(qsub(x, y))

# small exact constants (a + b*sqrt(2)) as pairs of Fractions
def Q(a_num, a_den=1, b_num=0, b_den=1):
    return (F(a_num, a_den), F(b_num, b_den))

ZERO = Q(0); ONE = Q(1)
SQRT2   = Q(0, 1, 1, 1)              # sqrt(2)
TWO_SQRT2 = Q(0, 1, 2, 1)            # 2*sqrt(2)

# ---------------------------------------------------------------------------
# PART B: the four witnesses (exact p-tables in Q(sqrt(2)), plus response tables).
#   rows of p are ab in order 00,01,10,11 ; columns are lambda = 1..N.
# ---------------------------------------------------------------------------
def build_witnesses():
    g   = Q(1, 4)                     # 1/4
    a_s = Q(2, 4, -1, 4)              # (2-sqrt2)/4
    b_s = Q(0, 1, 1, 4)               # sqrt2/4
    d   = Q(-1, 4, 1, 4)              # F*_N4 = (sqrt2-1)/4

    # N=4 signalling witness (RESULTS.md Sec.5 / paper main.tex Sec.6)
    sig_p = [
        [g, a_s, b_s, g],
        [a_s, g,   b_s, g],
        [b_s, g,   a_s, g],
        [b_s, g,   g,   a_s],
    ]
    sig_u = [[1, -1, -1, -1],
             [1,  1,  1, -1]]         # rows a=0,1 ; cols lambda=1..4
    sig_v = [[1, 1, -1, -1],
             [-1,-1, -1, -1]]         # rows b=0,1

    # N=4 NS witness (PROOFS.md Sec.9.5 construction at n=(1,1,1,1))
    ns_p = [
        [qsub(g, d), qadd(g, d), g, g],
        [qadd(g, d), qsub(g, d), g, g],
        [g, qadd(g, d), qsub(g, d), g],
        [qadd(g, d), g, g, qsub(g, d)],
    ]
    ns_u = [[1, 1, -1, -1],
            [-1, 1, 1, -1]]
    ns_v = [[-1, 1, -1, -1],
            [1, -1, -1, -1]]

    # N=2 round-robin (patterns {M_C,M_D}), F*=(sqrt2-1)/2
    h   = Q(1, 2)                     # 1/2
    F2  = Q(-1, 2, 1, 2)              # (sqrt2-1)/2
    rr2_p = [
        [h, h],
        [h, h],
        [qsub(h, F2), qadd(h, F2)],
        [qadd(h, F2), qsub(h, F2)],
    ]
    rr2_u = [[1, 1],
             [-1, 1]]
    rr2_v = [[1, 1],
             [1, 1]]

    # N=3 round-robin (patterns {M_A,M_C,M_D}), F*=(sqrt2-1)/3
    t   = Q(1, 3)                     # 1/3
    F3  = Q(-1, 3, 1, 3)              # (sqrt2-1)/3
    rr3_p = [
        [qsub(t, F3), qadd(t, F3), t],
        [t, t, t],
        [t, qsub(t, F3), qadd(t, F3)],
        [qadd(t, F3), t, qsub(t, F3)],
    ]
    rr3_u = [[1, 1, 1],
             [-1, -1, 1]]
    rr3_v = [[-1, 1, 1],
             [1, 1, 1]]

    # N=3 NS witness (ns_general.py Part B closed form; Theorem 7 construction):
    # per-row TV (d3, d3, d3/2, d3), F = d3 = F*. Deterministic +/-1 responses
    # loaded as integers from the saved arrays (exact).
    t   = Q(1, 3)                     # 1/3
    F3b = Q(-1, 3, 1, 3)              # (sqrt2-1)/3
    hd  = qdiv_int(F3b, 2)            # d3/2
    ns3_p = [
        [qsub(t, F3b), qadd(t, hd),   qadd(t, hd)],
        [qadd(t, hd),  qsub(t, F3b),  qadd(t, hd)],
        [qsub(t, hd),  qadd(t, hd),   t],
        [qadd(t, hd),  qadd(t, hd),   qsub(t, F3b)],
    ]
    import numpy as np
    WIT = _ROOT + "/witnesses/"
    ns3_u = [[int(x) for x in row] for row in np.load(WIT + "ns_solN3_u.npy")]
    ns3_v = [[int(x) for x in row] for row in np.load(WIT + "ns_solN3_v.npy")]

    return {
        "N4-signalling": dict(p=sig_p, u=sig_u, v=sig_v, Fstar=d,   N=4),
        "N4-NS":         dict(p=ns_p,  u=ns_u,  v=ns_v, Fstar=d,   N=4),
        "N2-roundrobin": dict(p=rr2_p, u=rr2_u, v=rr2_v, Fstar=F2, N=2),
        "N3-roundrobin": dict(p=rr3_p, u=rr3_u, v=rr3_v, Fstar=F3, N=3),
        "N3-NS":         dict(p=ns3_p, u=ns3_u, v=ns3_v, Fstar=F3b, N=3),
    }

# ---------------------------------------------------------------------------
# exact model metrics (Q(sqrt(2)))
# ---------------------------------------------------------------------------
def row_sum(row):
    s = ZERO
    for x in row:
        s = qadd(s, x)
    return s

def tv_row_uniform(row, N):
    """TV(p(.|ab) || uniform) = (1/2) sum |p - 1/N|, exact."""
    u = Q(1, N)
    s = ZERO
    for x in row:
        s = qadd(s, qabs(qsub(x, u)))
    return qdiv_int(s, 2)

def tv_pair(rowA, rowB):
    """TV(p||q) = (1/2) sum |p_i - q_i|, exact."""
    s = ZERO
    for x, y in zip(rowA, rowB):
        s = qadd(s, qabs(qsub(x, y)))
    return qdiv_int(s, 2)

def S_value(p, u, v):
    """S = E00+E01+E10-E11, exact. p rows ab in order 00,01,10,11."""
    sig = [1, 1, 1, -1]
    S = ZERO
    for ab in range(4):
        a, b = divmod(ab, 2)
        E = ZERO
        row = p[ab]
        for lam in range(len(row)):
            prod = u[a][lam] * v[b][lam]          # +-1
            term = qmul(row[lam], Q(prod))
            E = qadd(E, term)
        S = qadd(S, qmul(E, Q(sig[ab])))
    return S

def pairwise_tvs(p):
    """All 6 pairwise TVs (rows ab order 00,01,10,11)."""
    out = {}
    for i in range(4):
        for j in range(i + 1, 4):
            out[(i, j)] = tv_pair(p[i], p[j])
    return out

def ns_residuals(p, u, v, N):
    """Four NS residuals (exact). A-side: E[A|a,b] indep of b; B-side: E[B|a,b] indep of a.
       residual_a = |sum_lam p[0b]u[a] - sum_lam p[1b]u[a]| for each a  (b fixed -> compare)
       We return the four scalar residuals:
         rA_a = | (sum_lam p[ab=0?]) ... |.  Concretely:
         A-side for setting a:  M(a,b) = sum_lam p[ab, lam] u[a,lam];  residual = |M(a,0)-M(a,1)|.
         B-side for setting b:  N(a,b) = sum_lam p[ab, lam] v[b,lam];  residual = |N(0,b)-N(1,b)|."""
    res = []
    for a in range(2):
        def M(b):
            s = ZERO
            ab = a * 2 + b
            for lam in range(N):
                s = qadd(s, qmul(p[ab][lam], Q(u[a][lam])))
            return s
        res.append(qabs(qsub(M(0), M(1))))
    for b in range(2):
        def Nn(a):
            s = ZERO
            ab = a * 2 + b
            for lam in range(N):
                s = qadd(s, qmul(p[ab][lam], Q(v[b][lam])))
            return s
        res.append(qabs(qsub(Nn(0), Nn(1))))
    return res

def fmt(x):
    a, b = x
    if b == 0:
        return str(a)
    sign = "+" if b > 0 else "-"
    bb = abs(b)
    return f"({a} {sign} {bb}*sqrt2)"

# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    print("=" * 78)
    print("METRIC COMPARISON  (exact Q(sqrt2) parts + mpmath MI)")
    print("=" * 78)

    # ---- self-test of the Q(sqrt2) arithmetic -----------------------------
    assert qsign(Q(1, 2, -1)) == -1           # 1/2 - sqrt2 = -0.914 < 0
    assert qsign(qsub(Q(1), SQRT2)) == -1     # 1 - sqrt2 < 0
    assert qsign(qsub(SQRT2, Q(1))) == 1      # sqrt2 - 1 > 0
    assert qabs(qsub(Q(1), SQRT2)) == qsub(SQRT2, Q(1))
    assert qmul(SQRT2, SQRT2) == Q(2)         # (sqrt2)^2 = 2
    assert qadd(qsub(Q(3), SQRT2), SQRT2) == Q(3)
    print("[A] Q(sqrt2) arithmetic self-test: PASS")

    W = build_witnesses()
    names = ["N4-signalling", "N4-NS", "N2-roundrobin", "N3-roundrobin", "N3-NS"]
    # witnesses that saturate the triangle bound M_model == 2F* exactly:
    SATURATING = ["N4-signalling", "N4-NS", "N2-roundrobin", "N3-roundrobin"]

    results = {}
    for name in names:
        w = W[name]
        p, u, v, Fstar, N = w["p"], w["u"], w["v"], w["Fstar"], w["N"]

        # row sums == 1
        for ab in range(4):
            assert qcmp(row_sum(p[ab]), ONE) == 0, f"{name} row {ab} sum != 1"

        # per-row TV from uniform -> F
        tvs = [tv_row_uniform(p[ab], N) for ab in range(4)]
        Fmax = tvs[0]
        for t in tvs:
            if qcmp(t, Fmax) > 0:
                Fmax = t
        # S
        S = S_value(p, u, v)

        # pairwise TVs -> M_model (TV convention) and M_full (=2*M_model, Hall L1)
        pw = pairwise_tvs(p)
        Mmodel = list(pw.values())[0]
        for t in pw.values():
            if qcmp(t, Mmodel) > 0:
                Mmodel = t
        Mfull = qmul(Mmodel, Q(2))

        # NS residuals
        nsres = ns_residuals(p, u, v, N)
        is_ns = all(qcmp(r, ZERO) == 0 for r in nsres)

        results[name] = dict(F=Fmax, S=S, Mmodel=Mmodel, Mfull=Mfull,
                             tvs=tvs, pw=pw, nsres=nsres, is_ns=is_ns, Fstar=Fstar)

    # ---- print + assert the exact core facts ------------------------------
    for name in names:
        r = results[name]
        w = W[name]
        print(f"\n--- {name}  (N={w['N']}, {'NS' if r['is_ns'] else 'signalling'}) ---")
        print(f"  F* (target)      = {fmt(r['Fstar'])}")
        print(f"  per-row TV       = {[fmt(t) for t in r['tvs']]}")
        print(f"  F = max row TV   = {fmt(r['F'])}")
        print(f"  S                = {fmt(r['S'])}")
        pw = r["pw"]
        lab = {(0,1):"00-01",(0,2):"00-10",(0,3):"00-11",(1,2):"01-10",(1,3):"01-11",(2,3):"10-11"}
        print("  pairwise TV      = " + ", ".join(f"{lab[k]}:{fmt(pw[k])}" for k in lab))
        print(f"  M_model (max TV) = {fmt(r['Mmodel'])}")
        print(f"  M_full (=2*M_mod)= {fmt(r['Mfull'])}   [Hall L1 convention, in [0,2]]")
        print(f"  NS residuals     = {[fmt(x) for x in r['nsres']]}")

        # (i) F equals the headline F* exactly
        assert qcmp(r["F"], r["Fstar"]) == 0, f"{name}: F != F*"
        # (ii) S = 2 sqrt(2) exactly
        assert qcmp(r["S"], TWO_SQRT2) == 0, f"{name}: S != 2*sqrt2"
        # (iii) MAIN LEMMA: M_model <= 2F   (triangle inequality through rho)
        twoF = qmul(r["F"], Q(2))
        assert qcmp(r["Mmodel"], twoF) <= 0, f"{name}: M_model > 2F"
        # (iv) Hall's bound S <= min{2+3*M_full, 4}  (no NS required; see report)
        hall_bound = qadd(Q(2), qmul(r["Mfull"], Q(3)))
        assert qcmp(r["S"], hall_bound) <= 0, f"{name}: S > 2+3*M_full"
        assert qcmp(r["S"], Q(4)) <= 0

    # ---- the saturation fact: M_model == 2F* for the enumerated witnesses --
    # SCOPED claim: saturated exactly by the four closed-form witnesses below
    # (and, within search tolerance, by the saved bisection-produced solN3
    # array). NOT a universal statement over all F*-optimal witnesses — N3-NS
    # is an in-family counterexample with M_model/F* = 3/2 (checked next).
    print("\n[saturation] checking M_model == 2*F* exactly for the enumerated witnesses:")
    for name in SATURATING:
        r = results[name]
        twoFstar = qmul(r["Fstar"], Q(2))
        ok = qcmp(r["Mmodel"], twoFstar) == 0
        print(f"  {name}: M_model={fmt(r['Mmodel'])}  2*F*={fmt(twoFstar)}  -> {'EQUAL' if ok else 'NOT equal'}")
        assert ok, f"{name}: M_model != 2*F*"
    # N3-NS: does NOT saturate — exact value (sqrt2-1)/2 = (3/2) F*(N=3)
    r = results["N3-NS"]
    half = Q(-1, 2, 1, 2)              # (sqrt2-1)/2
    assert r["is_ns"], "N3-NS witness must have all four NS residuals exactly 0"
    assert qcmp(r["Mmodel"], half) == 0, f"N3-NS: M_model != (sqrt2-1)/2 exactly: {fmt(r['Mmodel'])}"
    print(f"  N3-NS: M_model={fmt(r['Mmodel'])} = (sqrt2-1)/2 exactly; "
          f"M_model/F* = 3/2 — NOT saturated (slack F*/2 in the triangle bound)")
    # saved bisection-produced solN3 array (float cross-check, not exact):
    import math as _math
    import numpy as _np
    WIT = _ROOT + "/witnesses/"
    p3f = _np.load(WIT + "solN3_p.npy")
    Mmodel_f = max(0.5 * float(_np.abs(p3f[i] - p3f[j]).sum()) for i in range(4) for j in range(i + 1, 4))
    Fstar_f = (_math.sqrt(2) - 1) / 3
    assert abs(Mmodel_f - 2 * Fstar_f) < 1e-4, \
        f"solN3 saved array: M_model {Mmodel_f} not within 1e-4 of 2*F*"
    print(f"  solN3 (saved bisection array): M_model={Mmodel_f:.9f} vs 2*F*={2 * Fstar_f:.9f} "
          f"(diff {abs(Mmodel_f - 2 * Fstar_f):.1e}; saturates within the 1e-5 search tolerance)")
    print("  => triangle bound M_model <= 2F is SATURATED by every enumerated closed-form witness "
          "(and the saved solN3 array, within tolerance); NOT by N3-NS (ratio 3/2). PASS")

    # ---- converse bounds --------------------------------------------------
    print("\n[converse] F vs M_model:")
    for name in names:
        w = W[name]; r = results[name]; p, N = w["p"], w["N"]
        # p_avg = (1/4) sum_ab p_ab
        cols = len(p[0])
        pavg = []
        for lam in range(cols):
            s = ZERO
            for ab in range(4):
                s = qadd(s, p[ab][lam])
            pavg.append(qdiv_int(s, 4))
        # TV(p_avg || uniform)  (the "center offset" c)
        c = tv_row_uniform(pavg, N)
        # F <= M_model + c   (triangle through p_avg)
        assert qcmp(r["F"], qadd(r["Mmodel"], c)) <= 0, f"{name}: F > M_model + center"
        # adaptive baseline: F_avg = max_ab TV(p_ab || p_avg) <= M_model
        favg_max = tv_pair(p[0], pavg)
        for ab in range(1, 4):
            t = tv_pair(p[ab], pavg)
            if qcmp(t, favg_max) > 0:
                favg_max = t
        assert qcmp(favg_max, r["Mmodel"]) <= 0, f"{name}: F_avg > M_model"
        print(f"  {name}: center TV(p_avg||rho)={fmt(c)}; "
              f"F={fmt(r['F'])} <= M_model+center={fmt(qadd(r['Mmodel'],c))}; "
              f"F_avg(max TV to p_avg)={fmt(favg_max)} <= M_model={fmt(r['Mmodel'])}")

    # counterexample: fixed baseline, all four rows = non-uniform q -> M_model=0, F>0
    q = [Q(1, 2), Q(1, 4), Q(1, 8), Q(1, 8)]   # sums to 1, not uniform (N=4)
    assert qcmp(row_sum(q), ONE) == 0
    allq = [list(q) for _ in range(4)]
    pw_c = pairwise_tvs(allq)
    Mmodel_c = max(pw_c.values(), key=lambda x: (x[0], x[1]))
    assert qcmp(Mmodel_c, ZERO) == 0, "counterexample M_model should be 0"
    F_c = tv_row_uniform(q, 4)
    assert qcmp(F_c, ZERO) > 0, "counterexample F should be > 0"
    print(f"  counterexample (all rows = {fmt(q[0])},...): M_model=0 but F={fmt(F_c)} > 0")
    print("  => NO bound F <= f(M_model) exists for a FIXED external baseline. PASS")

    # ---- Hall's bound interpretation --------------------------------------
    print("\n[Hall bound] S <= min{2+3*M_full, 4}:")
    for name in ["N4-signalling", "N4-NS"]:
        r = results[name]
        hb = qadd(Q(2), qmul(r["Mfull"], Q(3)))   # 2 + 3*(sqrt2-1) = 3*sqrt2 - 1
        print(f"  {name}: M_full={fmt(r['Mfull'])} -> 2+3*M_full={fmt(hb)}; "
              f"actual S=2*sqrt2={fmt(TWO_SQRT2)} (slack={fmt(qsub(hb, r['S']))})")

    # ---- PART G: mutual information (mpmath, high precision) --------------
    import mpmath as mp
    mp.mp.dps = 60
    def d(x, n=6):
        """fixed-decimal display of an mpf (values are O(1); float is fine for printing)."""
        return f"{float(x):.{n}f}"
    SQ = mp.sqrt(mp.mpf(2))
    def log2(x):
        return mp.log(x, 2)
    def to_mp(x):
        a, b = x
        return mp.mpf(a.numerator) / mp.mpf(a.denominator) + \
               (mp.mpf(b.numerator) / mp.mpf(b.denominator)) * SQ
    def Hbits(dist):
        h = mp.mpf(0)
        for p in dist:
            if p > 0:
                h -= p * log2(p)
        return h

    print("\n" + "=" * 78)
    print("PART G: I(X,Y;Lambda) in bits   (mpmath, dps=60; high-precision numerics)")
    print("=" * 78)

    def mutual_info(p):
        cols = len(p[0])
        plam = [sum(p[ab][lam] for ab in range(4)) / 4 for lam in range(cols)]
        Hl = Hbits(plam)
        Hcond = sum(Hbits(p[ab]) for ab in range(4)) / 4
        return Hl - Hcond, plam

    mi = {}
    for name in names:
        w = W[name]
        pmp = [[to_mp(x) for x in row] for row in w["p"]]
        I, plam = mutual_info(pmp)
        mi[name] = I
        print(f"  {name}:  I(X,Y;Lambda) = {d(I)} bits   "
              f"(H(Lambda)={d(Hbits(plam))}, avg H(Lambda|ab)={d(Hbits(plam) - I)})")

    # Hall-Branciard reference minima at S = 2 sqrt(2)
    S = mp.mpf(2) * SQ
    def hbin(p):
        return -p * log2(p) - (1 - p) * log2(1 - p) if 0 < p < 1 else mp.mpf(0)
    # retrocausal: I_R(S) = 2 - h((4-S)/8) - ((4+S)/8) log2 3      [H&B Eq (29)]
    IR = 2 - hbin((4 - S) / 8) - ((4 + S) / 8) * log2(3)
    # causal: I_C(S) = I1(S) = 2 - 2 h(sqrt((4-S)/8)) for S < S0   [H&B Eq (41),(45)]
    pX = mp.sqrt((4 - S) / 8)
    IC = 2 - 2 * hbin(pX)
    print(f"\n  Hall-Branciard reference at S=2*sqrt2={d(S)}:")
    print(f"    I_R (retrocausal min, Eq 29) = {d(IR)} bits")
    print(f"    I_C (causal min,      Eq 41) = {d(IC)} bits   (paper states ~0.080)")
    assert abs(IC - mp.mpf('0.080')) < mp.mpf('0.002'), "I_C should be ~0.080"
    assert abs(IR - mp.mpf('0.046')) < mp.mpf('0.002'), "I_R should be ~0.046"

    print("\n  comparison (our witnesses vs H&B minima):")
    for name in names:
        I = mi[name]
        print(f"    {name}: I={d(I)}  |  vs I_R={d(IR)}: {'ABOVE' if I > IR else 'BELOW'} "
              f"(diff {float(I-IR):+.6f})  |  vs I_C={d(IC)}: {'ABOVE' if I > IC else 'BELOW'} (diff {float(I-IC):+.6f})")

    # causal factorization check: p(ab|lam) rank-1  <=>  p00*p11 == p01*p10 (exact Q2)
    print("\n[causal factorization] does the witness satisfy p(x,y|lam)=p(x|lam)p(y|lam)?")
    for name in names:
        w = W[name]; p = w["p"]
        viol = 0
        for lam in range(len(p[0])):
            lhs = qmul(p[0][lam], p[3][lam])     # p00 * p11
            rhs = qmul(p[1][lam], p[2][lam])     # p01 * p10
            if qcmp(lhs, rhs) != 0:
                viol += 1
        status = "FACTORIZED (causal)" if viol == 0 else f"NOT factorized at {viol}/{len(p[0])} lambda (retrocausal-only)"
        print(f"  {name}: {status}")

    # cross-check MI against the saved .npy witnesses (bisection-tolerance)
    try:
        import numpy as np
        base = _ROOT + "/witnesses"
        print("\n[cross-check] MI from saved .npy witnesses (float, bisection tolerance):")
        for tag, name in [("solN2", "N2-roundrobin"), ("solN3", "N3-roundrobin"),
                          ("solN4", "N4-signalling"), ("ns_solN4", "N4-NS")]:
            p = np.load(f"{base}/{tag}_p.npy").astype(float)
            cols = p.shape[1]
            plam = p.sum(axis=0) / 4.0
            def Hf(d):
                d = d[d > 0]
                return float(-np.sum(d * np.log2(d)))
            I = Hf(plam) - sum(Hf(p[ab]) for ab in range(4)) / 4.0
            print(f"  {tag}: I={I:.6f} bits   (exact-construction value {d(mi[name])}; "
                  f"diff {float(abs(I - mi[name])):.2e})")
    except Exception as e:
        print(f"  [cross-check skipped: {e}]")

    print("\n" + "=" * 78)
    print("ALL ASSERTIONS PASSED")
    print("=" * 78)
    return 0

if __name__ == "__main__":
    sys.exit(main())
