#!/usr/bin/env python3
import os as _os; _ROOT = _os.path.normpath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # release: paths relative to repo root
"""Divergent explorer #2: exact S_max(N,F) curves for N = 9..14 + general-N closed form.

Reuses the structure of audit1.py Part A (exact Fraction enumeration over all
C(8+N-1, N) even-parity pattern multisets; per-cell affine lines and upper
envelopes). No floating point except the final F* cross-check (Part C).

What this script certifies (exit 0 iff every assertion passes):

PART A  For each N in 9..14: exhaustive exact curve.
        Key certificate: W_k := max_Q Smax_Q(k/N) equals the round-robin closed
        form cf(N, k/N) = 2 + sum_ab 2*min(k/N, n*_ab/N) at EVERY breakpoint
        k = 0..N. This certifies the whole piecewise-linear curve, not just the
        breakpoints: on each cell [k/N,(k+1)/N] both the envelope E(F) = max_Q
        line_Q(F) and cf(F) are affine in F (line_Q is the cell-k affine form of
        Smax_Q), Q_N's cell-k line IS cf on that cell, so D(F) := E(F)-cf(F) is
        convex with D(k/N)=D((k+1)/N)=0, hence D <= 0 on the whole cell; and
        E >= cf pointwise since Q_N is among the maximized multisets. So
        E = cf on every cell exactly.
        Also: explicit task assertions (Q_N lower bound == envelope on [0,1/N]
        and at the cap ceil(N/4)/N), per-cell winner census (which multisets
        attain the top line: balanced all-positive / non-balanced positive /
        with negative patterns), and legality checks for the exchange lemma.

PART B  Bend-structure table for N = 2..14 in terms of N mod 4, from the same
        exhaustive data (N=2..8 re-computed here; sizes <= 6435 multisets).

PART C  Float cross-check at F*=(sqrt(2)-1)/4: max_Q Smax_Q(F*) = 2*sqrt(2)
        for N = 9..14 (analog of audit1.py Part B).

PART D  Brute-force certificates (exact, all breakpoints k/N) of the two lemmas
        used in the general-N proof:
          E: exchange lemma — for all nC,nD in 1..N and all k=0..N:
             4/N + gk(nC-1)-gk(nC) + gk(nD-1)-gk(nD) >= 0,
             where gk(x) := g(F=k/N)(x) with the Lemma-1 edge convention
             (g = 0 at x in {0,N}); this is exactly the change of Smax_Q when a
             negative pattern -M_A is replaced by M_B (B != A): +4/N on S0 and
             -1 on exactly two wrong counts.
          B: balancing lemma — for all compositions p of N into 4 parts, all
             ordered pairs with p_i >= p_j+2, all k=0..N:
             gk(p_i)+gk(p_j) <= gk(p_i-1)+gk(p_j+1).
        Together with the (finite, N-independent) combinatorial argument in
        review/divergent/curves.md these certify the general-N theorem for all
        N; the brute force cross-checks the lemmas at every N = 2..14.

Output: verification/logs/curves_highN.out
"""
import itertools
import math
import time
from fractions import Fraction as Fr

LOG = _ROOT + "/logs/"
out = open(LOG + "curves_highN.out", "w")


def pr(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    out.write(s + "\n")

T_START = time.time()
WALL_LIMIT = 900.0  # hard guard: state what finished if exceeded (task: <~15 min/N)

SIGMA = (1, 1, 1, -1)  # ab order 00, 01, 10, 11
PATTERNS = [q for q in itertools.product((1, -1), repeat=4)
            if q[0] * q[1] * q[2] * q[3] == 1]
assert len(PATTERNS) == 8


def C(q):
    """CHSH value of one deterministic vertex (pattern)."""
    return q[0] + q[1] + q[2] - q[3]


CVAL = [C(q) for q in PATTERNS]
NEGIDX = [i for i in range(8) if CVAL[i] == -2]
POSIDX = [i for i in range(8) if CVAL[i] == 2]
assert len(NEGIDX) == 4 and len(POSIDX) == 4
# wrong-count contribution of each pattern (row ab: 1 iff q_ab != sigma_ab)
WRONG = [tuple(1 if q[ab] != SIGMA[ab] else 0 for ab in range(4)) for q in PATTERNS]


def count_vectors(N):
    """All 8-tuples of nonnegative integers summing to N (pattern multiplicities).

    Count is C(N+7,7); each tuple is a multiset Q of N even-parity patterns.
    """
    def rec(rem, d, cur):
        if d == 1:
            yield cur + (rem,)
        else:
            for x in range(rem + 1):
                for rest in rec(rem - x, d - 1, cur + (x,)):
                    yield rest
    yield from rec(N, 8, ())


def compositions4(N):
    """All 4-tuples of nonnegative integers summing to N (C(N+3,3) of them)."""
    for a in range(N + 1):
        for b in range(N - a + 1):
            for c in range(N - a - b + 1):
                yield (a, b, c, N - a - b - c)


def gk(k, N, x):
    """g(F=k/N)(x) = 2*min(k/N, x/N) for 0 < x < N; 0 at x in {0,N} (Lemma 1)."""
    if x <= 0 or x >= N:
        return Fr(0)
    return Fr(2 * min(k, x), N)


def nstar_of(N):
    """Round-robin count vector: state i takes miss-type M_{i mod 4}.

    Class j gets k+1 states iff j < r where N = 4k + r. (The closed form is
    symmetric in the assignment, so any permutation of these entries works.)
    """
    k, r = divmod(N, 4)
    return tuple(k + 1 if j < r else k for j in range(4))


def cf_value(N, F):
    """Conjectured (then proved) closed form: S_max(N,F) = 2 + sum_ab 2 min(F, n*/N)."""
    ns = nstar_of(N)
    return Fr(2) + sum(Fr(2 * min(F, Fr(x, N))) for x in ns)


def cf_pieces(N):
    """Piecewise form of cf from the k,r structure: list of (lo, hi, intercept, slope).

    N = 4k+r: [0,k/N] slope 8 (or 2*min(N,4)-ish for N<4); if r>0 a middle
    segment [k/N,(k+1)/N] with slope 2r and intercept 2 + 2k(4-r)/N; then 4.
    Derived independently of cf_value (cross-check target).
    """
    k, r = divmod(N, 4)
    ns = nstar_of(N)
    # generic derivation: slope on cell m is 2*#{j: m < n*_j < N}, intercept
    # 2 + sum_{j: 1<=n*_j<=m} 2 n*_j/N
    pieces = []
    for m in range(N):
        sl = 2 * sum(1 for x in ns if m < x < N)
        itc = Fr(2) + sum(Fr(2 * x, N) for x in ns if 1 <= x <= m)
        pieces.append((Fr(m, N), Fr(m + 1, N), itc, sl))
    # merge adjacent cells with identical (intercept, slope)
    merged = []
    for lo, hi, itc, sl in pieces:
        if merged and merged[-1][2] == itc and merged[-1][3] == sl:
            a, b, i2, s2 = merged.pop()
            merged.append((a, hi, itc, sl))
        else:
            merged.append((lo, hi, itc, sl))
    return merged


def analyze(N, detailed):
    """Exhaustive exact analysis for one N. Returns summary dict."""
    t0 = time.time()
    k, r = divmod(N, 4)
    ns = nstar_of(N)
    capnum = -(-N // 4)  # ceil(N/4)
    cap = Fr(capnum, N)

    W = [None] * (N + 1)       # max over Q of Smax_Q(k/N)
    Wneg = [None] * (N + 1)    # same, restricted to multisets with >=1 negative pattern
    # per-cell winner census: line of Q on cell m equals the top (cf) line
    tie_bal = [0] * N          # all-positive, balanced counts
    tie_nb = [0] * N           # all-positive, non-balanced
    tie_neg = [0] * N          # contains negative patterns
    legality_bad = 0
    bestF = 0.0                # float F* cross-check accumulator

    # cf cell lines (top line on each cell)
    m_cf = [0] * N
    b_cf = [Fr(0)] * N
    for m in range(N):
        sl = 2 * sum(1 for x in ns if m < x < N)
        itc = Fr(2) + sum(Fr(2 * x, N) for x in ns if 1 <= x <= m)
        m_cf[m] = sl
        b_cf[m] = itc

    Fstar = (math.sqrt(2) - 1) / 4.0
    total = 0
    for counts in count_vectors(N):
        total += 1
        snum = sum(CVAL[i] * counts[i] for i in range(8))   # S0 = snum/N
        T = sum(counts[i] for i in NEGIDX)                  # # negative patterns
        n = [0, 0, 0, 0]
        for i in range(8):
            c = counts[i]
            if c:
                w = WRONG[i]
                n[0] += w[0] * c
                n[1] += w[1] * c
                n[2] += w[2] * c
                n[3] += w[3] * c
        S0 = Fr(snum, N)

        # breakpoint values V_k = Smax_Q(k/N), k = 0..N  (one Fr add per k)
        for kk in range(N + 1):
            s = 0
            for j in range(4):
                nj = n[j]
                if 0 < nj < N:
                    s += min(kk, nj)
            v = S0 + Fr(2 * s, N)
            if W[kk] is None or v > W[kk]:
                W[kk] = v
            if T > 0 and (Wneg[kk] is None or v > Wneg[kk]):
                Wneg[kk] = v

        # cell lines + winner census + exchange-legality checks
        for m in range(N):
            sl = 0
            itc = S0
            for j in range(4):
                nj = n[j]
                if 1 <= nj <= m:
                    itc += Fr(2 * nj, N)
                elif m < nj < N:
                    sl += 2
            if sl == m_cf[m] and itc == b_cf[m]:
                if T == 0:
                    if sorted(n) == sorted(ns):
                        tie_bal[m] += 1
                    else:
                        tie_nb[m] += 1
                else:
                    tie_neg[m] += 1
        if T > 0:
            # legality of the exchange -M_A -> M_B (B != A): for each negative
            # type A present and each B != A, the two decreased rows C,D must
            # have n >= 1 (proved automatic in the report; certified here).
            for A in range(4):
                if counts[NEGIDX[A]] == 0:
                    continue
                for B in range(4):
                    if B == A:
                        continue
                    C_ = D_ = None
                    for j in range(4):
                        if j not in (A, B):
                            if C_ is None:
                                C_ = j
                            else:
                                D_ = j
                    if n[C_] < 1 or n[D_] < 1:
                        legality_bad += 1

        # float F* cross-check (independent evaluation, audit1 Part B style)
        s0f = float(snum) / N
        tot = s0f
        for j in range(4):
            nj = n[j]
            if 0 < nj < N:
                tot += 2.0 * min(Fstar, nj / N)
        if tot > bestF:
            bestF = tot

    assert total == math.comb(N + 7, 7), f"multiset count mismatch N={N}"
    assert legality_bad == 0, f"exchange-legality violation N={N}: {legality_bad}"

    # ---- the full-curve certificate: W_k == cf(k/N) for all k = 0..N --------
    ok = True
    for kk in range(N + 1):
        expect = cf_value(N, Fr(kk, N))
        if W[kk] != expect:
            ok = False
            pr(f"  MISMATCH N={N} k={kk}: envelope {W[kk]} != closed form {expect}")
    assert ok, f"closed-form curve mismatch at N={N}"

    # cross-check cf_value against the independently derived piecewise form
    pieces = cf_pieces(N)
    for lo, hi, itc, sl in pieces:
        assert itc + sl * lo == cf_value(N, lo), "piecewise/cf disagree (lo)"
        assert itc + sl * hi == cf_value(N, hi), "piecewise/cf disagree (hi)"

    # explicit task assertions: Q_N matches envelope on [0,1/N] and at the cap
    qn_first = (m_cf[0] == 2 * min(N, 4)) and (b_cf[0] == 2) \
        and (W[0] == 2) and (W[1] == cf_value(N, Fr(1, N)))
    assert qn_first, f"Q_N first-segment assertion failed N={N}"
    assert W[capnum] == 4, f"cap assertion failed N={N}: W[{capnum}] = {W[capnum]}"

    # negative-pattern tie probe at breakpoints (expect: strictly below top)
    neg_ties_bp = [kk for kk in range(N + 1)
                   if Wneg[kk] is not None and Wneg[kk] >= W[kk]]

    dt = time.time() - t0
    res = dict(N=N, k=k, r=r, ns=ns, cap=cap, W=W, pieces=pieces,
               tie_bal=tie_bal, tie_nb=tie_nb, tie_neg=tie_neg,
               neg_ties_bp=neg_ties_bp, bestF=bestF, dt=dt)

    if detailed:
        pr(f"--- N={N}: {total} pattern multisets (C({N+7},7)); {dt:.1f}s ---")
        pr(f"  round-robin counts n* = {ns};  cap F(S=4;N) = ceil({N}/4)/{N} = {cap}")
        pr("  S_max at F=k/N: " + ", ".join(str(v) for v in W))
        pr("  exact curve (certified equal to the exhaustive envelope on every cell):")
        for lo, hi, itc, sl in pieces:
            pr(f"    F in [{lo}, {hi}]:   S_max = {itc} + {sl}*F")
        pr("  winner census per cell (multisets attaining the top line on that cell):")
        for m in range(N):
            pr(f"    cell {m} [{Fr(m, N)}, {Fr(m+1, N)}]:  total="
               f"{tie_bal[m]+tie_nb[m]+tie_neg[m]}  balanced-pos={tie_bal[m]}  "
               f"nonbal-pos={tie_nb[m]}  with-neg={tie_neg[m]}")
        pr(f"  negative-pattern multisets tying the envelope at a breakpoint: "
           f"{neg_ties_bp if neg_ties_bp else 'none'}")
        pr(f"  Q_N lower bound == envelope on [0,1/N]: {qn_first};  "
           f"S_max at cap F={cap}: {W[capnum]} (expect 4)")

    return res


# ===========================================================================
pr("=== PART A: exact curves S_max(N,F) for N=9..14 (Fraction, exhaustive) ===")
pr("(structure reused from audit1.py Part A; full-curve certificate = breakpoint")
pr(" equality W_k == cf(k/N) at all k=0..N + convexity argument, see docstring)")
resA = {}
for N in range(9, 15):
    if time.time() - T_START > WALL_LIMIT:
        pr(f"  !! wall-clock guard hit before N={N}; finished: "
           f"{sorted(resA.keys())}")
        break
    resA[N] = analyze(N, detailed=True)

# ===========================================================================
pr("")
pr("=== PART B: bend-structure table N=2..14 (by N mod 4) ===")
resB = {}
for N in range(2, 9):
    resB[N] = analyze(N, detailed=False)
allres = dict(resB)
allres.update(resA)
pr("N  k r  n*            cap      segments (exact)")
for N in range(2, 15):
    R = allres[N]
    segs = " ; ".join(f"[{lo},{hi}] {itc}+{sl}F" for lo, hi, itc, sl in R["pieces"])
    pr(f"{N:<3}{R['k']:<2}{R['r']:<2}  {str(R['ns']):<14} {str(R['cap']):<8}  {segs}")
pr("")
pr("winner census summary (per cell: total / balanced-pos / nonbal-pos / with-neg):")
for N in range(2, 15):
    R = allres[N]
    cells = []
    for m in range(N):
        t = R["tie_bal"][m] + R["tie_nb"][m] + R["tie_neg"][m]
        cells.append(f"c{m}:{t}/{R['tie_bal'][m]}/{R['tie_nb'][m]}/{R['tie_neg'][m]}")
    pr(f"  N={N:<2}: " + "  ".join(cells))
# sub-cap structure assertions: on every cell fully contained in [0, cap] no
# negative-pattern multiset attains the top line; from the bend cell (k/N) up to
# the cap only balanced all-positive multisets attain it.
subcap_ok = True
for N in range(2, 15):
    R = allres[N]
    k, r = R["k"], R["r"]
    last_subcap_cell = -(-N // 4) - 1 if r == 0 else k   # cell index ending at cap
    for m in range(last_subcap_cell + 1):
        if R["tie_neg"][m] != 0:
            subcap_ok = False
            pr(f"  VIOLATION N={N} cell {m}: negative-pattern winner below cap")
    # bend cell (first cell with slope < 8 for r>0) must be balanced-only up to cap
    if r > 0:
        for m in range(k, last_subcap_cell + 1):
            if R["tie_nb"][m] != 0 or R["tie_neg"][m] != 0:
                subcap_ok = False
                pr(f"  VIOLATION N={N} cell {m}: non-balanced/negative winner "
                   f"on bend segment")
    # single-bend structure from the certified pieces
    segs = R["pieces"]
    slopes = [sl for _, _, _, sl in segs]
    if k == 0:      # N=2,3: one rising segment then the cap (no bend)
        assert slopes == [2 * min(N, 4), 0], f"N={N}: {slopes}"
    elif r == 0:    # N = 4k: 2+8F up to F=1/4, then 4 (no bend)
        assert slopes == [8, 0], f"N={N}: expected no bend, got {slopes}"
    else:           # N = 4k+r, k>=1, r>0: single bend at F=k/N, slope 2r after
        assert slopes == [8, 2 * r, 0], f"N={N}: bend structure {slopes}"
        lo_bend, hi_bend = segs[1][0], segs[1][1]
        assert (lo_bend, hi_bend) == (Fr(k, N), Fr(k + 1, N)), \
            f"N={N}: bend at wrong place {(lo_bend, hi_bend)}"
assert subcap_ok, "sub-cap winner-structure assertion failed"
pr("  sub-cap structure: for every N=2..14, no negative-pattern multiset attains")
pr("  S_max on any cell inside [0, cap]; from the bend cell k/N to the cap only")
pr("  balanced all-positive multisets attain it.  PASS")
neg_bp = {N: R["neg_ties_bp"] for N, R in allres.items() if R["neg_ties_bp"]}
pr(f"  (negative-pattern ties occur at breakpoints k/N >= ceil(N/4), i.e. at or")
pr(f"   beyond the cap, where S_max = 4: " +
   ", ".join(f"N={N}: k>={min(v)}" for N, v in sorted(neg_bp.items())) + ")")

# ===========================================================================
pr("")
pr("=== PART C: float cross-check at F*=(sqrt(2)-1)/4 for N=9..14 ===")
SQ2 = math.sqrt(2)
for N in range(9, 15):
    if N not in resA:
        continue
    b = resA[N]["bestF"]
    pr(f"  N={N}: max_Q Smax_Q(F*) = {b:.12f}   (2*sqrt2 = {2*SQ2:.12f}; "
       f"|diff| = {abs(b-2*SQ2):.2e})")
    assert abs(b - 2 * SQ2) < 1e-9, f"F* cross-check failed N={N}"

# ===========================================================================
pr("")
pr("=== PART D: brute-force certificates of the two proof lemmas (N=2..14) ===")
pr("Lemma E (exchange / negative removal): min over nC,nD in 1..N, k in 0..N of")
pr("  4/N + gk(nC-1)-gk(nC) + gk(nD-1)-gk(nD)   [must be >= 0 for all F]")
for N in range(2, 15):
    worst = Fr(0)
    worstarg = None
    for nC in range(1, N + 1):
        for nD in range(1, N + 1):
            for kk in range(N + 1):
                d = (Fr(4, N) + gk(kk, N, nC - 1) - gk(kk, N, nC)
                     + gk(kk, N, nD - 1) - gk(kk, N, nD))
                if d < worst:
                    worst = d
                    worstarg = (nC, nD, kk)
    assert worst >= 0, f"exchange lemma violated N={N} at {worstarg}: {worst}"
    pr(f"  N={N:<2}: min = {worst}   (at {worstarg})   PASS")
pr("Lemma B (balancing): for all compositions p of N into 4 parts, p_i >= p_j+2,")
pr("  gk(p_i)+gk(p_j) <= gk(p_i-1)+gk(p_j+1) for all k=0..N   [must hold]")
for N in range(2, 15):
    checked = 0
    bad = 0
    for p in compositions4(N):
        for i in range(4):
            for j in range(4):
                if i != j and p[i] >= p[j] + 2:
                    for kk in range(N + 1):
                        checked += 1
                        if (gk(kk, N, p[i]) + gk(kk, N, p[j])
                                > gk(kk, N, p[i] - 1) + gk(kk, N, p[j] + 1)):
                            bad += 1
    assert bad == 0, f"balancing lemma violated N={N}: {bad} counterexamples"
    pr(f"  N={N:<2}: {checked} checks, {bad} violations   PASS")

# ===========================================================================
pr("")
pr("=== SUMMARY ===")
pr(f"All assertions passed. Exact curves certified for N=9..14 (and re-certified")
pr(f"N=2..8) equal to the round-robin closed form S_max(N,F) = 2 + sum_ab 2*min(F, n*_ab/N).")
pr(f"For every N=2..14 no negative-pattern multiset attains the top line on any cell")
pr(f"inside [0, cap] (they tie only at/above the cap, where S_max = 4); from the bend")
pr(f"cell k/N to the cap only balanced all-positive multisets attain it.")
pr(f"Total wall time: {time.time()-T_START:.1f}s")
out.close()
print("curves_highN.out written to", LOG + "curves_highN.out")
