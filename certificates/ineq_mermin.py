#!/usr/bin/env python3
import os as _os; _ROOT = _os.path.normpath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # release: paths relative to repo root
"""
ineq_mermin.py — Divergent explorer #5: does the four-correlator saturation generalize?

Inequalities: Mermin (GHZ-parity) Bell inequalities for n = 3 and n = 5 parties,
2 settings each (X or Y), binary outcomes.

    M_n = sum_{s in {X,Y}^n, |S| even} (-1)^{|S|/2} E(s),   E(s) = prod_i o_i[s_i].

n=3:  M_3 = E(XXX) - E(XYY) - E(YXY) - E(YYX);  local bound 2, quantum value 4
      (GHZ_3 achieves the algebraic maximum via four perfect correlators).
      [Operator form per arXiv:1710.10802 (classical <= 2, quantum max 4);
       "four tripartite Mermin perfect correlators" per Hall PRA 84,022102 (2011),
       as cited in Alai arXiv:2608.00124.]
n=5:  M_5 = sum over 16 even-weight settings; local bound 4 (= 2^((n-1)/2)),
      quantum value 16 = algebraic maximum (GHZ_5).

Metric (exact generalization of the project's CHSH metric):
  N hidden states, uniform source rho = 1/N; conditional sources p[s, lam] >= 0,
  row-sum 1, for EVERY joint setting s in {X,Y}^n (only the even-weight rows enter
  M_n and are never biased at an optimum); deterministic responses o_i[t, lam] in
  {+/-1};  F = max_s TV(p(.|s) || rho),  TV = (1/2) sum_lam |p[s,lam] - 1/N|.

Exact finite-max formula (Theorem-1 analog; proof in review/divergent/ineq.md):
  S_max(N,F) = max_Q [ S0(Q) + sum_s g_s(F) ]
  Q = multiset of N distinct realizable patterns. A pattern is the even-weight
  row vector e_s(lam) = P * prod_{i: s_i=Y} r_i; the 2^(n+1) parameter pairs
  (P, r_1..r_n) have a 2-fold redundancy ((P,r) ~ (P,-r), invisible on
  even-weight rows), so there are 2^n distinct patterns.
  n_s(Q) = #{lam in Q : c(s) e_s = -1},  S0(Q) = (1/N) sum_lam C(pattern),
  g_s(F) = 2 min(F, n_s/N) if 0 < n_s < N, else 0.   (Lemma 1 applies: f in {+/-1}.)

All arithmetic exact (Fraction / integers). Exits 0 iff every assertion passes.
"""
import sys, time
from fractions import Fraction as Fr
from itertools import combinations_with_replacement

T0 = time.time()
NCHECK = 0
def ok(msg):
    global NCHECK; NCHECK += 1
    print(f"  [check {NCHECK:3d}] {msg}")

def gen_patterns(n):
    """All 2^n DISTINCT realizable Mermin-n patterns (even-weight row vectors).

    A deterministic strategy gives e_s = P * prod_{i: s_i=Y} r_i with
    (P, r_1..r_n) in {+/-1}^{n+1}; the 2^(n+1) parameter pairs carry a 2-fold
    redundancy, (P,r) ~ (P,-r), invisible on even-weight rows (global flip of r
    changes each even-weight product by (-1)^even = +1). We quotient it out by
    keeping only pairs with r_1 = +1 (bit 0 of rmask set). miss/C are invariant
    under the flip. Every distinct vector is realized by >= 1 strategy, and two
    strategies with the same vector are interchangeable for M_n, so this is the
    exact effective pattern space.
    """
    S = [m for m in range(1 << n) if bin(m).count("1") % 2 == 0]
    out = []
    for P in (1, -1):
        for rmask in range(1 << n):
            if not ((rmask >> 0) & 1):      # canonical representative: r_1 = +1
                continue
            r = tuple(1 if (rmask >> i) & 1 else -1 for i in range(n))
            miss = 0; Cval = 0
            for j, m in enumerate(S):
                v = P
                for i in range(n):
                    if (m >> i) & 1:
                        v *= r[i]
                c = 1 if (bin(m).count("1") // 2) % 2 == 0 else -1
                Cval += c * v
                if c * v == -1:
                    miss |= 1 << j
            out.append(dict(P=P, r=rmask, S=S, miss=miss, C=Cval))
    return out

def exhaustive(pats, N):
    """Exact S_max(N, k/N) for k=0..N by full multiset enumeration (integer inner loop).
    Also returns min over Q (with no fully-wrong row) of max_s n_s  -> F* numerator."""
    R = len(pats[0]["S"])
    missm = [p["miss"] for p in pats]; Cvals = [p["C"] for p in pats]
    Pn = len(pats)
    best_num = [-10**9] * (N + 1)
    min_maxns = None
    t_start = time.time()
    for combo in combinations_with_replacement(range(Pn), N):
        ns = [0] * R
        Csum = 0
        for p in combo:
            Csum += Cvals[p]
            m = missm[p]
            while m:
                lsb = m & -m
                ns[lsb.bit_length() - 1] += 1
                m ^= lsb
        mx = max(ns)
        if mx < N and (min_maxns is None or mx < min_maxns):
            min_maxns = mx
        for k in range(N + 1):
            t = 0
            for s in range(R):
                nsv = ns[s]
                if 0 < nsv < N:
                    t += k if k < nsv else nsv
            num = Csum + 2 * t
            if num > best_num[k]:
                best_num[k] = num
    return {k: Fr(best_num[k], N) for k in range(N + 1)}, min_maxns

# ---------------------------------------------------------------------------
print("=== PART A: n=3 pattern structure (exact) ===")
d3 = gen_patterns(3)
ok(f"n=3: {len(d3)} distinct realizable patterns (= 2^3; the 16 (P,r) pairs "
   f"carry a 2-fold redundancy (P,r)~(P,-r))")
assert len(d3) == 8
R3 = len(d3[0]["S"])
ok(f"n=3: {R3} active setting blocks (even-weight); algebraic max = {R3}")
miss_counts = sorted(bin(p["miss"]).count("1") for p in d3)
ok(f"n=3: min miss count = {min(miss_counts)} (parity obstruction: no 0-miss pattern)")
assert min(miss_counts) == 1
nm3 = [p for p in d3 if bin(p["miss"]).count("1") == 1]
ok(f"n=3: {len(nm3)} near-miss patterns (each misses exactly one condition)")
assert len(nm3) == 4
per_cond = {}
for p in nm3:
    j = (p["miss"].bit_length() - 1)
    per_cond.setdefault(j, []).append(p)
ok(f"n=3: near-miss patterns per condition = {[len(per_cond[j]) for j in range(R3)]}")
assert all(len(v) == 1 for v in per_cond.values()) and len(per_cond) == 4
B3 = max(p["C"] for p in d3)
ok(f"n=3: local bound B_3 = max pattern value = {B3} (Mermin: 2^((n-1)/2) = 2)")
assert B3 == 2
# quantum value on GHZ_3: <s> = i^{|S|} for even |S| -> M_3(GHZ) = sum c(s) i^{|S|}
q3 = sum((1 if bin(m).count("1") // 2 % 2 == 0 else -1) * (1 if bin(m).count("1") % 4 == 0 else -1)
         for m in d3[0]["S"])
ok(f"n=3: M_3(GHZ_3) = {q3} = algebraic max (quantum value = algebraic max)")
assert q3 == 4

# ---------------------------------------------------------------------------
print("\n=== PART B: n=3 exact curves S_max(N, k/N), N = 2..6 (full enumeration) ===")
res3 = {}
for N in range(2, 7):
    curve, mmns = exhaustive(d3, N)
    res3[N] = (curve, mmns)
    print(f"  N={N}: " + "  ".join(f"S({k}/{N})={curve[k]}" for k in range(N + 1)))
    # local bound at F=0
    assert curve[0] == Fr(2), f"N={N}: S_max(N,0) must be B_3=2"
    ok(f"n=3 N={N}: S_max(N,0) = 2 = B_3")
    # first segment: S_max(N, 1/N) = 2 + 2*min(N,4)/N
    assert curve[1] == Fr(2 * N + 2 * min(N, 4), N), f"N={N}: first-segment slope"
    ok(f"n=3 N={N}: S_max(N,1/N) = 2 + 2*min(N,4)/N  (first segment 2+2*min(N,4)*F)")
    # F* to reach algebraic max 4: min_Q max_s n_s / N  must equal ceil(N/4)/N
    cap_num = -(-N // 4)
    assert mmns == cap_num, f"N={N}: min_Q max n_s = {mmns} != ceil(N/4) = {cap_num}"
    ok(f"n=3 N={N}: F*(M_3=4) = min_Q max_s n_s/N = {cap_num}/{N} = ceil(N/4)/N")
    # curve consistency: value at cap breakpoint is exactly 4, below before
    assert curve[cap_num] == Fr(4)
    if cap_num > 0:
        assert curve[cap_num - 1] < Fr(4)
    ok(f"n=3 N={N}: S_max(N,{cap_num}/N) = 4 and S_max(N,{cap_num-1}/N) < 4")

# ---------------------------------------------------------------------------
print("\n=== PART C: n=3 general-N construction + explicit N=4 witness ===")
# round-robin over the 4 near-miss patterns: state i misses condition (i mod 4).
# Asserts max_s n_s <= ceil(N/4) and sum n_s = N for N<=64.
nm_list = [[] for _ in range(R3)]
for p in nm3:
    nm_list[p["miss"].bit_length() - 1].append(p)
assert all(len(v) == 1 for v in nm_list)
for N in range(2, 65):
    ns = [0] * R3
    for i in range(N):
        p = nm_list[i % 4][0]
        m = p["miss"]
        while m:
            lsb = m & -m
            ns[lsb.bit_length() - 1] += 1
            m ^= lsb
    cap_num = -(-N // 4)
    assert max(ns) <= cap_num and sum(ns) == N, f"construction failed at N={N}: {ns}"
ok("n=3: round-robin construction gives max_s n_s <= ceil(N/4), sum n_s = N for N = 2..64")
# Hence S_max(N, ceil(N/4)/N) >= 2 + 2*sum(min(F, n_s/N)) = 2 + 2*N/N = 4 exactly,
# and the per-pattern miss>=1 lower bound gives F* >= ceil(N/4)/N:  F*(n=3; N) = ceil(N/4)/N.

# explicit N=4 witness in Fractions: the round-robin Q_4 (one near-miss pattern
# per condition, n_s = 1 for all s); p[s, lam] = 1/3 on the 3 states correct at
# row s, per-row TV = 1/4, M_3 = 4.
states = [nm_list[j][0] for j in range(4)]
assert len(states) == 4
Fstar = Fr(1, 4)
M3val = Fr(0)
for sidx, m in enumerate(d3[0]["S"]):
    c = 1 if bin(m).count("1") // 2 % 2 == 0 else -1
    row = Fr(0); tv = Fr(0)
    for p in states:
        v = p["P"]
        for i in range(3):
            if (m >> i) & 1 and ((p["r"] >> i) & 1 == 0):
                v *= -1
        wgt = Fr(1, 3) if c * v == 1 else Fr(0)
        row += wgt * c * v
        tv += abs(wgt - Fr(1, 4))
    assert row == 1, f"row {sidx} value {row} != 1"
    assert Fr(1, 2) * tv == Fstar, f"row {sidx} TV {tv/2} != 1/4"
    M3val += row
ok(f"n=3 N=4 witness: all 4 rows at +1, per-row TV = 1/4 exactly, M_3 = {M3val}")
assert M3val == 4

# ---------------------------------------------------------------------------
print("\n=== PART D: n=5 pattern structure (exact) ===")
d5 = gen_patterns(5)
R5 = len(d5[0]["S"])
ok(f"n=5: {len(d5)} distinct realizable patterns; {R5} active setting blocks; "
   f"algebraic max = {R5}")
assert len(d5) == 32 and R5 == 16
mc5 = sorted(bin(p["miss"]).count("1") for p in d5)
m5min = mc5[0]
ok(f"n=5: min miss count = {m5min}  ('all-but-one' structure FAILS: m_5 > 1)")
assert m5min == 6
nm5 = [p for p in d5 if bin(p["miss"]).count("1") == m5min]
ok(f"n=5: {len(nm5)} min-miss patterns (of the 32 distinct patterns)")
assert len(nm5) == 16
B5 = max(p["C"] for p in d5)
ok(f"n=5: local bound B_5 = max pattern value = {B5} (Mermin: 2^((n-1)/2) = 4)")
assert B5 == 4
q5 = sum((1 if bin(m).count("1") // 2 % 2 == 0 else -1) * (1 if bin(m).count("1") % 4 == 0 else -1)
         for m in d5[0]["S"])
ok(f"n=5: M_5(GHZ_5) = {q5} = algebraic max (quantum value = algebraic max)")
assert q5 == 16

# ---------------------------------------------------------------------------
print("\n=== PART E: n=5 exact S_max(N, k/N), N = 2..4 (full enumeration) ===")
res5 = {}
for N in range(2, 5):
    curve, mmns = exhaustive(d5, N)
    res5[N] = (curve, mmns)
    print(f"  N={N}: " + "  ".join(f"S({k}/{N})={curve[k]}" for k in range(N + 1)))
    assert curve[0] == Fr(4), f"N={N}: S_max(N,0) must be B_5=4"
    ok(f"n=5 N={N}: S_max(N,0) = 4 = B_5")
    # full curve pinned to the exact enumeration (v55 §C/D). First segment is
    # 4 + 32/N for N >= 4 (slope 2*16: all 16 rows active on [0,1/N]); the
    # small-N values differ because not all rows can be jointly activated.
    EXPECT_CURVE = {2: [Fr(4), Fr(16), Fr(16)],
                    3: [Fr(4), Fr(12), Fr(16), Fr(16)],
                    4: [Fr(4), Fr(12), Fr(16), Fr(16), Fr(16)]}
    assert list(curve.values()) == EXPECT_CURVE[N], \
        f"N={N}: curve {list(curve.values())} != {EXPECT_CURVE[N]}"
    ok(f"n=5 N={N}: full curve pinned (S_max(N,1/N) = {EXPECT_CURVE[N][1]})")
    lb = -(-6 * N // 16)   # ceil(3N/8)
    assert mmns >= lb, f"N={N}: min max n_s = {mmns} < lower bound {lb}"
    # exact enumeration: the lower bound is TIGHT at N = 2,3,4 (residues 2,3,4
    # mod 16; only residues 0,1,2 have general-N constructions, Part F)
    assert mmns == lb, f"N={N}: min max n_s = {mmns} != ceil(3N/8) = {lb}"
    ok(f"n=5 N={N}: F*(M_5=16) = min_Q max n_s/N = {mmns}/{N} = ceil(3N/8)/N (tight)")

# ---------------------------------------------------------------------------
print("\n=== PART E2: A_max_B(N) — first-segment right slope at F=0 (exact) ===")
# On [0, 1/N] the finite-max formula is affine per Q: S0(Q) + 2F*A(Q), where
# A(Q) = #{s : 0 < n_s < N}. The right derivative at F=0 is therefore
# 2*A_max_B with A_max_B = max{A(Q) : S0(Q) = B_5} — the max over multisets
# whose patterns ALL achieve the local bound C = B_5 (the 16 min-miss
# patterns; each misses exactly m_min = 6 conditions, so sum_s n_s = 6N and
# A_max_B <= min(K, m_min*N)). Enumerated exactly over those 16 patterns.
cb_pats = [p for p in d5 if p["C"] == B5]
assert len(cb_pats) == 16, f"expected 16 bound-achieving patterns, got {len(cb_pats)}"
EXPECT_AMAXB = {2: 8, 3: 12, 4: 16, 5: 16, 6: 16}
# first-segment endpoint S_max(N, 1/N) from the pinned full curves (Part E
# here for N=2..4; ineq_mermin5_n5.out for N=5,6):
EXPECT_FIRST = {2: Fr(16), 3: Fr(12), 4: Fr(12), 5: Fr(52, 5), 6: Fr(28, 3)}
for N in range(2, 7):
    amaxb = 0
    for combo in combinations_with_replacement(range(len(cb_pats)), N):
        ns = [0] * R5
        for pidx in combo:
            m = cb_pats[pidx]["miss"]
            while m:
                lsb = m & -m
                ns[lsb.bit_length() - 1] += 1
                m ^= lsb
        a = sum(1 for x in ns if 0 < x < N)
        if a > amaxb:
            amaxb = a
    assert amaxb == EXPECT_AMAXB[N], \
        f"N={N}: A_max_B = {amaxb} != pinned {EXPECT_AMAXB[N]}"
    bound = min(R5, 6 * N)
    line_end = Fr(B5) + Fr(2 * amaxb, N)
    if N == 2:
        # chord (0,B_5)-(1/2, S(1/2)) is NOT supported by any bound-achieving Q:
        assert line_end < EXPECT_FIRST[N], \
            f"N=2: bound-achieving line {line_end} >= curve[1] {EXPECT_FIRST[N]}"
        ok(f"n=5 N=2: A_max_B = {amaxb} (bound min(K,6N) = {bound}, LOOSE); "
           f"first-segment slope 2*A_max_B = {2*amaxb}; bend in [0,1/2] "
           f"(chord to S(1/2)={EXPECT_FIRST[N]} needs A > A_max_B)")
    else:
        # single segment: the bound-achieving line reaches curve[1] exactly
        assert line_end == EXPECT_FIRST[N], \
            f"N={N}: line {line_end} != curve[1] {EXPECT_FIRST[N]}"
        ok(f"n=5 N={N}: A_max_B = {amaxb} (bound min(K,6N) = {bound}, "
           f"{'tight' if amaxb == bound else 'LOOSE'}); single segment "
           f"{B5}+{2*amaxb}*F on [0,1/N] (endpoint {line_end})")

# ---------------------------------------------------------------------------
print("\n=== PART F: n=5 general-N constructions (exact counting) ===")
# NOTE: an earlier draft claimed the 10 patterns (P=+1, wt(r)=2) were S_5-
# uniform with incidence 4 per condition. FALSE in two ways: that subfamily
# never misses condition 0 (the weight-0 string XXXXX), since e(XXXXX) = P = +1
# for every P=+1 pattern; and the true uniform family is the FULL set of 16
# DISTINCT min-miss patterns (incidence 6 per condition, asserted below).

# (F1) full 16 distinct min-miss family: uniform incidence 6 per condition.
inc = [0] * R5
for p in nm5:
    m = p["miss"]
    while m:
        lsb = m & -m
        inc[lsb.bit_length() - 1] += 1
        m ^= lsb
ok(f"n=5: full 16 min-miss family incidence per condition = {inc}")
assert all(x == 6 for x in inc) and sum(inc) == 16 * 6

# (F2) N = 16m: m copies of every distinct min-miss pattern -> n_s = 6m for
# ALL s, so max n_s = 6m; the lower bound ceil(3N/8)/N (min miss 6 => sum n_s
# >= 6N => max n_s >= ceil(6N/16); S=16 forces row s supported on its
# correct-sign states, TV_s >= n_s/N) then gives F*(M_5=16; 16m) = 3/8 exactly.
for m in range(1, 5):
    N = 16 * m
    assert 6 * m == -(-3 * N // 8), f"N={N}: 6m != ceil(3N/8)"
ok("n=5: F*(M_5=16; N=16m) = 3/8 exactly, m = 1..4 (construction + lower bound)")

# (F3) Remainder witnesses for r = 0,1,2 (hardcoded from exact search; each
# verified here by integer counting). Base: m copies of all 16 distinct
# min-miss patterns (incidence 6m per condition); add the witness multiset.
# If the witness keeps per-condition incidence <= ceil(3r/8), then
# F*(M_5=16; 16m+r) = ceil(3(16m+r)/8)/(16m+r) exactly (lower bound met).
REM = {0: [], 1: [0x71e7], 2: [0x2b42, 0xd4bd]}
# r=2 pair is the complementary flip pair: one pattern misses 6 conditions,
# its partner misses the other 10; each condition missed by exactly one.
assert REM[2][0] & REM[2][1] == 0 and REM[2][0] | REM[2][1] == (1 << R5) - 1
for r, ws in REM.items():
    cnt = [0] * R5
    for w in ws:
        for s in range(R5):
            if (w >> s) & 1:
                cnt[s] += 1
    cap = -(-3 * r // 8)
    assert max(cnt) <= cap, f"r={r}: witness incidence {max(cnt)} > ceil(3r/8) = {cap}"
ok("n=5: remainder witnesses r=0,1,2 keep per-condition incidence <= ceil(3r/8)")
for m in range(0, 3):
    for r in (0, 1, 2):
        N = 16 * m + r
        if N < 2:
            continue
        lb = -(-3 * N // 8)
        g = {0: 0, 1: 1, 2: 1}[r]   # max witness incidence for remainder r
        assert 6 * m + g == lb, f"N={N}: construction 6m+g != ceil(3N/8)"
ok("n=5: F*(M_5=16; N) = ceil(3N/8)/N exactly for all N ≡ 0,1,2 (mod 16), N >= 2")

# Coverage limit (stated in the report): for residues r not in {0,1,2} the
# exact value is unresolved; it is sandwiched between ceil(3N/8)/N and
# (6m + r)/(16m + r) with N = 16m + r (trivial remainder bound). In either
# case lim_{N->inf} F* = 3/8.

print(f"\nALL CHECKS PASSED ({NCHECK} checks).  elapsed {time.time()-T0:.1f}s")
sys.exit(0)
