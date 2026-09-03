#!/usr/bin/env python3
import os as _os; _ROOT = _os.path.normpath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # release: paths relative to repo root
"""
ineq_mermin5_n5.py — Divergent explorer #5 (supplement): Mermin n=5 exact
S_max(N, k/N) and F*(M_5 = 16; N) for N = 5 and N = 6 by FULL multiset
enumeration over C(32+N-1, N) multisets of the 32 distinct realizable
patterns (ineq_mermin.py Part E covered N = 2..4 only).

Purpose: test whether the PROVED lower bound F* >= ceil(3N/8)/N (min miss
count 6 => sum_s n_s >= 6N => max_s n_s >= ceil(6N/16) = ceil(3N/8); S = 16
forces each row fully bent, TV_s >= n_s/N) is TIGHT at residues r = 5 and
r = 6 (mod 16), which ineq_mermin.py left unresolved (its constructions
cover N ≡ 0,1,2 mod 16; exact enumeration covered N = 2,3,4).

Pattern generation is copied verbatim from ineq_mermin.py (gen_patterns);
the exhaustive loop mirrors ineq_mermin.py's `exhaustive` (integer inner
loop, Fraction only at output). All arithmetic exact. Exits 0 iff every
assertion passes.
"""
import sys, time
from fractions import Fraction as Fr
from itertools import combinations_with_replacement

T0 = time.time()
NCHECK = 0
def ok(msg):
    global NCHECK; NCHECK += 1
    print(f"  [check {NCHECK:3d}] {msg}")

# --- copied verbatim from ineq_mermin.py (gen_patterns) ----------------------
def gen_patterns(n):
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
    """Exact S_max(N, k/N) for k=0..N by full multiset enumeration."""
    R = len(pats[0]["S"])
    missm = [p["miss"] for p in pats]; Cvals = [p["C"] for p in pats]
    Pn = len(pats)
    best_num = [-10**9] * (N + 1)
    min_maxns = None
    ncombos = 0
    t_start = time.time()
    for combo in combinations_with_replacement(range(Pn), N):
        ncombos += 1
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
    return {k: Fr(best_num[k], N) for k in range(N + 1)}, min_maxns, ncombos

d5 = gen_patterns(5)
assert len(d5) == 32 and len(d5[0]["S"]) == 16
R5 = 16

print("=== Mermin n=5: exact S_max(N, k/N), N = 5, 6 (full enumeration) ===")
for N in (5, 6):
    t1 = time.time()
    curve, mmns, ncombos = exhaustive(d5, N)
    el = time.time() - t1
    print(f"  N={N}: " + "  ".join(f"S({k}/{N})={curve[k]}" for k in range(N + 1))
          + f"   [{ncombos} multisets, {el:.1f}s]")
    assert curve[0] == Fr(4), f"N={N}: S_max(N,0) must be B_5 = 4"
    ok(f"n=5 N={N}: S_max(N,0) = 4 = B_5")
    # full curve pinned to the exact enumeration (v55 §C/D); first segment is
    # 4 + 32/N = slope 2*16 (all 16 rows active on [0,1/N]) for N >= 4
    EXPECT_CURVE = {5: [Fr(4), Fr(52, 5), Fr(72, 5), Fr(16), Fr(16), Fr(16)],
                    6: [Fr(4), Fr(28, 3), Fr(14), Fr(16), Fr(16), Fr(16), Fr(16)]}
    assert list(curve.values()) == EXPECT_CURVE[N], \
        f"N={N}: curve {list(curve.values())} != {EXPECT_CURVE[N]}"
    ok(f"n=5 N={N}: full curve pinned (S_max(N,1/N) = 4 + 32/N = {Fr(4 * N + 32, N)})")
    lb = -(-3 * N // 8)   # ceil(3N/8), PROVED lower bound (ineq_mermin.py Part E logic)
    assert mmns >= lb, f"N={N}: min max n_s = {mmns} < lower bound {lb}"
    if mmns == lb:
        ok(f"n=5 N={N}: F*(M_5=16) = min_Q max n_s/N = {mmns}/{N} = ceil(3N/8)/N "
           f"— lower bound TIGHT at residue {N % 16}")
    else:
        ok(f"n=5 N={N}: F*(M_5=16) = min_Q max n_s/N = {mmns}/{N} > ceil(3N/8)/N "
           f"= {lb}/{N} — lower bound NOT tight at residue {N % 16}")
    # curve consistency: value at the F* breakpoint is exactly 16, below before
    assert curve[mmns] == Fr(16), f"N={N}: S_max(N, mmns/N) must be 16"
    if mmns > 0:
        assert curve[mmns - 1] < Fr(16), f"N={N}: curve must cross 16 at mmns/N"
    ok(f"n=5 N={N}: S_max(N,{mmns}/N) = 16 and S_max(N,{mmns-1}/N) < 16")

# --- Part B: audit of ineq_mermin.py Part F hardcoded remainder witnesses ----
# GAP CLOSURE: ineq_mermin.py Part F asserts the per-condition incidence of the
# hardcoded REM masks but never verifies that each mask is the miss set of a
# REALIZABLE pattern. This part closes that gap (exact, from gen_patterns).
print("\n=== Part B: REM witness audit (realizability + incidence) ===")
REM = {0: [], 1: [0x71e7], 2: [0x2b42, 0xd4bd]}   # copied from ineq_mermin.py Part F
allmasks = set(p["miss"] for p in d5)
for r, ws in REM.items():
    for w in ws:
        assert w in allmasks, f"r={r}: mask {w:#06x} is NOT a realizable miss mask"
    ok(f"n=5 r={r}: all {len(ws)} hardcoded witness mask(s) are realizable pattern miss masks")
assert REM[2][0] & REM[2][1] == 0 and REM[2][0] | REM[2][1] == (1 << R5) - 1
ok("n=5 r=2: witness pair disjoint and complementary (each condition missed by exactly one)")
for r, ws in REM.items():
    cnt = [0] * R5
    for w in ws:
        for s in range(R5):
            if (w >> s) & 1:
                cnt[s] += 1
    cap = -(-3 * r // 8)
    assert max(cnt) <= cap, f"r={r}: witness incidence {max(cnt)} > ceil(3r/8) = {cap}"
ok("n=5: per-condition witness incidence <= ceil(3r/8) for r = 0,1,2 (Part F claim re-verified)")
# the r=2 construction is not a single hardcoded accident: count ALL disjoint
# realizable mask pairs (each gives a valid remainder-2 witness pair).
mlist = sorted(allmasks)
npairs = sum(1 for i in range(len(mlist)) for b in mlist[i + 1:] if mlist[i] & b == 0)
ok(f"n=5: {npairs} disjoint realizable miss-mask pairs exist (r=2 construction family)")
assert npairs >= 1

print(f"\nALL CHECKS PASSED ({NCHECK} checks).  elapsed {time.time()-T0:.1f}s")
sys.exit(0)
