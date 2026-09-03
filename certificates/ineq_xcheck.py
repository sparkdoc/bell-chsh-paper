#!/usr/bin/env python3
import os as _os; _ROOT = _os.path.normpath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # release: paths relative to repo root
"""
ineq_xcheck.py — Divergent explorer #5: INDEPENDENT cross-check of ineq_mermin.py
and ineq_cglmp.py by a different method.

(1) Mermin n=3, n=5: brute force over ALL 4^n deterministic strategies directly
    (outcomes o_i[X], o_i[Y] in {+/-1} per party), WITHOUT the (P,r) pattern
    parameterization used in ineq_mermin.py. Checks: local bound = max M_n,
    number of distinct correlator vectors, min miss count, near-miss counts.
(2) GHZ quantum values: explicit exact complex-arithmetic (Q(i), Fraction-based)
    computation of <GHZ_n | P_1 x ... x P_n | GHZ_n> for every even-weight Pauli
    string with X/Y entries — no cos(m*pi/2) shortcut. Checks M_n(GHZ_n) = 2^(n-1)
    (algebraic maximum) for n=3,5.
(3) CGLMP d=3: brute force over all 81 deterministic strategies (A1,A2,B1,B2 in
    Z_3), WITHOUT the difference-pattern reduction of ineq_cglmp.py. Checks:
    local bound 2; modular constraint; 4 near-misses (one per row), each value 2.

All arithmetic exact (integers / Fraction). Exits 0 iff every assertion passes.
"""
import sys, time
from fractions import Fraction as Fr
from itertools import product

T0 = time.time()
NCHECK = 0
def ok(msg):
    global NCHECK; NCHECK += 1
    print(f"  [check {NCHECK:3d}] {msg}")

# ---------------------------------------------------------------------------
print("=== (1) Mermin n=3 and n=5 by direct strategy brute force ===")
def mermin_direct(n):
    """All 4^n deterministic strategies; return list of (value, missmask, evector)."""
    S = [m for m in range(1 << n) if bin(m).count("1") % 2 == 0]
    out = []
    # strategy: per party i, outcomes (oX_i, oY_i) in {+1,-1}^2
    for combo in product([(-1, 1), (1, -1), (-1, -1), (1, 1)], repeat=n):
        ev = []
        val = 0; miss = 0
        for j, m in enumerate(S):
            v = 1
            for i in range(n):
                o = combo[i][0] if ((m >> i) & 1) == 0 else combo[i][1]
                v *= o
            c = 1 if (bin(m).count("1") // 2) % 2 == 0 else -1
            ev.append(c * v)
            val += c * v
            if c * v == -1:
                miss |= 1 << j
        out.append((val, miss, tuple(ev)))
    return S, out

for n in (3, 5):
    S, rows = mermin_direct(n)
    R = len(S)
    B = max(r[0] for r in rows)
    evs = {r[2] for r in rows}
    missmin = min(bin(r[1]).count("1") for r in rows)
    nm = [r for r in rows if bin(r[1]).count("1") == missmin]
    per_cond = [sum(1 for r in nm if (r[1] >> j) & 1) for j in range(R)]
    print(f"  n={n}: strategies={len(rows)} distinct evectors={len(evs)} "
          f"local bound={B} min miss={missmin} #min-miss={len(nm)} per-cond={per_cond}")
    if n == 3:
        assert len(rows) == 64 and len(evs) == 8, "n=3: pattern count mismatch"
        assert B == 2, "n=3: local bound != 2"
        assert missmin == 1 and len(nm) == 32, "n=3: near-miss structure mismatch"
        assert all(x == 8 for x in per_cond), "n=3: not 8 near-miss strategies per condition"
        ok("n=3 direct brute force: B=2, 8 distinct patterns (64 strategies), min miss 1, "
           "one near-miss pattern per condition — matches ineq_mermin.py Part A")
    else:
        assert len(rows) == 1024 and len(evs) == 32, "n=5: pattern count mismatch"
        assert B == 4, "n=5: local bound != 4"
        assert missmin == 6 and len(nm) == 16 * 32, "n=5: min-miss structure mismatch"
        ok("n=5 direct brute force: B=4, 32 distinct patterns (1024 strategies), min miss 6, "
           "16 min-miss patterns — matches ineq_mermin.py Part D")

# ---------------------------------------------------------------------------
print("\n=== (2) GHZ quantum values by exact Q(i) matrix computation ===")
# complex numbers as (re, im) pairs of Fractions
def cadd(a, b): return (a[0] + b[0], a[1] + b[1])
def cmul(a, b): return (a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def cconj(a): return (a[0], -a[1])
X = ((Fr(0), Fr(0)), (Fr(1), Fr(0))), ((Fr(1), Fr(0)), (Fr(0), Fr(0)))
Y = ((Fr(0), Fr(0)), (Fr(0), Fr(-1))), ((Fr(0), Fr(1)), (Fr(0), Fr(0)))
Z = ((Fr(1), Fr(0)), (Fr(0), Fr(0))), ((Fr(0), Fr(0)), (-Fr(1), Fr(0)))

def sum_tuple(ts):
    acc = (Fr(0), Fr(0))
    for t in ts:
        acc = cadd(acc, t)
    return acc

def kron(A, B):
    rA, cA = len(A), len(A[0]); rB, cB = len(B), len(B[0])
    out = [[(Fr(0), Fr(0)) for _ in range(cA * cB)] for _ in range(rA * rB)]
    for i in range(rA):
        for k2 in range(cA):          # A's column index (block-column)
            if A[i][k2] == (Fr(0), Fr(0)):
                continue
            for k in range(rB):       # B's row index (row offset within block)
                rowB = B[k]
                for j in range(cB):   # B's column index
                    out[i*rB + k][k2*cB + j] = cadd(out[i*rB + k][k2*cB + j],
                                                    cmul(A[i][k2], rowB[j]))
    return tuple(tuple(row) for row in out)

def matvec(M, v):
    r = len(M); c = len(M[0])
    assert c == len(v)
    return tuple(sum_tuple([cmul(M[i][j], v[j]) for j in range(c)]) for i in range(r))

def inner(a, b):
    return sum_tuple([cmul(cconj(x), y) for x, y in zip(a, b)])

def ghz(n):
    """UNNORMALIZED GHZ vector |0..0> + |1..1> (norm^2 = 2); expectation values
    are computed as <v|P|v> / <v|v>, staying in Q(i)."""
    dim = 1 << n
    v = [(Fr(0), Fr(0))] * dim
    v[0] = (Fr(1), Fr(0)); v[dim-1] = (Fr(1), Fr(0))
    return tuple(v)

def pauli_string(n, m):
    """Tensor product over parties: bit i of m selects Y if set else X."""
    M = (((Fr(1), Fr(0)),),)  # 1x1 identity (kron unit)
    for i in range(n):
        M = kron(M, Y if ((m >> i) & 1) else X)
    return M

for n in (3, 5):
    S = [m for m in range(1 << n) if bin(m).count("1") % 2 == 0]
    psi = ghz(n)
    norm2 = inner(psi, psi)
    assert norm2 == (Fr(2), Fr(0)), "GHZ norm^2 must be 2"
    total = Fr(0)
    vals = {}
    for m in S:
        P = pauli_string(n, m)
        e = inner(psi, matvec(P, psi))
        assert e[1] == 0, "correlator not real"
        e = (e[0] / norm2[0], Fr(0))   # divide by norm^2 = 2
        c = 1 if (bin(m).count("1") // 2) % 2 == 0 else -1
        vals[m] = e[0]
        total += c * e[0]
    ok(f"n={n}: M_n(GHZ_n) = {total} = 2^(n-1) = algebraic max (exact Q(i) computation)")
    assert total == Fr(1 << (n - 1))
    # and each correlator is exactly +/-1 with sign (-1)^{|S|/2}
    for m in S:
        expect = Fr(1) if bin(m).count("1") % 4 == 0 else Fr(-1)
        assert vals[m] == expect, f"n={n} s={m:05b}: E(s) = {vals[m]} != expected {expect}"
    ok(f"n={n}: every even-weight correlator on GHZ_n is exactly (-1)^{{|S|/2}} (no shortcut used)")

# ---------------------------------------------------------------------------
print("\n=== (3) CGLMP d=3 by direct strategy brute force ===")
W = [
    lambda d: 1 if d == 0 else (-1 if d == 2 else 0),   # row (A1,B1)
    lambda d: 1 if d == 0 else (-1 if d == 1 else 0),   # row (A1,B2)
    lambda d: 1 if d == 2 else (-1 if d == 0 else 0),   # row (A2,B1)
    lambda d: 1 if d == 0 else (-1 if d == 2 else 0),   # row (A2,B2)
]
rows = []
for A1, A2, B1, B2 in product(range(3), repeat=4):
    d = ((A1 - B1) % 3, (A1 - B2) % 3, (A2 - B1) % 3, (A2 - B2) % 3)
    # modular constraint check (should hold for every strategy)
    assert (d[0] + d[3]) % 3 == (d[1] + d[2]) % 3
    w = tuple(W[r](d[r]) for r in range(4))
    rows.append((sum(w), d, w))
B = max(r[0] for r in rows)
ok(f"CGLMP: {len(rows)} strategies all satisfy d00+d11 == d01+d10 (mod 3); local bound B = {B}")
assert len(rows) == 81 and B == 2
distinct_d = {r[1] for r in rows}
ok(f"CGLMP: {len(distinct_d)} distinct difference patterns (27 expected)")
assert len(distinct_d) == 27
# algebraic-max target unattainable
assert (0, 0, 2, 0) not in distinct_d
ok("CGLMP: algebraic-max target (d00,d01,d10,d11)=(0,0,2,0) is not a valid pattern")
# near-misses: block exactly one row (w_r < 1), value 2, one per row
nm = [(r[1], [q for q in range(4) if r[2][q] < 1]) for r in rows if sum(r[2]) == B]
nm_uniq = {}
for d, b in nm:
    if len(b) == 1:
        nm_uniq.setdefault(b[0], set()).add(d)
ok(f"CGLMP: near-miss difference patterns per blocked row = {[len(nm_uniq.get(j, ())) for j in range(4)]}")
assert sorted(nm_uniq.keys()) == [0, 1, 2, 3] and all(len(v) >= 1 for v in nm_uniq.values())
# every pattern blocks at least one row
assert all(any(w < 1 for w in r[2]) for r in rows)
ok("CGLMP: every pattern blocks at least one row (cap lower-bound ingredient)")

print(f"\nALL CHECKS PASSED ({NCHECK} checks).  elapsed {time.time()-T0:.1f}s")
sys.exit(0)
