#!/usr/bin/env python3
"""Correspondence between Puetz et al. (PRL 113, 190402, 2014) Eq. (11), CHSH_max = 4(1-2l) for MDL models
with no-signalling and uniform inputs (l = min_{xy,lambda} P(xy|lambda)), and the main paper's S <= 2+8F with
F = max_xy TV(p(.|xy) || uniform).  Exact arithmetic in Q(sqrt2) via Fractions + symbolic sqrt2 where needed.

Checks:
 A. Algebra: 4(1-2l) with l = 1/4 - F equals 2+8F; Puetz threshold l0=(2-sqrt2)/4 equals 1/4 - (sqrt2-1)/4.
 B. Puetz-type extremal model = round-robin over the four CHSH miss-strategies with P(xy|lambda) = l on the
    strategy's one wrong input pair and (1-l)/3 on the other three, rho(lambda)=1/4:  uniform inputs, S=4(1-2l),
    and in the TV metric every row has TV = 1/4 - l exactly  -> the two extremal families coincide at N=4.
 C. The converse map is not an equivalence: for the paper's own N=4 NS witness, l = (1/4-d)/(1+d) != 1/4-d.
 D. A no-signalling response gauge exists for the Puetz-type model (so it is also an NS witness in our sense).
"""
from fractions import Fraction as Fr
import itertools, math

SIG = (1, 1, 1, -1)                       # sigma_ab for ab = 00,01,10,11
MISS = [(-1, 1, 1, -1), (1, -1, 1, -1), (1, 1, -1, -1), (1, 1, 1, 1)]   # M_A..M_D, each wrong on exactly one ab
wrong = [next(ab for ab in range(4) if MISS[j][ab] != SIG[ab]) for j in range(4)]   # = [0,1,2,3]

# ---- A (float sanity + exact statement)
for F in (Fr(0), Fr(1, 10), Fr(1, 4)):
    l = Fr(1, 4) - F
    assert 4*(1-2*l) == 2 + 8*F
print("A: 4(1-2l) == 2+8F under l = 1/4 - F   (exact, identity in F)")
print(f"   Puetz threshold l0 = (2-sqrt2)/4 = {(2-math.sqrt(2))/4:.6f};  1/4 - F* = 1/4 - (sqrt2-1)/4 = {0.25-(math.sqrt(2)-1)/4:.6f}")

# ---- B: Puetz-type model in the paper's variables, exact for rational l
def puetz_model(l):
    """rows p[ab, lambda] = P(lambda|ab) = P(ab|lambda) rho(lambda) / P(ab) = P(ab|lambda) (uniform inputs, rho=1/4)"""
    P = [[None]*4 for _ in range(4)]
    for lam in range(4):
        for ab in range(4):
            P[ab][lam] = l if wrong[lam] == ab else (1-l)/3
    return P
for l in (Fr(1, 4), Fr(1, 5), Fr(1, 8), Fr(0)):
    P = puetz_model(l)
    assert all(sum(row) == 1 for row in P)                                   # rows are distributions
    inputs = [sum(P[ab][lam] for ab in range(4))/4 for lam in range(4)]       # P(lambda) = sum_ab P(ab) P(lambda|ab)
    assert all(x == Fr(1, 4) for x in inputs)                                 # marginal = uniform baseline: uniform inputs
    E = [sum(P[ab][lam]*MISS[lam][ab] for lam in range(4)) for ab in range(4)]
    S = sum(SIG[ab]*E[ab] for ab in range(4))
    TV = [sum(abs(P[ab][lam]-Fr(1, 4)) for lam in range(4))/2 for ab in range(4)]
    assert S == 4*(1-2*l) and all(t == Fr(1, 4)-l for t in TV)
    print(f"B: l={str(l):4}  S = {S} = 4(1-2l);  per-row TV = {TV[0]} = 1/4 - l;  min P(ab|lam) = {min(min(r) for r in P)} = l")

# ---- C: the paper's N=4 NS witness (Thm 7 construction) in Puetz's variable, symbolic in d
# rows: p[ab, lam] = 1/4 + x[ab,lam]*d with x from Eq. (ns4rows): class counts n_j = 1 at N=4
X = {(0, 0): -1, (0, 1): +1, (1, 0): +1, (1, 1): -1, (2, 1): +1, (2, 2): -1, (3, 0): +1, (3, 3): -1}
def ns_witness_rows(d):   # d as Fraction (numeric probe)
    return [[Fr(1, 4) + d*X.get((ab, lam), 0) for lam in range(4)] for ab in range(4)]
d = Fr(1, 10)
P = ns_witness_rows(d)
rho_true = [sum(P[ab][lam] for ab in range(4))/4 for lam in range(4)]           # actual marginal (uniform inputs)
ratios = [P[ab][lam]*Fr(1, 4)/rho_true[lam] for ab in range(4) for lam in range(4)]   # P(ab|lam) = P(lam|ab)P(ab)/P(lam)
print(f"C: paper's NS witness at d={d}: true marginal rho = {[str(r) for r in rho_true]} (NOT uniform);")
print(f"   l = min P(ab|lam) = {min(ratios)} = (1/4-d)/(1+d) = {(Fr(1,4)-d)/(1+d)};  1/4 - d = {Fr(1,4)-d}  -> different point in Puetz's variable")

# ---- D: NS gauge for the Puetz-type model (search all 2^16 response tables realising the four miss patterns)
def realises(u0, u1, v0, v1):
    return all((u0[lam]*v0[lam], u0[lam]*v1[lam], u1[lam]*v0[lam], u1[lam]*v1[lam]) == MISS[lam] for lam in range(4))
def ns_ok(P, u0, u1, v0, v1):
    U = (u0, u1); V = (v0, v1)
    for a in range(2):
        m = [sum(P[2*a+b][lam]*U[a][lam] for lam in range(4)) for b in range(2)]
        if m[0] != m[1]: return False
    for b in range(2):
        m = [sum(P[2*a+b][lam]*V[b][lam] for lam in range(4)) for a in range(2)]
        if m[0] != m[1]: return False
    return True
P = puetz_model(Fr(1, 8)); found = []
for bits in itertools.product((1, -1), repeat=16):
    u0, u1, v0, v1 = bits[:4], bits[4:8], bits[8:12], bits[12:]
    if realises(u0, u1, v0, v1) and ns_ok(P, u0, u1, v0, v1): found.append((u0, u1, v0, v1))
print(f"D: response gauges realising M_A..M_D that make the Puetz-type model no-signalling: {len(found)} of 16"
      + (f"; e.g. u0={found[0][0]} u1={found[0][1]} v0={found[0][2]} v1={found[0][3]}" if found else ""))
print("\nCONCLUSION: Puetz et al. 2014 Eq. (11) and the paper's N=4 headline describe the same one-parameter family of "
      "extremal models (round-robin over the four miss strategies), related by l = 1/4 - F; the linear law and the "
      "Tsirelson threshold coincide.  The metrics are not equivalent in general (part C), so the TV bound for arbitrary "
      "models and N is formally a different statement, but the headline is a re-parametrisation of a 2014 result.")
