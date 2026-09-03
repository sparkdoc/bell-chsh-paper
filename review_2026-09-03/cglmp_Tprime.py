#!/usr/bin/env python3
"""F*(T';N) for CGLMP d=3 with the ADGL target T' = 1+sqrt(11/3), N=2,3,4, by exhaustive floating-point
enumeration over all multisets of the 27 difference patterns (per-multiset exact piecewise-linear curve
via the generalized row lemma evaluated at grid points k/N).  Floating point: a cross-check, not an exact
certificate."""
import itertools, math
from fractions import Fraction as Fr
def w00(d): return (d == 0) - (d == 2)
def w01(d): return (d == 0) - (d == 1)
def w10(d): return (d == 2) - (d == 0)
def w11(d): return (d == 0) - (d == 2)
pats = sorted({((A1-B1) % 3, (A1-B2) % 3, (A2-B1) % 3, (A2-B2) % 3) for A1, A2, B1, B2 in itertools.product(range(3), repeat=4)})
W = [tuple(f(d[i]) for i, f in enumerate((w00, w01, w10, w11))) for d in pats]
vals = [sum(w) for w in W]
def h_gain(w, N, k):
    ws = sorted(w)
    if k >= N: return Fr(max(w)) - Fr(sum(w), N)
    return Fr(k*max(w) - sum(ws[:k]), N)
def Fstar(T, N):
    best = 9
    for Q in itertools.combinations_with_replacement(range(27), N):
        V0 = sum(vals[i] for i in Q)/N
        ws = [[W[i][r] for i in Q] for r in range(4)]
        g = [[float(h_gain(ws[r], N, k)) for k in range(N+1)] for r in range(4)]
        def V(F):
            k = min(int(F*N), N-1); fr = F*N - k
            return V0 + sum(g[r][k] + fr*(g[r][k+1]-g[r][k]) for r in range(4))
        if V(1.0) < T: continue
        lo, hi = 0.0, 1.0
        for _ in range(60):
            mid = (lo+hi)/2
            if V(mid) >= T: hi = mid
            else: lo = mid
        best = min(best, hi)
    return best
for name, T in (("T_ME=(12+8sqrt3)/9", (12+8*math.sqrt(3))/9), ("T_Q=1+sqrt(11/3)", 1+math.sqrt(11/3))):
    for N in (2, 3, 4):
        got = Fstar(T, N); exp = (T-2)/(2*min(N, 4))
        print(f"{name}  N={N}: F* = {got:.9f}   (T-2)/(2 min(N,4)) = {exp:.9f}   {'MATCH' if abs(got-exp) < 1e-9 else 'DIFFER'}")
