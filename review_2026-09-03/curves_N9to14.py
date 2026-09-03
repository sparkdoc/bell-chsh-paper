#!/usr/bin/env python3
"""Independent exhaustive check that S_max(N,F) equals the round-robin closed form
min{2 + sum_ab 2 min(F, n_ab(Q_N)/N), 4} for N = 9..14, at all grid points k/N and all cell midpoints
(convexity on each cell then gives the whole curve).  Exact Fraction arithmetic; no project code."""
import itertools, sys, time
from fractions import Fraction as Fr
SIG = (1, 1, 1, -1)
ALL = sorted({(a*c, a*d, b*c, b*d) for a, b, c, d in itertools.product((1, -1), repeat=4)})
CV = [q[0]+q[1]+q[2]-q[3] for q in ALL]
WRONG = [[int(q[ab] != SIG[ab]) for ab in range(4)] for q in ALL]
ok_all = True
for N in range(9, 15):
    t0 = time.time()
    grid = [Fr(k, 2*N) for k in range(2*N+1)]
    best = [None]*len(grid)
    for Q in itertools.combinations_with_replacement(range(8), N):
        S0 = Fr(sum(CV[i] for i in Q), N)
        n = [sum(WRONG[i][ab] for i in Q) for ab in range(4)]
        act = [Fr(c, N) for c in n if 0 < c < N]
        for j, F in enumerate(grid):
            v = S0 + sum(2*min(F, c) for c in act)
            if best[j] is None or v > best[j]: best[j] = v
    ns = [N//4 + (1 if j < N % 4 else 0) for j in range(4)]
    rr = [min(2 + sum(2*min(F, Fr(c, N)) for c in ns), Fr(4)) for F in grid]
    ok = best == rr; ok_all &= ok
    print(f"N={N}: {'PASS' if ok else 'FAIL'}  exhaustive envelope == round-robin closed form at {len(grid)} points "
          f"(cap {Fr(-(-N//4), N)}); {time.time()-t0:.0f}s", flush=True)
print("SUMMARY:", "ALL PASS" if ok_all else "FAIL"); sys.exit(0 if ok_all else 1)
