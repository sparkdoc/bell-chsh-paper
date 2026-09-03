#!/usr/bin/env python3
"""
Independent re-verification of the companion paper (paper/companion.tex), written from scratch
2026-09-03.  No project code reused.  Exit 0 iff all assertions pass; NOTE: lines are observations.
Run: python3 indep_check_companion.py
"""
import itertools, math, sys, time
from fractions import Fraction as Fr
import numpy as np

FAILS = []
def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")
    if not cond: FAILS.append(name)
def ceil_div(a, b): return -(-a//b)

# =====================================================================  Mermin n=3 and n=5
def mermin_patterns(n):
    """realizable correlator vectors over even-weight subsets S (Y-settings), value = sum (-1)^{|S|/2} e_S"""
    subsets = [S for k in range(0, n+1, 2) for S in itertools.combinations(range(n), k)]
    pats = set()
    for x in itertools.product((1, -1), repeat=n):
        for y in itertools.product((1, -1), repeat=n):
            e = tuple(math.prod(y[i] if i in S else x[i] for i in range(n)) for S in subsets)
            pats.add(e)
    target = tuple((-1)**(len(S)//2) for S in subsets)
    return subsets, sorted(pats), target
for n, K, npat, B, mmin, nmin, inc in ((3, 4, 8, 2, 1, 4, 1), (5, 16, 32, 4, 6, 16, 6)):
    subsets, pats, target = mermin_patterns(n)
    check(f"Mermin n={n}: K={K} even-weight rows", len(subsets) == K)
    check(f"Mermin n={n}: {npat} realizable patterns", len(pats) == npat, len(pats))
    vals = [sum(t*e for t, e in zip(target, p)) for p in pats]
    check(f"Mermin n={n}: local bound B={B}", max(vals) == B, max(vals))
    miss = [sum(1 for t, e in zip(target, p) if t != e) for p in pats]
    check(f"Mermin n={n}: m_min={mmin}", min(miss) == mmin, min(miss))
    minmiss = [p for p in pats if sum(1 for t, e in zip(target, p) if t != e) == mmin]
    check(f"Mermin n={n}: {nmin} min-miss patterns", len(minmiss) == nmin, len(minmiss))
    incid = [sum(1 for p in minmiss if p[s] != target[s]) for s in range(K)]
    check(f"Mermin n={n}: incidence {inc} per condition", set(incid) == {inc}, incid)
    if n == 5:
        # A_maxB: max #active rows among multisets of bound-achieving patterns, N=2..6
        masks = [tuple(int(p[s] != target[s]) for s in range(K)) for p in minmiss]
        for N in range(2, 7):
            best = 0
            for Q in itertools.combinations_with_replacement(range(len(masks)), N):
                cnt = [sum(masks[i][s] for i in Q) for s in range(K)]
                best = max(best, sum(1 for c in cnt if 0 < c < N))
            check(f"Mermin5: A_maxB(N={N}) = {(8,12,16,16,16)[N-2]}", best == (8, 12, 16, 16, 16)[N-2], best)
        # exact curve V(k/N) and F*(16;N) over ALL 32 patterns, N=2..6
        allmask = [tuple(int(p[s] != target[s]) for s in range(K)) for p in pats]
        allval = [sum(t*e for t, e in zip(target, p)) for p in pats]
        table = {2: ([4, 16], Fr(1, 2)), 3: ([4, 12, 16], Fr(2, 3)), 4: ([4, 12, 16], Fr(1, 2)),
                 5: ([4, Fr(52, 5), Fr(72, 5), 16], Fr(3, 5)), 6: ([4, Fr(28, 3), 14, 16], Fr(1, 2))}
        for N in range(2, 7):
            t0 = time.time()
            ks = list(range(0, len(table[N][0])))
            bestV = [None]*len(ks); bestcap = N
            for Q in itertools.combinations_with_replacement(range(32), N):
                cnt = [sum(allmask[i][s] for i in Q) for s in range(K)]
                V0 = Fr(sum(allval[i] for i in Q), N)
                bestcap = min(bestcap, max(cnt))
                for j, k in enumerate(ks):
                    F = Fr(k, N)
                    v = V0 + sum(2*min(F, Fr(c, N)) for c in cnt if 0 < c < N)
                    if bestV[j] is None or v > bestV[j]: bestV[j] = v
            check(f"Mermin5 N={N}: V(k/N) table", bestV == [Fr(x) for x in table[N][0]], f"{bestV} ({time.time()-t0:.0f}s)")
            check(f"Mermin5 N={N}: F*(16;N) = min_Q max_s n_s/N = {table[N][1]}", Fr(bestcap, N) == table[N][1], Fr(bestcap, N))
        # N=2: envelope bends at F=1/4 (V(1/4)=8, not 10): grid values alone do not give the curve
        F = Fr(1, 4); best = None
        for Q in itertools.combinations_with_replacement(range(32), 2):
            cnt = [sum(allmask[i][s] for i in Q) for s in range(K)]
            v = Fr(sum(allval[i] for i in Q), 2) + sum(2*min(F, Fr(c, 2)) for c in cnt if 0 < c < 2)
            best = v if best is None else max(best, v)
        print(f"NOTE: Mermin5 N=2: V(1/4) = {best} vs linear interpolation of grid values (4,16) = 10  -> the envelope "
              f"is NOT determined by V(k/N) alone (contradicts Sec. II.C sentence 'listing V(k/N) gives the complete curve')")
        # remainder-witness combination: equality also for r=3,4,6 (given exact N=3,4,6 values)
        print("NOTE: since N=3,4,6 attain ceil(3N/8)/N, the construction 'm copies of all 16 + Q_r' gives equality "
              "F* = ceil(3N/8)/N also for N = 3,4,6 (mod 16) -- Theorem states equality only for r=0,1,2 (weaker than provable).")

# GHZ quantum values (exact via complex integer arithmetic scaled)
X = np.array([[0, 1], [1, 0]], complex); Y = np.array([[0, -1j], [1j, 0]]); I2 = np.eye(2)
def kron(*ms):
    out = np.array([[1]], complex)
    for m in ms: out = np.kron(out, m)
    return out
for n in (3, 5):
    ghz = np.zeros(2**n, complex); ghz[0] = ghz[-1] = 1/math.sqrt(2)
    subsets = [S for k in range(0, n+1, 2) for S in itertools.combinations(range(n), k)]
    M = 0
    for S in subsets:
        op = kron(*[Y if i in S else X for i in range(n)])
        e = (ghz.conj() @ op @ ghz).real
        M += (-1)**(len(S)//2) * e
        assert abs(e - (-1)**(len(S)//2)) < 1e-12
    check(f"GHZ_{n}: Mermin value = {2**(n-1)} (algebraic max), all correlators perfect", abs(M - 2**(n-1)) < 1e-12, M)

# =====================================================================  CGLMP d=3
def w00(d): return (d == 0) - (d == 2)
def w01(d): return (d == 0) - (d == 1)
def w10(d): return (d == 2) - (d == 0)
def w11(d): return (d == 0) - (d == 2)
pats = set()
for A1, A2, B1, B2 in itertools.product(range(3), repeat=4):
    d = ((A1-B1) % 3, (A1-B2) % 3, (A2-B1) % 3, (A2-B2) % 3)
    pats.add(d)
pats = sorted(pats)
check("CGLMP: 27 realizable difference patterns", len(pats) == 27)
check("CGLMP: modular constraint d00+d11 = d01+d10 (mod 3)", all((d[0]+d[3]-d[1]-d[2]) % 3 == 0 for d in pats))
W = [tuple(f(d[i]) for i, f in enumerate((w00, w01, w10, w11))) for d in pats]
vals = [sum(w) for w in W]
check("CGLMP: local bound 2", max(vals) == 2)
target = (0, 0, 2, 0)   # algebraic-maximum conditions d00=0, d01=0, d10=2, d11=0
misses = [sum(1 for i in range(4) if d[i] != target[i]) for d in pats]
near = [pats[i] for i in range(27) if misses[i] == 1]
check("CGLMP: no pattern satisfies all four conditions (m_min=1)", min(misses) == 1)
check("CGLMP: four near-miss (exactly-one-miss) patterns as listed", set(near) == {(2, 0, 2, 0), (0, 1, 2, 0), (0, 0, 0, 0), (0, 0, 2, 2)}, near)
check("CGLMP: near-miss = three rows +1, one row -1", all(sorted(W[pats.index(p)]) == [-1, 1, 1, 1] for p in near))
bound_attaining = [pats[i] for i in range(27) if vals[i] == 2]
print(f"NOTE: CGLMP: {len(bound_attaining)} patterns attain the local bound 2 (4 near-miss of type (+1,+1,+1,-1) and "
      f"{len(bound_attaining)-4} of type (+1,+1,0,0)); the companion text mentions only the four near-miss patterns.")
# direct check of I_3 for the 81 deterministic strategies (independent of the difference reduction)
def I3_det(A1, A2, B1, B2):
    P = lambda cond: 1.0 if cond else 0.0
    return (P(A1 == B1) - P(A1 == (B1-1) % 3) + P(B1 == (A2+1) % 3) - P(B1 == A2)
            + P(A2 == B2) - P(A2 == (B2-1) % 3) + P(B2 == A1) - P(B2 == (A1-1) % 3))
check("CGLMP: max over 81 strategies = 2", max(I3_det(*s) for s in itertools.product(range(3), repeat=4)) == 2)
# generalized row lemma h(k/N) for weights in {-1,0,1} vs LP
from scipy.optimize import linprog
def lp_row(w, N, F):
    c = np.concatenate([-np.array(w, float), np.zeros(N)]); A, b = [], []
    for i in range(N):
        r = np.zeros(2*N); r[i] = 1; r[N+i] = -1; A.append(r); b.append(1/N)
        r = np.zeros(2*N); r[i] = -1; r[N+i] = -1; A.append(r); b.append(-1/N)
    r = np.zeros(2*N); r[N:] = 1; A.append(r); b.append(2*F)
    res = linprog(c, A_ub=np.array(A), b_ub=b, A_eq=np.array([np.concatenate([np.ones(N), np.zeros(N)])]), b_eq=[1],
                  bounds=[(0, None)]*(2*N), method="highs"); assert res.status == 0
    return -res.fun
def h_formula(w, N, k):
    ws = sorted(w); mu = Fr(sum(w), N)
    if k >= N: return Fr(max(w))
    return mu + Fr(k*max(w) - sum(ws[:k]), N)
dev = 0
for N in range(2, 7):
    for w in itertools.product((-1, 0, 1), repeat=N):
        for k in range(0, N+1):
            dev = max(dev, abs(lp_row(w, N, k/N) - float(h_formula(w, N, k))))
check("CGLMP generalized row lemma (Lemma 1 of companion) vs LP at F=k/N", dev < 1e-9, f"max dev {dev:.1e}")
# F*(T;N) for N=2,3 by float enumeration over all multisets (per-multiset piecewise linear crossing)
T = (12 + 8*math.sqrt(3))/9
def cglmp_Fstar(N):
    best = 9
    for Q in itertools.combinations_with_replacement(range(27), N):
        V0 = sum(vals[i] for i in Q)/N
        # curve = V0 + sum_r gain_r(F), gain_r piecewise linear from h_formula at grid; evaluate on fine grid
        ws = [[W[i][r] for i in Q] for r in range(4)]
        def V(F):
            k = int(F*N); frac = F*N - k
            tot = V0
            for r in range(4):
                lo = float(h_formula(ws[r], N, k) - Fr(sum(ws[r]), N))
                hi = float(h_formula(ws[r], N, min(k+1, N)) - Fr(sum(ws[r]), N))
                tot += lo + frac*(hi-lo)
            return tot
        lo, hi = 0.0, 1.0
        if V(1.0) < T: continue
        for _ in range(60):
            mid = (lo+hi)/2
            if V(mid) >= T: hi = mid
            else: lo = mid
        best = min(best, hi)
    return best
for N in (2, 3, 4):
    exp = (T-2)/(2*min(N, 4))
    got = cglmp_Fstar(N)
    check(f"CGLMP N={N}: F*(T;N) = (T-2)/(2 min(N,4)) = {exp:.6f}", abs(got-exp) < 1e-9, f"got {got:.9f}")
# quantum value: CGLMP maximally entangled state and Fourier measurements; and ADGL non-max entangled state
def cglmp_value(gamma):
    d = 3
    psi = np.zeros(9, complex)
    for j in range(3): psi[j*3+j] = [1, gamma, 1][j]
    psi /= np.linalg.norm(psi)
    om = np.exp(2j*math.pi/3)
    def meas(alpha, party):
        # Fourier-type: |k>_a = 1/sqrt3 sum_j exp(i 2pi j (k+alpha)/3)|j>  (Bob with -j convention)
        vecs = []
        for k in range(3):
            sgn = 1 if party == 0 else -1   # Bob: exp(i 2pi j (-k + beta)/3)  (CGLMP convention)
            v = np.array([np.exp(1j*2*math.pi*j*(sgn*k+alpha)/3) for j in range(3)])/math.sqrt(3)
            vecs.append(v)
        return vecs
    A = [meas(0.0, 0), meas(0.5, 0)]; B = [meas(0.25, 1), meas(-0.25, 1)]
    def P(a, b, ka, kb):
        proj = np.kron(np.outer(A[a][ka], A[a][ka].conj()), np.outer(B[b][kb], B[b][kb].conj()))
        return (psi.conj() @ proj @ psi).real
    def Peq(a, b, shift):  # P(A_a = B_b + shift mod 3)
        return sum(P(a, b, (kb+shift) % 3, kb) for kb in range(3))
    # I3 = [P(A1=B1) - P(A1=B1-1)] + [P(B1=A2+1) - P(B1=A2)] + [P(A2=B2) - P(A2=B2-1)] + [P(B2=A1) - P(B2=A1-1)]
    return (Peq(0, 0, 0) - Peq(0, 0, -1) + Peq(1, 0, -1) - Peq(1, 0, 0) + Peq(1, 1, 0) - Peq(1, 1, -1)
            + Peq(0, 1, 0) - Peq(0, 1, 1))
v1 = cglmp_value(1.0)
gam = (math.sqrt(11) - math.sqrt(3))/2
v2 = cglmp_value(gam)
print(f"NOTE: CGLMP I3 numerics: max-entangled state -> {v1:.6f} (paper's T={T:.6f}); "
      f"ADGL state gamma={gam:.4f} -> {v2:.6f} (1+sqrt(11/3) = {1+math.sqrt(11/3):.6f})")
check("CGLMP: maximally entangled Fourier value equals T=(12+8sqrt3)/9", abs(v1 - T) < 1e-9, f"{v1:.9f}")
if v2 > T + 1e-6:
    print(f"NOTE: a LARGER qutrit quantum value than the paper's target exists (ADGL 2002): I3 = {v2:.6f} > T; "
          f"the companion's headline F* for N>=4 would then be (I3-2)/8 = {(v2-2)/8:.6f}, not {(T-2)/8:.6f}")

# =====================================================================  Peres-Mermin square
sx = np.array([[0, 1], [1, 0]]); sy = np.array([[0, -1j], [1j, 0]]); sz = np.array([[1, 0], [0, -1]]); i2 = np.eye(2)
grid = [[np.kron(sx, sx), np.kron(sx, i2), np.kron(i2, sx)],
        [np.kron(sy, sy), np.kron(sx, sz), np.kron(sz, sx)],
        [np.kron(sz, sz), np.kron(i2, sz), np.kron(sz, i2)]]
lines = [[(r, c) for c in range(3)] for r in range(3)] + [[(r, c) for r in range(3)] for c in range(3)]
prods = []
for L in lines:
    M = np.eye(4, dtype=complex)
    for (r, c) in L: M = M @ grid[r][c]
    prods.append(M)
    comm = all(np.allclose(grid[a][b] @ grid[e][f], grid[e][f] @ grid[a][b]) for (a, b) in L for (e, f) in L)
    assert comm
signs = [1 if np.allclose(M, np.eye(4)) else (-1 if np.allclose(M, -np.eye(4)) else 0) for M in prods]
check("PM square: line products (R0,R1,R2,C0,C1,C2) = (+,+,+,-,+,+), all lines commuting", signs == [1, 1, 1, -1, 1, 1], signs)
t = signs
misssets = set()
for v in itertools.product((1, -1), repeat=9):
    s = tuple(t[c]*math.prod(v[3*r+cc] for (r, cc) in lines[c]) for c in range(6))
    misssets.add(tuple(int(x == -1) for x in s))
check("PM: every assignment misses an odd number of lines", all(sum(m) % 2 == 1 for m in misssets))
check("PM: exactly 32 odd miss-sets occur (6 singles, 20 triples, 6 quintuples)",
      len(misssets) == 32 and sorted(sum(m) for m in misssets).count(1) == 6)
check("PM: a near-miss type exists for every line", all(tuple(int(i == c) for i in range(6)) in misssets for c in range(6)))
# exact curves by DP over achievable w-vectors (N up to 10)
def pm_curve(N):
    layer = {(0,)*6}
    for _ in range(N):
        layer = {tuple(w[i]+m[i] for i in range(6)) for w in layer for m in misssets}
    def G(F):
        return max(Fr(1) - Fr(sum(w), 6*N) + Fr(1, 6)*sum(min(F, Fr(wc, N)) for wc in w) for w in layer)
    return G, layer
paper = {2: lambda F: Fr(5, 6)+Fr(1, 3)*F, 3: lambda F: Fr(5, 6)+Fr(1, 2)*F, 4: lambda F: Fr(5, 6)+Fr(2, 3)*F,
         5: lambda F: Fr(5, 6)+Fr(5, 6)*F, 6: lambda F: Fr(5, 6)+F,
         7: lambda F: min(Fr(5, 6)+F, Fr(41, 42)+Fr(1, 6)*(F-Fr(1, 7))),
         8: lambda F: min(Fr(5, 6)+F, Fr(23, 24)+Fr(1, 3)*(F-Fr(1, 8))),
         9: lambda F: min(Fr(5, 6)+F, Fr(17, 18)+Fr(1, 2)*(F-Fr(1, 9))),
         10: lambda F: min(Fr(5, 6)+F, Fr(14, 15)+Fr(2, 3)*(F-Fr(1, 10)))}
for N in range(2, 11):
    t0 = time.time()
    G, layer = pm_curve(N)
    cap = Fr(min(max(w) for w in layer), N)
    check(f"PM N={N}: cap F*_PM = ceil(N/6)/N", cap == Fr(ceil_div(N, 6), N), cap)
    Fs = Fr(ceil_div(N, 6), N)
    grid = sorted({Fr(k, 2*N) for k in range(0, 2*N+1) if Fr(k, 2*N) <= Fs})
    bad = [(F, G(F), min(paper[N](F), 1)) for F in grid if G(F) != min(paper[N](F), Fr(1))]
    check(f"PM N={N}: exact curve equals paper table on grid+midpoints up to F*", not bad, f"{bad[:2]} ({time.time()-t0:.0f}s, {len(layer)} w-vectors)")

print("\nSUMMARY:", "ALL PASS" if not FAILS else f"{len(FAILS)} FAIL(S): {FAILS}")
sys.exit(0 if not FAILS else 1)
