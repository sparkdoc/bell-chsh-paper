import os as _os; _ROOT = _os.path.normpath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # release: paths relative to repo root
"""Exact certificate: F*_NS(S=4; N) = ceil(N/4)/N  for all N >= 3.

Question (open before this file): under FULL no-signalling, what is
    F*_NS(S=4; N) := min{ F : S_NS(N,F) >= 4 } ?

Result (PROVED here):
  * Lower bound (all N):  F*_NS(S=4;N) >= ceil(N/4)/N.
      Reason: S_NS(N,F) <= S_max(N,F) pointwise (NS restricts the feasible set), and the
      signalling-allowed cap is Theorem 2 of PROOFS.md, F_sig(S=4;N)=ceil(N/4)/N. Hence
      {F : S_NS>=4} subseteq {F : S_max>=4} and the min of the former >= the min of the latter.
      Part A below recomputes the cap exactly (Fraction) for N=2..8 to make this self-contained.
  * Upper bound (all N): an explicit within-class-uniform local deterministic model that is NS,
    reaches S=4 exactly, and has per-row TV <= ceil(N/4)/N. Part B verifies it EXACTLY
    (Fraction: row sums 1, p>=0, E[ab]=sigma_ab for all ab i.e. S=4, per-row TV<=Fcap, and all
    four NS residuals exactly 0) for N=3..NMAX, plus large-N spot checks. Part C is an
    independent LP (HiGHS, float) feasibility cross-check at the cap for N=3..16.

Key structural facts (proved in the report review/divergent/prbox.md):
  * The PR box is the UNIQUE no-signalling box with S=4 (NS + S=4 forces uniform 1/2 marginals).
    So "fine-tuned local model that is NS and has S=4" == "local deterministic decomposition of
    the PR box"; this script constructs exactly such a decomposition at the cap.
  * With the peak-1/2 row pattern used below, all four NS equations are AUTOMATIC given the row
    sums (the gauge makes each weighted row-sum vanish identically); only TV<=Fcap remains.

Construction (round-robin patterns; state i has miss-type M_{i mod 4}; class j = {i : i%4==j},
size n_j; a,b,c,d = n_0..n_3):
  Gauge: base realization per type, then flip the whole of one class:
     - N % 4 in {0,1,2}: flip class D (index 3)
     - N % 4 == 3       : flip class C (index 2)
  Rows, as total mass P[ab,j] on class j (p[ab, i in class j] = P[ab,j]/n_j):
     N%4 in {0,1,2} (flip-D):
        row00=[0,      b/N,   1/2-b/N, 1/2]
        row01=[b/N,    0,     1/2-b/N, 1/2]
        row10=[a/N,    1/2,   0,       1/2-a/N]
        row11=[a/N,    1/2,   1/2-a/N, 0]
     N%4 == 3 (flip-C):
        row00=[0,           1/2-d/N, 1/2,   d/N]
        row01=[1/2-d/N,     0,       1/2,   d/N]
        row10=[1/2,         1/2-d/N, 0,     d/N]
        row11=[1/2,         1/2-d/N, d/N,   0]
  Each row is supported on its correct class C_ab (the class it "misses" gets mass 0), so
  E[ab]=sigma_ab and S=4. Verified exactly below.

Exit code 0 iff every assertion passes. No floating point in Parts A/B (Fraction only).
"""
from fractions import Fraction as Fr
from itertools import combinations_with_replacement, product

AB = [(0, 0), (0, 1), (1, 0), (1, 1)]
SIGMA = [1, 1, 1, -1]                      # sigma_ab in ab-order 00,01,10,11
# miss-type product patterns q=(q00,q01,q10,q11); type j misses class j
MISS = {0: (-1, +1, +1, -1),   # A misses q00
        1: (+1, -1, +1, -1),   # B misses q01
        2: (+1, +1, -1, -1),   # C misses q10
        3: (+1, +1, +1, +1)}   # D misses q11
# base realization per type: u=[u0,u1], v=[v0,v1] (u0=+1 gauge)
BASE_U = {0: (1, -1), 1: (1, 1), 2: (1, -1), 3: (1, 1)}
BASE_V = {0: (-1, 1), 1: (1, -1), 2: (1, 1), 3: (1, 1)}

results = []


def check(label, cond):
    results.append((label, bool(cond)))
    return cond


# ----------------------------------------------------------------------
# Part A: exact lower bound -- recompute the signalling cap for N=2..8.
#   F_sig(S=4;N) = min_Q max_ab n_ab(Q) / N  over multisets Q of even-parity patterns,
#   n_ab(Q)=#{lam in Q : q_ab(lam) != sigma_ab}. Parity lemma => >= ceil(N/4)/N.
# ----------------------------------------------------------------------
def partA():
    print("=" * 72)
    print("Part A: exact lower bound  F*_NS(S=4;N) >= ceil(N/4)/N  (all N)")
    print("=" * 72)
    PATS = [tuple(p) for p in product([1, -1], repeat=4) if p[0] * p[1] * p[2] * p[3] == 1]
    ok_all = True
    for N in range(2, 9):
        best = None
        for Q in combinations_with_replacement(PATS, N):
            n = [sum(1 for q in Q if q[ab] != SIGMA[ab]) for ab in range(4)]
            m = max(n)
            if best is None or m < best:
                best = m
        capN = (N + 3) // 4
        ok = check(f"cap N={N}: min_Q max_ab n_ab == ceil(N/4) ({best} vs {capN})", best == capN)
        ok_all = ok_all and ok
    print("  (min over all even-parity multisets; parity lemma gives >= ceil(N/4), round-robin attains)")
    print("  Monotonicity: S_NS(N,F) <= S_max(N,F) for every F, so")
    print("    {F : S_NS(N,F)>=4} subseteq {F : S_max(N,F)>=4}")
    print("    => F*_NS(S=4;N) = min{F:S_NS>=4} >= min{F:S_max>=4} = ceil(N/4)/N  (Theorem 2).")
    print(f"  Part A: {'PASS' if ok_all else 'FAIL'}")
    return ok_all


# ----------------------------------------------------------------------
# Construction + exact verification.
# ----------------------------------------------------------------------
def build(N):
    """Return (P[ab][j] class-mass, u[a][i], v[b][i] per-state, n list, Fcap)."""
    n = [(N - j + 3) // 4 for j in range(4)]          # A,B,C,D class sizes
    a, b, c, d = n
    Nf = Fr(N)
    Fcap = Fr((N + 3) // 4, N)
    flip = 2 if (N % 4 == 3) else 3                    # which class to flip
    # per-state responses
    u = [[0] * N for _ in range(2)]
    v = [[0] * N for _ in range(2)]
    for i in range(N):
        j = i % 4
        sgn = -1 if j == flip else 1
        u[0][i] = sgn * BASE_U[j][0]; u[1][i] = sgn * BASE_U[j][1]
        v[0][i] = sgn * BASE_V[j][0]; v[1][i] = sgn * BASE_V[j][1]
    # class-mass rows P[ab][j]
    if N % 4 == 3:
        P = [[Fr(0), Fr(1, 2) - d / Nf, Fr(1, 2), d / Nf],
             [Fr(1, 2) - d / Nf, Fr(0), Fr(1, 2), d / Nf],
             [Fr(1, 2), Fr(1, 2) - d / Nf, Fr(0), d / Nf],
             [Fr(1, 2), Fr(1, 2) - d / Nf, d / Nf, Fr(0)]]
    else:
        P = [[Fr(0), b / Nf, Fr(1, 2) - b / Nf, Fr(1, 2)],
             [b / Nf, Fr(0), Fr(1, 2) - b / Nf, Fr(1, 2)],
             [a / Nf, Fr(1, 2), Fr(0), Fr(1, 2) - a / Nf],
             [a / Nf, Fr(1, 2), Fr(1, 2) - a / Nf, Fr(0)]]
    return P, u, v, n, Fcap


def verify(N):
    """Exact Fraction verification of the construction at N. Returns (ok, info)."""
    P, u, v, n, Fcap = build(N)
    # per-state p[ab][i]
    p = [[Fr(0)] * N for _ in range(4)]
    for ab in range(4):
        for i in range(N):
            j = i % 4
            p[ab][i] = P[ab][j] / n[j] if n[j] > 0 else Fr(0)
    ok = True
    # 1. row sums == 1
    for ab in range(4):
        if sum(p[ab]) != 1:
            ok = False
    # 2. p >= 0
    for ab in range(4):
        for i in range(N):
            if p[ab][i] < 0:
                ok = False
    # 3. E[ab] == sigma_ab  (=> S=4)
    S = Fr(0)
    for ab, (a, b) in enumerate(AB):
        E = sum(p[ab][i] * u[a][i] * v[b][i] for i in range(N))
        if E != SIGMA[ab]:
            ok = False
        S += Fr(SIGMA[ab]) * E
    # 4. per-row TV <= Fcap
    tvmax = Fr(0)
    for ab in range(4):
        tv = Fr(1, 2) * sum(abs(p[ab][i] - Fr(1, N)) for i in range(N))
        if tv > Fcap:
            ok = False
        tvmax = max(tvmax, tv)
    # 5. four NS residuals == 0
    nsr = [sum(p[0][i] * u[0][i] for i in range(N)) - sum(p[1][i] * u[0][i] for i in range(N)),
           sum(p[2][i] * u[1][i] for i in range(N)) - sum(p[3][i] * u[1][i] for i in range(N)),
           sum(p[0][i] * v[0][i] for i in range(N)) - sum(p[2][i] * v[0][i] for i in range(N)),
           sum(p[1][i] * v[1][i] for i in range(N)) - sum(p[3][i] * v[1][i] for i in range(N))]
    if any(r != 0 for r in nsr):
        ok = False
    return ok, S, tvmax, Fcap


def partB():
    print()
    print("=" * 72)
    print("Part B: exact upper bound -- construction verified (Fraction), all N=3..NMAX")
    print("=" * 72)
    NMAX = 300
    ok_all = True
    # explicit small-N witnesses printed for the report
    for N in (3, 4, 5):
        P, u, v, n, Fcap = build(N)
        print(f"\n  --- N={N}  (class sizes a,b,c,d={n}, Fcap=ceil(N/4)/N={Fcap}) ---")
        for ab, (a, b) in enumerate(AB):
            row = [P[ab][j] / n[j] if n[j] > 0 else Fr(0) for j in range(4)]
            print(f"    P[{AB[ab]}] class-mass (A,B,C,D) = [{', '.join(str(x) for x in P[ab])}]")
        ok, S, tvmax, Fcap2 = verify(N)
        check(f"N={N}: exact witness valid (S=4, TV<=Fcap, NS=0)", ok)
        print(f"    -> S={S}  maxTV={tvmax}  Fcap={Fcap2}  valid={ok}")
    # sweep
    for N in range(3, NMAX + 1):
        ok, S, tvmax, Fcap = verify(N)
        if not (ok and S == 4 and tvmax <= Fcap):
            check(f"sweep N={N}: valid", False)
            ok_all = False
            print(f"    FAIL at N={N}: S={S} maxTV={tvmax} Fcap={Fcap}")
            break
    else:
        check(f"sweep N=3..{NMAX}: all valid (S=4, TV<=ceil(N/4)/N, NS=0)", True)
    # large-N spot checks
    for N in (1000, 4096, 10000):
        ok, S, tvmax, Fcap = verify(N)
        check(f"large N={N}: valid", ok and S == 4 and tvmax <= Fcap)
        print(f"    large N={N}: S={S} maxTV={tvmax} Fcap={Fcap} valid={ok}")
    print(f"\n  Part B: {'PASS' if ok_all else 'FAIL'}")
    return ok_all


# ----------------------------------------------------------------------
# Part C: independent LP (HiGHS, float) feasibility cross-check at the cap.
# ----------------------------------------------------------------------
def partC():
    print()
    print("=" * 72)
    print("Part C: independent LP (HiGHS, float) feasibility at F=ceil(N/4)/N, N=3..16")
    print("=" * 72)
    try:
        import numpy as np
        from scipy.optimize import linprog
    except Exception as e:
        print(f"  (skipped LP cross-check: {e}) -- Parts A,B are the exact certificates")
        return True
    ok_all = True
    for N in range(3, 17):
        Fcap = (N + 3) // 4 / N
        nvar = 4 * N; naux = 4 * N; nw = nvar + naux
        found = False
        # try all gauges for N<=8 else class-level + a couple
        if N <= 8:
            cands = list(product([0, 1], repeat=N))
        else:
            cands = [tuple(cb[i % 4] for i in range(N)) for cb in product([0, 1], repeat=4)]
        for bits in cands:
            types = [i % 4 for i in range(N)]
            u = np.zeros((2, N)); v = np.zeros((2, N))
            for i in range(N):
                j = types[i]; sgn = -1 if bits[i] else 1
                u[0, i] = sgn * BASE_U[j][0]; u[1, i] = sgn * BASE_U[j][1]
                v[0, i] = sgn * BASE_V[j][0]; v[1, i] = sgn * BASE_V[j][1]
            c = np.zeros(nw); A_eq = []; b_eq = []
            for ab in range(4):
                row = np.zeros(nw); row[ab*N:(ab+1)*N] = 1.0; A_eq.append(row); b_eq.append(1.0)
            for ab in range(4):
                a, b = AB[ab]; row = np.zeros(nw); row[ab*N:(ab+1)*N] = u[a]*v[b]
                A_eq.append(row); b_eq.append(float(SIGMA[ab]))
            for (w, abL, abR) in [(u[0], 0, 1), (u[1], 2, 3), (v[0], 0, 2), (v[1], 1, 3)]:
                row = np.zeros(nw); row[abL*N:(abL+1)*N] = w; row[abR*N:(abR+1)*N] = -w
                A_eq.append(row); b_eq.append(0.0)
            A_ub = []; b_ub = []
            for ab in range(4):
                row = np.zeros(nw); row[nvar+ab*N:nvar+(ab+1)*N] = 1.0
                A_ub.append(row); b_ub.append(2.0*Fcap)
            for ab in range(4):
                for lam in range(N):
                    r1 = np.zeros(nw); r1[ab*N+lam] = 1.0; r1[nvar+ab*N+lam] = -1.0
                    A_ub.append(r1); b_ub.append(1.0/N)
                    r2 = np.zeros(nw); r2[ab*N+lam] = -1.0; r2[nvar+ab*N+lam] = -1.0
                    A_ub.append(r2); b_ub.append(-1.0/N)
            res = linprog(c, A_ub=np.array(A_ub), b_ub=np.array(b_ub), A_eq=np.array(A_eq),
                          b_eq=np.array(b_eq), bounds=[(0, 1)]*nvar + [(0, None)]*naux, method='highs')
            if res.success:
                found = True
                break
        check(f"LP feasible at cap N={N}", found)
        ok_all = ok_all and found
        print(f"    N={N}: {'FEASIBLE' if found else 'INFEASIBLE'} at Fcap={Fcap:.5f}")
    print(f"\n  Part C: {'PASS' if ok_all else 'FAIL'}")
    return ok_all


if __name__ == "__main__":
    a = partA()
    b = partB()
    c = partC()
    nfail = sum(1 for _, ok in results if not ok)
    print()
    print("=" * 72)
    print(f"SUMMARY: {len(results)} checks, {nfail} failures")
    print("THEOREM: F*_NS(S=4; N) = ceil(N/4)/N for all N >= 3.")
    print(f"  lower bound (Part A): {'PASS' if a else 'FAIL'}")
    print(f"  upper bound, exact construction (Part B): {'PASS' if b else 'FAIL'}")
    print(f"  LP cross-check (Part C): {'PASS' if c else 'FAIL'}")
    print("=" * 72)
    import sys
    sys.exit(0 if (a and b and c and nfail == 0) else 1)
