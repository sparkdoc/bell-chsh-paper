import os as _os; _ROOT = _os.path.normpath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # release: paths relative to repo root
"""NS general-N check via state-splitting embedding of the exact N=4 NS witness.

Take the saved N=4 no-signalling witness (ns_solN4_{p,u,v}.npy) achieving
S = 2*sqrt(2) at F* = (sqrt(2)-1)/4 with all four rows at TV = F*.  Split state
lambda into m_lambda copies (m in {floor(N/4), ceil(N/4)}, sum m = N), keep u,v
identical on copies, and set p'[ab,(lambda,j)] = p4[ab,lambda]/m_lambda.

Then:
  E'[ab]   = sum_lambda p4[ab,lambda] q_ab(lambda)          (unchanged -> S = 2*sqrt(2))
  NS sums  = sum_lambda p4[ab,lambda] u[a,lambda] etc.      (unchanged -> NS holds)
  TV'_ab   = (1/2) sum_lambda |p4[ab,lambda] - m_lambda/N|  (the only thing that changes)

So the embedding is a valid NS model at N with S = 2*sqrt(2) iff max_ab TV'_ab <= F*.
We try all 4 choices of which states receive the ceil copies (rem = N mod 4 of them).
"""
import numpy as np

Fstar = (np.sqrt(2.0) - 1.0) / 4.0
p4 = np.load(_ROOT + "/witnesses/ns_solN4_p.npy")     # (4,4): rows ab in order 00,01,10,11
u4 = np.load(_ROOT + "/witnesses/ns_solN4_u.npy")     # (2,4)
v4 = np.load(_ROOT + "/witnesses/ns_solN4_v.npy")     # (2,4)

# sanity: the witness itself
E4 = np.array([(p4[ab] * u4[a] * v4[b]).sum() for ab, (a, b) in enumerate([(0,0),(0,1),(1,0),(1,1)])])
S4 = E4[0] + E4[1] + E4[2] - E4[3]
tv4 = np.array([0.5 * np.abs(p4[ab] - 0.25).sum() for ab in range(4)])
print(f"sanity N=4 witness: S = {S4:.12f} (2*sqrt2 = {2*np.sqrt(2):.12f}), per-row TV = {np.round(tv4, 12)}")
assert abs(S4 - 2*np.sqrt(2)) < 1e-9 and tv4.max() <= Fstar + 1e-12

print(f"\nF* = (sqrt2-1)/4 = {Fstar:.12f}")
print(f"{'N':>3} {'rem':>3} {'best split':>10} {'max TV_ab':>14} {'margin F*-TV':>14}  verdict")
results = {}
for N in range(4, 25):
    base, rem = divmod(N, 4)
    best = None
    from itertools import combinations
    for chosen in combinations(range(4), rem):
        m = np.array([base + (1 if l in chosen else 0) for l in range(4)], dtype=float)
        tvp = 0.5 * np.abs(p4 - m / N).sum(axis=1)   # per ab row
        mx = tvp.max()
        if best is None or mx < best[0]:
            best = (mx, chosen, tvp)
    mx, chosen, tvp = best
    ok = mx <= Fstar + 1e-12
    results[N] = (ok, mx - Fstar)
    print(f"{N:>3} {rem:>3} {str(chosen):>10} {mx:>14.12f} {Fstar-mx:>+14.6e}  {'PASS' if ok else 'FAIL'}")

npass = sum(1 for ok, _ in results.values() if ok)
print(f"\n{npass}/{len(results)} values of N in [4,24] embed with max TV <= F*")
fails = [N for N, (ok, _) in results.items() if not ok]
if fails:
    print("failing N:", fails)
