import os as _os; _ROOT = _os.path.normpath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # release: paths relative to repo root
"""TASK 2: direct NS tightness evaluations at F = (sqrt(2)-1)/N for N=3,4.

Writes ns_tight_N3.out and ns_tight_N4.out. Context (date-independent):
  - NS variant: phase3.exact_Smax_ns(N,F) enumerates all 2^(4N) response vertices and solves
    the LP with the four no-signalling equality rows added (phase3.py lines 43-53).
  - Since S_NS(N,F) <= S_max(N,F) always, and S_max(3,(sqrt2-1)/3)=S_max(4,(sqrt2-1)/4)=2*sqrt(2)
    exactly (PROOFS.md Theorem 3), a direct evaluation with S >= 2*sqrt(2) - 1e-5 settles
    F*_NS(N) = (sqrt2-1)/N EXACTLY (the curve is strictly increasing up to F*).
  - Earlier bisections (13 iters from [0,0.5], resolution ~6.1e-5 in F) reported
    F*_NS(3)~0.13812, F*_NS(4)~0.10358 -- slightly above (sqrt2-1)/3=0.13807119 and
    (sqrt2-1)/4=0.10355339. If the direct S here equals 2*sqrt(2) to ~1e-6, those offsets
    were pure bisection-resolution artifacts.
  - Fallback if S is below 2*sqrt(2) by more than 1e-5 (a REAL NS effect): bisect F*_NS(N)
    on [F_N, F_N+0.001] for 40 iterations and print the final bracket.
  - LP tolerance: HiGHS ~1e-9..1e-7; decision threshold 1e-5 as specified.
"""
import sys, time
import numpy as np

sys.path.insert(0, _ROOT + "/certificates")
from phase3 import exact_Smax_ns

D = _ROOT + "/logs/"
NPROC = 96
TGT = 2.0 * np.sqrt(2.0)
THRESH = 1e-5

def run_N(N):
    F = (np.sqrt(2.0) - 1.0) / float(N)
    out = open(D + f"ns_tight_N{N}.out", "w")
    def p(*a):
        print(*a, file=out, flush=True)

    p(f"=== TASK 2: direct NS tightness evaluation, N={N} ===")
    p(f"F = (sqrt(2)-1)/{N} = {F:.15f}")
    p(f"S_NS(N,F) via phase3.exact_Smax_ns (Pool nproc={NPROC}, HiGHS); "
      f"target 2*sqrt(2) = {TGT:.15f}")
    p("")
    t0 = time.time()
    S = exact_Smax_ns(N, float(F), nproc=NPROC)
    dt = time.time() - t0
    diff = S - TGT
    p(f"S_NS({N}, F) = {S:.12f}")
    p(f"S - 2*sqrt(2) = {diff:+.12e}")
    p(f"runtime = {dt:.1f} s   ({2**(4*N)} vertices, one NS LP each)")
    p("")
    if S >= TGT - THRESH:
        # S_NS <= S_max = 2*sqrt(2) exactly, so equality holds to within LP tolerance
        p(f"VERDICT: S_NS({N},(sqrt2-1)/{N}) = 2*sqrt(2) to within LP tolerance "
          f"(>= target - {THRESH:g}).")
        p(f"=> F*_NS({N}) = (sqrt2-1)/{N} EXACTLY. The earlier bisection offsets were pure")
        p("   bisection-resolution artifacts (13 iters from [0,0.5], resolution ~6.1e-5).")
        p(f"RESULT: PASS")
    else:
        p(f"VERDICT: S is below 2*sqrt(2) by {-diff:.3e} > {THRESH:g}: REAL NS effect.")
        p("Bisecting F*_NS(N) on [F_N, F_N+0.001], 40 iterations:")
        lo, hi = float(F), float(F) + 0.001
        for it in range(40):
            mid = 0.5 * (lo + hi)
            S_mid = exact_Smax_ns(N, mid, nproc=NPROC)
            p(f"  it{it:2d} F={mid:.12f}  S={S_mid:.12f}  {'hi=' if S_mid >= TGT - 1e-6 else 'lo='}",
              flush=True)
            if S_mid >= TGT - 1e-6:
                hi = mid
            else:
                lo = mid
        p(f"final bracket: F*_NS({N}) in [{lo:.12f}, {hi:.12f}]")
        p("RESULT: FAIL vs exact (sqrt2-1)/N -- genuine NS gap; see bracket above.")
    out.close()
    print(f"N={N}: S={S:.12f}  diff={diff:+.3e}  runtime={dt:.1f}s", flush=True)

if __name__ == "__main__":
    run_N(3)
    run_N(4)
