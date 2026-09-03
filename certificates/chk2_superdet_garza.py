#!/usr/bin/python3
import os as _os; _ROOT = _os.path.normpath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # release: paths relative to repo root
"""chk2_superdet_garza.py — apply the project's Proposition (S <= 2+8F for ANY local CHSH model,
any source rho, deterministic or stochastic responses; PROOFS.md Sec.5) to the EXTERNAL explicit
model of Garza & Hance 2025 (arXiv:2511.06043), "Quantum-Like Correlations from Local
Hidden-Variable Theories Under Conservation Law".

Their model (verified against the paper, literature/garza_hance_2025.pdf):
  - lambda = angle in [0, pi]; baseline rho(lambda) = 1/pi uniform, explicitly independent of
    settings ("rho(lambda) remains independent of the measurement settings", Sec. IV).
  - Conservation-law measurement error excludes a band of half-width b = sin(theta)/(2*dL)
    centered at pi-theta from Bob's lambda domain (theta = angle between the two settings);
    observed density is uniform over the remainder: p(lambda|ab) = 1/N on allowed set, 0 on the
    band, N = pi - sin(theta)/dL.  Deterministic local responses (sign functions).
  - Correlator (their Eq. 13):  E(theta; dL) = dL*(2*theta - pi)/(pi*dL - sin(theta)), theta in [0,pi].

Setting geometry: two parameterizations are used and both reported.
  (P1) LITERAL box reading of the paper (alpha, beta in [0,pi]): settings a0,a1,b0,b1 in [0,pi],
       pairwise angle theta_ij = |ai - bj|.
  (P2) Oriented-directions (circle) reading: a0 = 0 (rotational invariance), a1,b0,b1 in [0,2pi),
       pairwise angle theta_ij = circular distance in [0,pi].  This admits the standard
       Tsirelson configuration {45,45,45,135} deg, which is infeasible under (P1).

This is a floating-point ANALYSIS of an external model (not a certificate of the project's
theorems): theta values are transcendental, so exact Fraction/Q(sqrt(2)) arithmetic does not
apply.  All checks use double precision with margins >> 1e-9, and each analytic identity is
cross-checked by independent numerical integration.

Exit 0 iff every assertion passes.
"""
import math
import sys

try:
    import numpy as np
except ImportError:
    print("numpy required")
    sys.exit(2)

LOG = open(_ROOT + "/logs/out_superdet_garza.txt", "w")


def out(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.write(s + "\n")


PI = math.pi


def E(theta, dL):
    """Garza & Hance Eq. (13)."""
    return dL * (2.0 * theta - PI) / (PI * dL - np.sin(theta))


def circ_dist(x, y):
    """Circular distance in [0, pi] between angles x, y in [0, 2pi)."""
    d = np.abs(x - y) % (2.0 * PI)
    return np.minimum(d, 2.0 * PI - d)


def tv_single_band_numeric(theta, dL, n=2_000_001):
    """Direct numerical TV(p(lambda|ab)||rho) for ONE excluded band of half-width
    b = sin(theta)/(2 dL) centered at pi-theta; p uniform 1/N on the complement."""
    s = math.sin(theta) / dL          # band length
    N = PI - s                        # allowed measure
    if N <= 0:
        return None
    b = s / 2.0
    c = PI - theta                    # band center
    x = (np.arange(n) + 0.5) * (PI / n)
    p = np.where((x > c - b) & (x < c + b), 0.0, 1.0 / N)
    return float(0.5 * np.sum(np.abs(p - 1.0 / PI)) * (PI / n))


def tv_mixture_numeric(theta, dL, n=2_000_001):
    """TV for the equal-weight mixture of the two outcome-conditional bands (centered at
    pi-theta [Alice +] and theta [Alice -])."""
    s = math.sin(theta) / dL
    N = PI - s
    if N <= 0:
        return None
    b = s / 2.0
    c1, c2 = PI - theta, theta
    x = (np.arange(n) + 0.5) * (PI / n)
    in1 = (x > c1 - b) & (x < c1 + b)
    in2 = (x > c2 - b) & (x < c2 + b)
    p = 0.5 * ((np.where(in1, 0.0, 1.0 / N)) + (np.where(in2, 0.0, 1.0 / N)))
    return float(0.5 * np.sum(np.abs(p - 1.0 / PI)) * (PI / n))


def main():
    ok = True

    # ---- 1. Analytic single-band TV identity: TV = sin(theta)/(pi dL) -------------------
    out("== 1. Single-band TV identity (analytic vs direct numerical integration) ==")
    for theta, dL in [(PI / 4, 0.77), (3 * PI / 4, 0.77), (PI / 3, 1.0)]:
        num = tv_single_band_numeric(theta, dL)
        ana = math.sin(theta) / (PI * dL)
        err = abs(num - ana)
        out(f"theta={math.degrees(theta):6.1f} deg, dL={dL}: TV_numeric={num:.10f}  "
            f"TV_analytic=sin(th)/(pi dL)={ana:.10f}  err={err:.2e}")
        if err > 1e-6:
            ok = False
    out()

    # ---- 2. Mixture TV at the standard CHSH angles, dL = 0.77 --------------------------
    out("== 2. Mixture (outcome-averaged) TV at dL=0.77 ==")
    for theta in (PI / 4, 3 * PI / 4):
        m = tv_mixture_numeric(theta, 0.77)
        ub = math.sin(theta) / (PI * 0.77)
        out(f"theta={math.degrees(theta):6.1f} deg: TV_mixture={m:.10f}   "
            f"single-band upper bound sin/(pi dL)={ub:.10f}")
    # analytic closed form for the two-band mixture, valid when the bands are
    # disjoint (|pi - 2 theta| >= s) and s < pi/2: with band length
    # s = sin(theta)/dL and N = pi - s, p = 1/N outside both bands and 1/(2N)
    # inside exactly one, so
    #   TV = (1/2)[ 2s |1/(2N) - 1/pi| + (pi - 2s) |1/N - 1/pi| ]
    #      = (1/2)[ 2s (pi - 2s)/(2 pi N) + (pi - 2s) s/(pi N) ]
    #      = s (pi - 2s) / (pi (pi - s)).
    for theta in (PI / 4, 3 * PI / 4):
        dL = 0.77
        s = math.sin(theta) / dL
        if not (s < PI / 2 and abs(PI - 2.0 * theta) >= s):
            out(f"theta={math.degrees(theta):6.1f} deg: CLOSED-FORM PRECONDITIONS FAIL "
                f"(s={s:.8f}) — mixture TV identity not applicable")
            ok = False
            continue
        ana_mix = s * (PI - 2.0 * s) / (PI * (PI - s))
        num_mix = tv_mixture_numeric(theta, dL)
        err = abs(num_mix - ana_mix)
        out(f"theta={math.degrees(theta):6.1f} deg: TV_mixture_analytic="
            f"s(pi-2s)/(pi(pi-s))={ana_mix:.10f}  err_vs_numeric={err:.2e}")
        if err > 1e-6:
            ok = False
    out()

    # ---- 3. Self-consistency: |E(theta)| <= 1 for all theta? ---------------------------
    out("== 3. Self-consistency of their correlator: max_theta |E(theta; dL)| ==")
    th = np.linspace(0.0, PI, 200_001)
    for dL in (0.3, 0.4, 0.45, 0.5, 0.77, 1.0, 2.0):
        Em = E(th, dL)
        out(f"dL={dL:5.2f}: max|E|={np.abs(Em).max():.8f} "
            f"({'OK (|E|<=1)' if np.abs(Em).max() <= 1 + 1e-9 else 'UNPHYSICAL (|E|>1 somewhere)'}); "
            f"E(45deg)={E(np.array([PI/4]), dL)[0]:+.6f}, E(135deg)={E(np.array([3*PI/4]), dL)[0]:+.6f}")
    out()

    # ---- 4. Standard Tsirelson configuration {45,45,45,135} (P2 geometry) vs dL --------
    out("== 4. S at the standard Tsirelson config (a0=0,a1=90,b0=45,b1=-45 deg; "
        "thetas {45,45,45,135}) ==")
    a0v, a1v, b0v, b1v = 0.0, PI / 2, PI / 4, -PI / 4
    ths = [circ_dist(a0v, b0v), circ_dist(a0v, b1v), circ_dist(a1v, b0v), circ_dist(a1v, b1v)]
    for dL in (0.5, 0.77, 1.0, 2.0, 5.0):
        S = E(ths[0], dL) + E(ths[1], dL) + E(ths[2], dL) - E(ths[3], dL)
        Fub = max(np.sin(t) for t in ths) / (PI * dL)
        margin = 2.0 + 8.0 * Fub - abs(S)
        out(f"dL={dL:5.2f}: S={S:+.8f}  |S|={abs(S):.8f}  F(upper bound)={Fub:.8f}  "
            f"margin(2+8F-|S|)={margin:.6f}  |S|-2sqrt2={abs(S)-2*math.sqrt(2):+.6f}")
        if margin < -1e-9:
            ok = False
    out()

    # ---- 5. (P1) box grid: max |S| and Proposition check, theta_ij = |ai-bj| -----------
    out("== 5. (P1) LITERAL box [0,pi]^4, theta=|ai-bj|: max |S| + Proposition check ==")
    n = 181
    ang = np.linspace(0.0, PI, n)
    worst_margin = float("inf")
    for dL in (0.5, 0.77, 1.0, 2.0):
        best_S, best_cfg = -1.0, None
        # separable: S = [f|a0-b0|+f|a1-b0|] + [f|a0-b1|-f|a1-b1|]; also a0<->a1 relabel S'
        for i in range(n):
            d0 = np.abs(ang[i] - ang)
            f0 = E(d0, dL)
            for j in range(i + 1, n):
                d1 = np.abs(ang[j] - ang)
                f1 = E(d1, dL)
                gmax, gmin = (f0 + f1).max(), (f0 + f1).min()
                hmax, hmin = (f0 - f1).max(), (f0 - f1).min()
                cands = ((gmax + hmax, 0, 0), (gmin + hmin, 1, 1),
                         (gmax - hmin, 0, 1), (gmin - hmax, 1, 0))
                for Sij, mk, ml in cands:
                    if abs(Sij) > best_S:
                        best_S = abs(Sij)
                        k = int(np.argmax(f0 + f1)) if mk == 0 else int(np.argmin(f0 + f1))
                        l = int(np.argmax(f0 - f1)) if ml == 0 else int(np.argmin(f0 - f1))
                        best_cfg = (ang[i], ang[j], ang[k], ang[l])
        a0v, a1v, b0v, b1v = best_cfg
        F = max(math.sin(abs(a0v - b0v)), math.sin(abs(a0v - b1v)),
                math.sin(abs(a1v - b0v)), math.sin(abs(a1v - b1v))) / (PI * dL)
        cfgdeg = tuple(round(math.degrees(x), 2) for x in best_cfg)
        out(f"dL={dL:5.2f}: max|S|={best_S:.8f} at (a0,a1,b0,b1)={cfgdeg} deg; "
            f"F={F:.8f}; floor(S-2)/8={(best_S-2.0)/8:.8f}")
        # full-grid Proposition check, chunked over (a0,a1) pairs
        for i in range(n):
            d0 = np.sin(np.abs(ang[i] - ang))
            f0 = E(np.abs(ang[i] - ang), dL)
            for j in range(i + 1, n):
                d1 = np.sin(np.abs(ang[j] - ang))
                f1 = E(np.abs(ang[j] - ang), dL)
                Smat = f0[:, None] + f0[None, :] + f1[:, None] - f1[None, :]
                # a0<->a1 relabel (covers a0>a1 quadruples): S(a1,a0,b0,b1)
                Smatp = (f0 + f1)[:, None] + (f1 - f0)[None, :]
                Fmat = np.maximum(np.maximum(d0[:, None], d0[None, :]),
                                  np.maximum(d1[:, None], d1[None, :])) / (PI * dL)
                m = min((2.0 + 8.0 * Fmat - np.abs(Smat)).min(),
                        (2.0 + 8.0 * Fmat - np.abs(Smatp)).min())
                if m < worst_margin:
                    worst_margin = m
        out(f"   P1 Proposition margin over full grid at dL={dL}: {worst_margin:.6e}")
    out(f"WORST P1 margin so far: {worst_margin:.6e} (must be >= -1e-9)")
    if worst_margin < -1e-9:
        ok = False
    out()

    # ---- 6. (P2) circle grid: max |S| and Proposition check ----------------------------
    out("== 6. (P2) oriented directions on circle, a0=0, thetas=circular distance ==")
    n2 = 91
    angc = np.linspace(0.0, 2.0 * PI, n2, endpoint=False)
    worst_margin2 = float("inf")
    for dL in (0.5, 0.77, 1.0):
        best_S, best_cfg = -1.0, None
        for i in range(n2):          # a1
            for k in range(n2):      # b0
                th_a1b0 = float(circ_dist(angc[i], angc[k]))
                f_a1b0 = E(np.array([th_a1b0]), dL)[0]
                d_ab0 = circ_dist(0.0, angc[k])             # theta(0,b0) over b1? no: fixed b0
                th_0b0 = float(circ_dist(0.0, angc[k]))
                f_0b0 = E(np.array([th_0b0]), dL)[0]
                db1 = circ_dist(angc[i], angc)              # |a1 - b1| over b1
                f_a1b1 = E(db1, dL)
                d0b1 = circ_dist(0.0, angc)                 # |0 - b1| over b1
                f_0b1 = E(d0b1, dL)
                Svec = (f_0b0 + f_a1b0) + f_0b1 - f_a1b1    # over b1, a0=0 first
                Svec2 = (f_0b0 + f_a1b0) - f_0b1 + f_a1b1   # Alice labels swapped
                for Sv in (Svec, Svec2):
                    imax, imin = int(np.argmax(Sv)), int(np.argmin(Sv))
                    for idx in (imax, imin):
                        Sij = float(Sv[idx])
                        if abs(Sij) > best_S:
                            best_S = abs(Sij)
                            best_cfg = (0.0, angc[i], angc[k], angc[idx])
        a0v, a1v, b0v, b1v = best_cfg
        ths = [float(circ_dist(a0v, b0v)), float(circ_dist(a0v, b1v)),
               float(circ_dist(a1v, b0v)), float(circ_dist(a1v, b1v))]
        F = max(math.sin(t) for t in ths) / (PI * dL)
        cfgdeg = tuple(round(math.degrees(x % (2 * PI)), 2) for x in best_cfg)
        thdeg = tuple(round(math.degrees(t), 2) for t in ths)
        out(f"dL={dL:5.2f}: max|S|={best_S:.8f} at (a0,a1,b0,b1)={cfgdeg} deg, "
            f"thetas={thdeg}; F={F:.8f}; floor(S-2)/8={(best_S-2.0)/8:.8f}")
        # full-grid Proposition check on the n2^3 grid (a1,b0,b1), vectorized over b1
        for i in range(n2):
            db1 = circ_dist(angc[i], angc)
            f_a1b1 = E(db1, dL)
            s_a1b1 = np.sin(db1)
            d0b1 = circ_dist(0.0, angc)
            f_0b1 = E(d0b1, dL)
            s_0b1 = np.sin(d0b1)
            for k in range(n2):
                th_a1b0 = float(circ_dist(angc[i], angc[k]))
                f_a1b0 = E(np.array([th_a1b0]), dL)[0]
                s_a1b0 = math.sin(th_a1b0)
                th_0b0 = float(circ_dist(0.0, angc[k]))
                f_0b0 = E(np.array([th_0b0]), dL)[0]
                s_0b0 = math.sin(th_0b0)
                Svec = (f_0b0 + f_a1b0) + f_0b1 - f_a1b1
                Svec2 = (f_0b0 + f_a1b0) - f_0b1 + f_a1b1   # Alice labels swapped
                Fvec = np.maximum(np.maximum(s_0b0, s_a1b0), np.maximum(s_0b1, s_a1b1)) / (PI * dL)
                m = min((2.0 + 8.0 * Fvec - np.abs(Svec)).min(),
                        (2.0 + 8.0 * Fvec - np.abs(Svec2)).min())
                if m < worst_margin2:
                    worst_margin2 = m
        out(f"   P2 Proposition margin over full grid at dL={dL}: {worst_margin2:.6e}")
    out(f"WORST P2 margin: {worst_margin2:.6e} (must be >= -1e-9)")
    if worst_margin2 < -1e-9:
        ok = False
    out()

    # ---- 7. Headline numbers at their best-fit dL ~ 0.77 --------------------------------
    out("== 7. Headline at dL=0.77 (their best fit to QM, max disagreement ~0.03) ==")
    a0v, a1v, b0v, b1v = 0.0, PI / 2, PI / 4, -PI / 4
    ths = [float(circ_dist(a0v, b0v)), float(circ_dist(a0v, b1v)),
           float(circ_dist(a1v, b0v)), float(circ_dist(a1v, b1v))]
    S = E(ths[0], 0.77) + E(ths[1], 0.77) + E(ths[2], 0.77) - E(ths[3], 0.77)
    Fub = max(math.sin(t) for t in ths) / (PI * 0.77)
    Fmix45 = tv_mixture_numeric(PI / 4, 0.77)
    Fmix135 = tv_mixture_numeric(3 * PI / 4, 0.77)
    out(f"thetas (deg) = {[round(math.degrees(t),2) for t in ths]}")
    out(f"S={S:+.8f}   (|S|={abs(S):.8f};  2sqrt2={2*math.sqrt(2):.8f})")
    out(f"F (single-band upper bound) = {Fub:.8f}  = {100*Fub:.2f}%")
    out(f"F (true mixture, 45-deg block) = {Fmix45:.8f};  (135-deg block) = {Fmix135:.8f}")
    Ftrue = max(Fmix45, Fmix135)
    out(f"Proposition margin with TRUE mixture F: 2+8F-|S| = {2+8*Ftrue-abs(S):.6f}")
    if 2 + 8 * Ftrue - abs(S) < -1e-9:
        ok = False
    floor = (math.sqrt(2) - 1) / 4
    out(f"project floor F*=(sqrt2-1)/4 = {floor:.8f} = {100*floor:.2f}%")
    out(f"excess fine-tuning: F/F* = {Ftrue/floor:.3f} (mixture), {Fub/floor:.3f} (upper bound)")
    out()

    # ---- 8. Model validity range: N = pi - sin(theta)/dL must stay positive -------------
    out("== 8. Validity: allowed measure N = pi - sin(theta)/dL > 0 ==")
    for dL in (0.2, 0.23, 0.3, 0.5, 0.77):
        N45 = PI - math.sin(PI / 4) / dL
        N135 = PI - math.sin(3 * PI / 4) / dL
        out(f"dL={dL}: N(45deg)={N45:.6f}, N(135deg)={N135:.6f} "
            f"({'OK' if min(N45, N135) > 0 else 'BREAKS DOWN'})")
    out()

    out("ALL PASS" if ok else "FAILURES PRESENT")
    LOG.close()
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
