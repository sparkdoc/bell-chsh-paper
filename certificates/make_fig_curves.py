#!/usr/bin/env python3
import os as _os; _ROOT = _os.path.normpath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # release: paths relative to repo root
"""Generate paper/fig_curves.pdf -- the seven exact S_max(N,F) curves (Theorem 5).

The curve definitions below are taken verbatim from Eq. (eq:curves) of
main.tex, which is certified by audit1.out Part A (exact rational enumeration
over all multisets of even-parity patterns) and cross-checked by the LP
pipeline (capregion.out). Before plotting, this script re-derives every marked
point in exact Fraction arithmetic and asserts it against those certified
values. Exit code 0 iff every assertion passes.

Marked points:
  - Bell point (0, 2) on every curve;
  - bends at F=1/N for N=5,6,7 with certified values 18/5, 10/3, 22/7
    (supp.tex Table II / audit1.out Part A);
  - caps F(S=4;N) = ceil(N/4)/N where each curve reaches S=4 (Theorem 3);
  - Tsirelson crossings F*_N = (sqrt(2)-1)/min(N,4) with S = 2*sqrt(2)
    (Theorem 4 / headline).
"""
import math
import sys
from fractions import Fraction as Fr

SQRT2 = math.sqrt(2)

# Segments as (F_start, intercept, slope, F_end): on [a,b], S = c + m*F.
# Verbatim from main.tex eq:curves; "then 4" beyond the cap is not drawn
# (each curve terminates at its cap point, marked with a diamond).
CURVES = {
    2: [(0, Fr(2), 4, Fr(1, 2))],
    3: [(0, Fr(2), 6, Fr(1, 3))],
    4: [(0, Fr(2), 8, Fr(1, 4))],
    5: [(0, Fr(2), 8, Fr(1, 5)), (Fr(1, 5), Fr(16, 5), 2, Fr(2, 5))],
    6: [(0, Fr(2), 8, Fr(1, 6)), (Fr(1, 6), Fr(8, 3), 4, Fr(1, 3))],
    7: [(0, Fr(2), 8, Fr(1, 7)), (Fr(1, 7), Fr(16, 7), 6, Fr(2, 7))],
    8: [(0, Fr(2), 8, Fr(1, 4))],
}


def S_exact(N, f):
    """S_max(N, f) in exact Fraction arithmetic (f a Fraction)."""
    for (a, c, m, b) in CURVES[N]:
        if a <= f <= b:
            return c + m * f
    raise ValueError("F out of drawn range")


def main():
    checks = []

    def check(name, cond):
        checks.append((name, bool(cond)))
        if not cond:
            print(f"FAIL: {name}")

    # 1. Bell point and continuity at internal junctions (exact).
    for N in range(2, 9):
        check(f"S({N},0)=2", S_exact(N, Fr(0)) == 2)
        segs = CURVES[N]
        for i in range(len(segs) - 1):
            a1, c1, m1, b1 = segs[i]
            a2, c2, m2, b2 = segs[i + 1]
            check(f"S({N},.) continuous at F={b1}",
                  c1 + m1 * b1 == c2 + m2 * a2 and b1 == a2)

    # 2. Caps: stored end of last segment == ceil(N/4)/N, and S = 4 there.
    for N in range(2, 9):
        cap = Fr(-(-N // 4), N)  # ceil(N/4)/N
        check(f"cap({N})={cap} matches Theorem 3", CURVES[N][-1][3] == cap)
        check(f"S({N},cap)=4", S_exact(N, cap) == 4)

    # 3. Certified bend values (supp.tex Table II / audit1.out Part A).
    for N, f, val in [(5, Fr(1, 5), Fr(18, 5)),
                      (6, Fr(1, 6), Fr(10, 3)),
                      (7, Fr(1, 7), Fr(22, 7))]:
        check(f"bend S({N},{f})={val}", S_exact(N, f) == val)

    # 4. Printed cross-check numbers from the paper/logs:
    #    main.tex cap note: S_max(5,1/4)=16/5+2*(1/4)=37/10=3.7;
    #    capregion.out:     S_max(6,0.2)=8/3+4*0.2=52/15.
    check("S(5,1/4)=37/10 (main.tex cap note)", S_exact(5, Fr(1, 4)) == Fr(37, 10))
    check("S(6,1/5)=52/15 (capregion.out)", S_exact(6, Fr(1, 5)) == Fr(52, 15))

    # 5. Tsirelson crossings: F*_N=(sqrt2-1)/min(N,4) lies on the first
    #    segment and S(F*_N)=2*sqrt(2). Exact algebra: on the first segment
    #    S=2+2*min(N,4)*F, so S(F*_N)=2+2*(sqrt2-1)=2*sqrt(2) identically.
    for N in range(2, 9):
        m = min(N, 4)
        fstar = (SQRT2 - 1) / m
        a, c, sl, b = CURVES[N][0]
        check(f"F*({N}) on first segment", a <= fstar <= float(b))
        s_val = float(c + sl * Fr(str(fstar)))  # dense-enough float check
        check(f"S({N},F*)=2*sqrt(2) to 1e-12", abs(s_val - 2 * SQRT2) < 1e-12)

    n_fail = sum(1 for _, ok in checks if not ok)
    print(f"{len(checks)} exact-arithmetic checks, {n_fail} failures")
    if n_fail:
        sys.exit(1)

    # ---- plot ----------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "pdf.fonttype": 42,   # TrueType fonts (journal requirement)
        "ps.fonttype": 42,
        "font.size": 8,
        "axes.labelsize": 9,
        "legend.fontsize": 7,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
    })

    fig, ax = plt.subplots(figsize=(3.5, 2.9))

    for N in range(2, 9):
        xs, ys = [0.0], [2.0]
        for (a, c, m, b) in CURVES[N]:
            xs += [float(a), float(b)]
            ys += [float(c + m * a), float(c + m * b)]
        ax.plot(xs, ys, color=f"C{N-2}", lw=(2.4 if N == 4 else 1.3), ls=("--" if N == 4 else "-"), zorder=(1 if N == 4 else 2), label=f"$N={N}$")
        # cap marker (curve flattens at S=4 beyond this point)
        cap = CURVES[N][-1][3]
        ax.plot([float(cap)], [4.0], marker="D", ms=3.2, mfc=f"C{N-2}",
                mec="k", mew=0.4, ls="none")

    # Tsirelson line
    ax.axhline(2 * SQRT2, color="0.35", lw=0.8, ls="--", zorder=1)
    ax.text(0.495, 2 * SQRT2 + 0.045, "$2\\sqrt{2}$ (Tsirelson)",
            ha="right", va="bottom", fontsize=7)

    # F* crossings: N=2 and N=3 distinct; all N>=4 coincide at (sqrt2-1)/4.
    for N in (2, 3):
        fstar = (SQRT2 - 1) / N
        ax.plot([fstar], [2 * SQRT2], "o", ms=3.0, mfc=f"C{N-2}",
                mec="k", mew=0.4, ls="none")
    fstar4 = (SQRT2 - 1) / 4
    ax.plot([fstar4], [2 * SQRT2], "*", ms=7, mfc="0.15", mec="k",
            mew=0.4, ls="none", zorder=5)
    ax.annotate("$F^*=(\\sqrt{2}-1)/4$ (all N $\\geq$ 4)",
                xy=(fstar4, 2 * SQRT2), xytext=(fstar4 + 0.035, 2.45),
                fontsize=7, arrowprops=dict(arrowstyle="-", lw=0.5))

    # bend markers (open circles) for N=5,6,7
    for N, f, val in [(5, Fr(1, 5), Fr(18, 5)),
                      (6, Fr(1, 6), Fr(10, 3)),
                      (7, Fr(1, 7), Fr(22, 7))]:
        ax.plot([float(f)], [float(val)], "o", ms=3.4, mfc="w",
                mec=f"C{N-2}", mew=0.9, ls="none")
    ax.annotate("bends", xy=(1 / 5, 18 / 5), xytext=(0.21, 3.0),
                fontsize=7, arrowprops=dict(arrowstyle="-", lw=0.5))

    # Bell point annotation
    ax.text(0.005, 2.06, "Bell bound $S=2$", ha="left", va="bottom",
            fontsize=7)

    ax.set_xlim(0, 0.5)
    ax.set_ylim(1.95, 4.18)
    ax.set_xlabel("$F$ (worst-case TV fine-tuning)")
    ax.set_ylabel("$S_{\\max}(N,F)$")
    ax.grid(True, lw=0.3, alpha=0.5)
    ax.legend(loc="upper left", framealpha=0.9, handlelength=1.6)

    fig.tight_layout()
    fig.savefig(_ROOT + "/paper/fig_curves.pdf")


if __name__ == "__main__":
    main()
