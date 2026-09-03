import os as _os; _ROOT = _os.path.normpath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # release: paths relative to repo root
"""TASK 3 (construction half): explicit NS witnesses at F = (sqrt(2)-1)/N, N in {3,4}.

Writes ns_solN{3,4}_{p,u,v}.npy (float64; p shape (4,N) rows ab=00,01,10,11; u,v shape (2,N),
entries +/-1 -- same layout as existing solN*.npy, confirmed by np.load before saving) and the
construction section of ns_witness.out.

Method: loop over all 2^(4N) response vertices; for each solve the NS LP of
phase3.lp_max_S_ns (copied VERBATIM below, additionally returning res.x). Keep the vertex with
maximal S; decompose x: first 4N entries are p[ab,lambda] = x[ab*N:(ab+1)*N]; remaining 8N are
TV auxiliary vars (ignored).

Re-verification from the saved .npy files ALONE is done by a separate fresh process,
ns_verify.py, which appends its PASS/FAIL section to ns_witness.out.
"""
import sys, time
import numpy as np
from scipy.optimize import linprog
from multiprocessing import Pool

WIT = _ROOT + "/witnesses/"
LOG = _ROOT + "/logs/"
sys.path.insert(0, _ROOT + "/certificates")

SIGMA = np.array([[1.0, 1.0], [1.0, -1.0]])
AB = [(a, b) for a in range(2) for b in range(2)]
NPROC = 96


def lp_max_S_ns_x(N, u, v, F0):
    """Verbatim copy of phase3.lp_max_S_ns, additionally returning res.x and res.success."""
    rho = np.full(N, 1.0 / N)
    c = np.zeros(4 * N)
    for ab in range(4):
        a, b = AB[ab]
        c[ab * N:(ab + 1) * N] = SIGMA[a, b] * u[a] * v[b]
    n = 12 * N
    c_full = np.zeros(n); c_full[:4 * N] = c
    nrows = 4 + 4 * N + 4        # row-sums, p-s+t=rho, four NS rows
    A_eq = np.zeros((nrows, n)); b_eq = np.zeros(nrows)
    A_ub = np.zeros((4, n)); b_ub = np.full(4, 2.0 * F0)
    for ab in range(4):
        A_eq[ab, ab * N:(ab + 1) * N] = 1.0; b_eq[ab] = 1.0
        base = 4 + ab * N
        for lam in range(N):
            A_eq[base + lam, ab * N + lam] = 1.0
            A_eq[base + lam, 4 * N + ab * N + lam] = -1.0
            A_eq[base + lam, 8 * N + ab * N + lam] = 1.0
            b_eq[base + lam] = rho[lam]
        A_ub[ab, 4 * N + ab * N:4 * N + (ab + 1) * N] = 1.0
        A_ub[ab, 8 * N + ab * N:8 * N + (ab + 1) * N] = 1.0
    r = 4 + 4 * N
    # AB blocks: 0=(0,0) 1=(0,1) 2=(1,0) 3=(1,1)
    ns_rows = [
        (u[0], 0, 1),   # a=0: sum p[00]u[0] = sum p[01]u[0]
        (u[1], 2, 3),   # a=1: sum p[10]u[1] = sum p[11]u[1]
        (v[0], 0, 2),   # b=0: sum p[00]v[0] = sum p[10]v[0]
        (v[1], 1, 3),   # b=1: sum p[01]v[1] = sum p[11]v[1]
    ]
    for (w, abL, abR) in ns_rows:
        A_eq[r, abL * N:(abL + 1) * N] = w
        A_eq[r, abR * N:(abR + 1) * N] = -w
        b_eq[r] = 0.0
        r += 1
    res = linprog(-c_full, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                  bounds=[(0, 1)] * (4 * N) + [(0, None)] * (8 * N), method='highs')
    if not res.success:
        return -np.inf, None
    return -res.fun, res.x


def vertex_sv(m, N):
    bits = ((m >> np.arange(4 * N)) & 1)
    u = np.where(bits[:2 * N].reshape(2, N) == 0, -1.0, 1.0)
    v = np.where(bits[2 * N:].reshape(2, N) == 0, -1.0, 1.0)
    return u, v


def _work(args):
    N, F0, m0, m1 = args
    best, best_m = -np.inf, None
    for m in range(m0, m1):
        u, v = vertex_sv(m, N)
        s, _ = lp_max_S_ns_x(N, u, v, F0)
        if s > best:
            best, best_m = s, m
    return best, best_m


def ns_residuals(p, u, v):
    """The four NS equality residuals (phase3 ns_rows)."""
    return [
        float(p[0] @ u[0]) - float(p[1] @ u[0]),
        float(p[2] @ u[1]) - float(p[3] @ u[1]),
        float(p[0] @ v[0]) - float(p[2] @ v[0]),
        float(p[1] @ v[1]) - float(p[3] @ v[1]),
    ]


def diagnostics(N, p, u, v, F):
    S = 0.0; tvs = []; rs = []
    for ab in range(4):
        a, b = AB[ab]
        S += SIGMA[a, b] * float(p[ab] @ (u[a] * v[b]))
        tvs.append(0.5 * float(np.abs(p[ab] - 1.0 / N).sum()))
        rs.append(float(p[ab].sum()))
    return S, tvs, rs, ns_residuals(p, u, v)


def run_N(N):
    F = (np.sqrt(2.0) - 1.0) / float(N)
    n = 2 ** (4 * N)
    chunks = np.linspace(0, n, NPROC + 1, dtype=np.int64)
    tasks = [(N, float(F), int(chunks[i]), int(chunks[i + 1])) for i in range(NPROC)
             if chunks[i] < chunks[i + 1]]
    t0 = time.time()
    with Pool(NPROC) as pool:
        results = pool.map(_work, tasks)
    best_S, best_m = max(results, key=lambda r: r[0])
    dt_enum = time.time() - t0

    u, v = vertex_sv(best_m, N)
    S_lp, x = lp_max_S_ns_x(N, u, v, float(F))
    p_raw = np.array(x[:4 * N].reshape(4, N), dtype=np.float64)
    dt_final = time.time() - t0

    # --- decide whether the raw HiGHS point already meets the save tolerances ---
    S_raw, tvs_raw, rs_raw, nsr_raw = diagnostics(N, p_raw, u, v, F)
    tol_tv = F + 1e-9; tol_rs = 1e-9; tol_ns = 1e-8
    raw_ok = (float(p_raw.min()) > -1e-12 and all(t <= tol_tv for t in tvs_raw)
              and all(abs(r - 1.0) <= tol_rs for r in rs_raw)
              and all(abs(q) < tol_ns for q in nsr_raw))

    note = "raw HiGHS point saved as-is (no post-processing)"
    p_save = p_raw
    if not raw_ok:
        # minimal rounding: clip tiny negatives, renormalize rows to sum exactly 1
        p_save = np.clip(p_raw, 0.0, None)
        p_save = p_save / p_save.sum(axis=1, keepdims=True)
        note = ("raw point outside save tolerances -> clipped negatives and renormalized "
                "rows (change ~1e-12); saved the rounded point")

    S_chk, tvs, rs, nsr = diagnostics(N, p_save, u, v, F)
    ok_tv = all(t <= tol_tv for t in tvs)
    ok_rs = all(abs(r - 1.0) <= tol_rs for r in rs)
    ok_ns = all(abs(q) < tol_ns for q in nsr)
    TGT = 2.0 * np.sqrt(2.0)

    with open(LOG + "ns_witness.out", "a") as out:
        def p(*a):
            print(*a, file=out, flush=True)
        p(f"=== TASK 3 (construction): NS witness N={N}, F=(sqrt2-1)/{N}={F:.15f} ===")
        p(f"enumerated all {n} vertices, one NS LP each (Pool nproc={NPROC}); "
          f"enumeration {dt_enum:.1f}s, total {dt_final:.1f}s")
        p(f"winning vertex m={best_m}  S_LP={S_lp:.12f}  (target 2*sqrt(2)={TGT:.12f})")
        p(f"u =\n{np.array2string(u, precision=0)}")
        p(f"v =\n{np.array2string(v, precision=0)}")
        p(f"p[ab,lam] (rows 00,01,10,11) =\n{np.array2string(p_save, precision=12)}")
        p(f"note: {note}")
        p(f"S recomputed from saved arrays = {S_chk:.12f}   S - 2*sqrt(2) = {S_chk-TGT:+.3e}")
        for ab in range(4):
            a, b = AB[ab]
            E = float(p_save[ab] @ (u[a] * v[b]))
            p(f"row {ab} ({a}{b}): E={E:+.12f}  TV={tvs[ab]:.12f} (<= F+1e-9: "
              f"{'ok' if tvs[ab] <= tol_tv else 'VIOLATION'})  rowsum={rs[ab]:.15f}")
        p(f"NS residuals = {[f'{q:+.3e}' for q in nsr]}  (all < 1e-8: "
          f"{'ok' if ok_ns else 'VIOLATION'})")
        max_tv = max(tvs)
        p(f"max per-row TV = {max_tv:.12f} vs F = {F:.12f}  (at least one row ~ F: "
          f"{'yes' if any(abs(t - F) < 1e-6 for t in tvs) else 'no'})")
        p(f"construction PASS/FAIL: S~2sqrt2={abs(S_chk - TGT) < 1e-5}  TV<=F+1e-9={ok_tv}  "
          f"rows=1+/-1e-9={ok_rs}  NS<1e-8={ok_ns}")
        p("")

    np.save(WIT + f"ns_solN{N}_u.npy", u.astype(np.float64))
    np.save(WIT + f"ns_solN{N}_v.npy", v.astype(np.float64))
    np.save(WIT + f"ns_solN{N}_p.npy", p_save.astype(np.float64))
    print(f"N={N}: S={S_chk:.12f} m={best_m} runtime={dt_final:.1f}s "
          f"raw_ok={raw_ok}", flush=True)


if __name__ == "__main__":
    run_N(3)
    run_N(4)
