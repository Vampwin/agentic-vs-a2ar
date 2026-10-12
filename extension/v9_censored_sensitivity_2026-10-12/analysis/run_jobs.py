"""ADR-AIIM-412: run every job of the pre-registration in parallel processes (one thread each, except the three-thread determinism job).
Each job writes results/<name>.npz and results/<name>.json; a finished job is skipped on a second call."""
from __future__ import annotations

import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

sys.path.insert(0, str(Path(__file__).parent))
import numpy as np                       # noqa: E402
import pandas as pd                      # noqa: E402
from sklearn.ensemble import RandomForestRegressor   # noqa: E402
from sklearn.metrics import mean_squared_error, r2_score   # noqa: E402

import lib                               # noqa: E402

RES = lib.OUT / "results"
ANCHOR = "D1sorted_D2trainfirst_D3divert"


def save(name, arrays, meta):
    RES.mkdir(exist_ok=True)
    np.savez_compressed(RES / f"{name}.npz", **arrays)
    (RES / f"{name}.json").write_text(json.dumps(meta, indent=1), encoding="utf-8")


def job_em(name, gi, init, max_iter, tol, n_jobs):
    d = lib.data()
    g = lib.groups()[gi]
    tr, _ = lib.split_for(g[0])
    t0 = time.time()
    r = lib.em_passes(d["X"], d["y0"], d["kind"], tr, d["X_A"], init, max_iter, tol, n_jobs)
    P = np.vstack([x["pred"] for x in r["recs"]])
    meta = {"group": g, "init": init, "max_iter": max_iter, "tol": tol, "n_jobs": n_jobs, "n_pass": r["n_pass"],
            "stopped_by_tol": r["stopped_by_tol"], "final_change": r["final_change"],
            "max_change": [x["max_change"] for x in r["recs"]], "mean_change": [x["mean_change"] for x in r["recs"]],
            "n_moved": [x["n_moved_gt_0.01"] for x in r["recs"]], "seconds": round(time.time() - t0, 1)}
    save(name, {"pass_pred": P, "final_pred": r["final_pred"]}, meta)
    return name, meta["seconds"]


def job_tobit(name, gi):
    d = lib.data()
    g = lib.groups()[gi]
    tr, te = lib.split_for(g[0])
    X, y0, kind = d["X"], d["y0"], d["kind"]
    codes = pd.factorize(pd.Series(d["scaf"][tr]))[0]
    t0 = time.time()
    m = lib.RidgeTobit().fit(X[tr], y0[tr], kind[tr], codes)
    te_pt, te_up, te_lo = te[kind[te] == 0], te[kind[te] == 1], te[kind[te] == -1]
    p_pt = m.predict(X[te_pt])
    pu = m.predict(X[te_up]) if len(te_up) else np.array([])      # a split can hold no censored record of one kind
    pl = m.predict(X[te_lo]) if len(te_lo) else np.array([])
    pool = m.predict(d["X_A"])
    meta = {"group": g, "lambda": m.lam, "cv_nll": {str(k): v for k, v in m.cv.items()}, "sigma": m.sigma,
            "converged": bool(m.res.success), "n_iter": int(m.res.nit), "message": str(m.res.message),
            "n_features": int(m.d), "n_train": int(len(tr)),
            "heldout_point_spearman": lib.sp(y0[te_pt], p_pt),
            "heldout_point_r2": float(r2_score(y0[te_pt], p_pt)),
            "heldout_point_rmse": float(np.sqrt(mean_squared_error(y0[te_pt], p_pt))),
            "upper_n": int(len(te_up)), "upper_pass_frac": float((pu <= y0[te_up]).mean()) if len(te_up) else None,
            "lower_n": int(len(te_lo)), "lower_pass_frac": float((pl >= y0[te_lo]).mean()) if len(te_lo) else None,
            "pool_min": float(pool.min()), "pool_max": float(pool.max()), "seconds": round(time.time() - t0, 1)}
    save(name, {"pool_pred": pool}, meta)
    return name, meta["seconds"]


def job_sim(name, truth, seed):
    d = lib.data()
    X, y0, kind = d["X"], d["y0"], d["kind"]
    n = len(y0)
    rng = np.random.default_rng(1000 + seed)
    t0 = time.time()
    if truth == "lin":
        var_cols = np.where(X.var(axis=0) >= 0.01)[0]
        cols = rng.choice(var_cols, 100, replace=False)
        w = rng.normal(size=100)
        z = X[:, cols] @ w
    else:
        pts = np.where(kind == 0)[0]
        rf = RandomForestRegressor(random_state=seed, n_jobs=1).fit(X[pts], y0[pts])
        z = rf.predict(X)
    lat = 6.3 + 1.3 * (z - z.mean()) / z.std()
    v = lat + rng.normal(0, 0.35, n)
    n_up, n_lo = int((kind == 1).sum()), int((kind == -1).sum())
    lo_b = float(np.percentile(v, 100 * n_up / n))
    hi_b = float(np.percentile(v, 100 * (1 - n_lo / n)))
    ks, ys = np.zeros(n, dtype=int), v.copy()
    ks[v < lo_b], ys[v < lo_b] = 1, lo_b
    ks[v > hi_b], ys[v > hi_b] = -1, hi_b
    tr, te = lib.split_for(ANCHOR)
    out = {}
    trp = tr[ks[tr] == 0]
    out["drop"] = lib.pipe(1).fit(X[trp], ys[trp]).predict(X[te])
    out["bound"] = lib.pipe(1).fit(X[tr], ys[tr]).predict(X[te])
    out["em"] = lib.em_passes(X, ys, ks, tr, X[te], "bound", 10, 0.01, 1)["final_pred"]
    codes = pd.factorize(pd.Series(d["scaf"][tr]))[0]
    tb = lib.RidgeTobit().fit(X[tr], ys[tr], ks[tr], codes)
    out["tobit"] = tb.predict(X[te])
    below = lat[te] < lo_b
    meta = {"truth": truth, "seed": seed, "lo_bound": lo_b, "hi_bound": hi_b, "n_te": int(len(te)),
            "n_te_below": int(below.sum()), "tobit_lambda": tb.lam, "seconds": round(time.time() - t0, 1)}
    for k, p in out.items():
        meta[k] = {"rmse": lib.rmse(p, lat[te]), "spearman": lib.sp(p, lat[te]),
                   "min_pred_among_below": float(p[below].min()) if below.any() else None,
                   "share_below_bound_predicted_below": float((p[below] < lo_b).mean()) if below.any() else None}
    save(name, {f"pred_{k}": p for k, p in out.items()} | {"lat_te": lat[te]}, meta)
    return name, meta["seconds"]


def jobs():
    J = []
    for gi in range(5):                                           # C2 and the anchor: 30 passes, no tolerance stop
        J.append((f"c2_g{gi}_bound_30", job_em, (f"c2_g{gi}_bound_30", gi, "bound", 30, 0.0, 1)))
    for truth in ("lin", "rf"):                                   # C1
        for seed in range(1, 6):
            J.append((f"c1_{truth}_{seed}", job_sim, (f"c1_{truth}_{seed}", truth, seed)))
    for gi in range(5):                                           # C3, registered passes and tolerance
        for init in ("inside1", "random", "mean"):
            J.append((f"c3_g{gi}_{init}", job_em, (f"c3_g{gi}_{init}", gi, init, 10, 0.01, 1)))
    for tag, nj in (("a", 1), ("b", 1), ("c", 3)):                # C5
        J.append((f"c5_{tag}", job_em, (f"c5_{tag}", 4, "bound", 10, 0.01, nj)))
    for gi in range(5):                                           # C4
        J.append((f"c4_g{gi}", job_tobit, (f"c4_g{gi}", gi)))
    return J


def main():
    RES.mkdir(exist_ok=True)
    todo = [(n, f, a) for n, f, a in jobs() if not (RES / f"{n}.json").exists()]
    print(f"{len(todo)} jobs to run", flush=True)
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=16) as ex:
        futs = {ex.submit(f, *a): n for n, f, a in todo}
        for fu in as_completed(futs):
            try:
                n, s = fu.result()
                print(f"[{time.time()-t0:6.0f}s] done {n} ({s}s)", flush=True)
            except Exception as e:           # a failed job is reported, not hidden
                print(f"[{time.time()-t0:6.0f}s] FAILED {futs[fu]}: {e!r}", flush=True)
    print("all submitted jobs finished", flush=True)


if __name__ == "__main__":
    main()
