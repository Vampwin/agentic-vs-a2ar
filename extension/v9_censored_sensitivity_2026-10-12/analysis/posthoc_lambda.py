"""ADR-AIIM-413, POST HOC (not pre-registered): the pre-registered penalty grid of C4 ended at its upper edge (10,000 chosen in all
five distinct splits). The grid is extended to 1e5, 1e6 and 1e7 and the C4 quantities recomputed; reported as post hoc."""
from __future__ import annotations
import json, os, sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
sys.path.insert(0, str(Path(__file__).parent))
import numpy as np, pandas as pd
from sklearn.metrics import r2_score
import lib

def job(gi):
    d = lib.data(); g = lib.groups()[gi]; tr, te = lib.split_for(g[0])
    X, y0, kind = d["X"], d["y0"], d["kind"]
    codes = pd.factorize(pd.Series(d["scaf"][tr]))[0]
    m = lib.RidgeTobit(lambdas=(1e3, 1e4, 1e5, 1e6, 1e7)).fit(X[tr], y0[tr], kind[tr], codes)
    pt, up = te[kind[te] == 0], te[kind[te] == 1]
    p = m.predict(X[pt]); pu = m.predict(X[up])
    return gi, {"group": g, "lambda": m.lam, "cv_nll": {str(k): v for k, v in m.cv.items()}, "sigma": m.sigma,
                "converged": bool(m.res.success), "heldout_point_spearman": lib.sp(y0[pt], p),
                "heldout_point_r2": float(r2_score(y0[pt], p)), "upper_pass_frac": float((pu <= y0[up]).mean()),
                "pool_pred": m.predict(d["X_A"]).tolist()}

if __name__ == "__main__":
    with ProcessPoolExecutor(5) as ex:
        res = dict(ex.map(job, range(5)))
    out = {"note": "POST HOC, not pre-registered", "per_split": {}}
    by = {}
    for gi, r in res.items():
        by[lib.groups()[gi][0]] = np.array(r["pool_pred"])
        out["per_split"][lib.groups()[gi][0]] = {k: v for k, v in r.items() if k != "pool_pred"}
    d = lib.data(); uid, mur = list(d["uid"]), list(d["murcko"])
    rk, sl = lib.select(lib.expand(by), uid, mur)
    prev = json.loads((lib.RUNS / "em_survivors_2026-10-06" / "em_survivors.json").read_text("utf-8"))
    surv15 = set(pd.read_csv(lib.RUNS / "em_survivors_2026-10-06" / "survivors.csv")["SNP_ID"])
    reg = [None]
    from collections import OrderedDict
    forest_mean = pd.read_csv(lib.RUNS / "snplib_dock_2026-10-06" / "pool69_predictions.csv").set_index("SNP_ID").loc[uid]["mean_pred_pki"].to_numpy()
    tob_mean = np.vstack(lib.expand(by)).mean(axis=0)
    out.update({"survivors": rk, "shortlist": sl, "shortlist_in_the_eight": len(set(sl) & set(prev["shortlist"])),
                "survivors_in_the_fifteen": len(set(rk) & surv15), "pool_spearman_vs_forest_mean": lib.sp(forest_mean, tob_mean)})
    (lib.OUT / "posthoc_lambda.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    for k, v in out["per_split"].items():
        print(k[:32], "lambda", v["lambda"], "sp %.4f" % v["heldout_point_spearman"], "up-pass %.3f" % v["upper_pass_frac"], "conv", v["converged"])
    print("survivors", len(rk), "in15", out["survivors_in_the_fifteen"], "shortlist in8", out["shortlist_in_the_eight"], "pool sp %.3f" % out["pool_spearman_vs_forest_mean"])
