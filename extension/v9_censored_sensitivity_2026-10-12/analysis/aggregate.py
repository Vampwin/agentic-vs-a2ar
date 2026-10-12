"""ADR-AIIM-412: read the job results, check the anchor, grade the sixteen predictions as pre-registered, write results.json."""
from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
sys.path.insert(0, str(Path(__file__).parent))
import numpy as np
import pandas as pd

import lib

RES = lib.OUT / "results"
ANCHOR_CONV = "D1sorted_D2trainfirst_D3divert"


def ld(name):
    z = np.load(RES / f"{name}.npz")
    m = json.loads((RES / f"{name}.json").read_text("utf-8"))
    return {k: z[k] for k in z.files}, m


def main():
    d = lib.data()
    uid, murcko = list(d["uid"]), list(d["murcko"])
    gs = lib.groups()
    em = lib.em_module()
    stored = pd.read_csv(lib.RUNS / "snplib_dock_2026-10-06" / "pool69_predictions.csv").set_index("SNP_ID").loc[uid]
    prev = json.loads((lib.RUNS / "em_survivors_2026-10-06" / "em_survivors.json").read_text("utf-8"))
    surv_prev = pd.read_csv(lib.RUNS / "em_survivors_2026-10-06" / "survivors.csv")["SNP_ID"].tolist()
    out = {"adr": "ADR-AIIM-412 (pre-registration) / ADR-AIIM-413 (results)",
           "prereg_sha256": "e161e8592971386bc04207a705dc5fd9f9094f0079f116ff798e0e498dca7420"}

    # ---- C-A anchor: pass 10 of the 30-pass run equals the registered one-thread result
    ca = {}
    reg_by_group = {}
    maxdiff = 0.0
    for gi, g in enumerate(gs):
        a, m = ld(f"c2_g{gi}_bound_30")
        p10 = a["pass_pred"][9]
        reg_by_group[g[0]] = p10
        for nm in g:
            col = stored[f"pred_{nm}"].to_numpy()
            ok = ~np.isnan(col)
            maxdiff = max(maxdiff, float(np.max(np.abs(col[ok] - p10[ok]))))
    ranked, shortlist = lib.select(lib.expand(reg_by_group), uid, murcko)
    ca = {"max_abs_diff_to_stored_4dp": maxdiff, "stored_matches": bool(maxdiff < 5e-5),
          "survivors_reproduced_in_order": ranked == surv_prev, "n_survivors": len(ranked),
          "shortlist_reproduced_in_order": shortlist == prev["shortlist"], "shortlist": shortlist}
    ca["PASS"] = bool(ca["stored_matches"] and ca["survivors_reproduced_in_order"] and ca["shortlist_reproduced_in_order"])
    out["C_A_anchor"] = ca

    # ---- C2 convergence
    c2 = {"per_split": {}}
    P5 = P6 = P7 = 0
    for gi, g in enumerate(gs):
        a, m = ld(f"c2_g{gi}_bound_30")
        mc = np.array(m["max_change"])
        first = next((i + 1 for i, v in enumerate(mc) if v < 0.01), None)
        p10, p30 = a["pass_pred"][9], a["pass_pred"][29]
        t10, t30 = set(np.argsort(-p10)[:20]), set(np.argsort(-p30)[:20])
        mean_c = np.array(m["mean_change"])
        mono = bool(np.all(np.diff(mean_c[1:]) <= 1e-12))
        rec = {"group": g, "first_pass_max_change_below_0.01": first, "max_change_pass10": float(mc[9]),
               "max_change_pass30": float(mc[29]), "mean_change_pass10": float(mean_c[9]), "mean_change_pass30": float(mean_c[29]),
               "n_moved_gt_0.01_pass10": m["n_moved"][9], "n_moved_gt_0.01_pass30": m["n_moved"][29],
               "spearman_pass10_vs_pass30": lib.sp(p10, p30), "top20_overlap_pass10_pass30": len(t10 & t30),
               "mean_change_non_increasing_from_pass2": mono,
               "max_change_by_pass": [round(float(x), 5) for x in mc]}
        c2["per_split"][g[0]] = rec
        P5 += first is not None
        P6 += (rec["spearman_pass10_vs_pass30"] >= 0.99) and (rec["top20_overlap_pass10_pass30"] >= 18)
        P7 += mono
    c2["grading"] = {"P5": {"splits_below_tol_within_30": P5, "need": 4, "verdict": "MET" if P5 >= 4 else "MISSED"},
                     "P6": {"splits_meeting": P6, "need": 5, "verdict": "MET" if P6 == 5 else "MISSED"},
                     "P7": {"splits_monotone": P7, "need": 5, "verdict": "MET" if P7 == 5 else "MISSED"}}
    out["C2_convergence"] = c2

    # ---- C3 initialisation
    c3 = {"inits": {}}
    sp_min, jac_min, sl_min = 1.0, 1.0, 8
    base_surv = set(ranked)
    base_sl = set(shortlist)
    sets_prev15, sets_prev8 = set(surv_prev), set(prev["shortlist"])
    for init in ("inside1", "random", "mean"):
        by_group = {}
        per = {}
        for gi, g in enumerate(gs):
            a, m = ld(f"c3_g{gi}_{init}")
            by_group[g[0]] = a["final_pred"]
            ref = reg_by_group[g[0]]
            s = lib.sp(a["final_pred"], ref)
            ov = len(set(np.argsort(-a["final_pred"])[:20]) & set(np.argsort(-ref)[:20]))
            per[g[0]] = {"spearman_vs_start_at_bound": s, "top20_overlap": ov, "n_pass": m["n_pass"],
                         "stopped_by_tol": m["stopped_by_tol"], "final_change": m["final_change"]}
            sp_min = min(sp_min, s)
        rk, sl = lib.select(lib.expand(by_group), uid, murcko)
        jac15 = len(set(rk) & sets_prev15) / len(set(rk) | sets_prev15)
        jac8 = len(set(sl) & sets_prev8) / len(set(sl) | sets_prev8)
        kept8 = len(set(sl) & sets_prev8)
        jac_min, sl_min = min(jac_min, jac15), min(sl_min, kept8)
        c3["inits"][init] = {"per_split": per, "n_survivors": len(rk), "survivors": rk, "shortlist": sl,
                             "jaccard_vs_fifteen": jac15, "jaccard_vs_eight": jac8, "eight_kept_in_shortlist": kept8}
    c3["grading"] = {"P8": {"min_spearman": sp_min, "need": 0.98, "verdict": "MET" if sp_min >= 0.98 else "MISSED"},
                     "P9": {"min_jaccard_survivors": jac_min, "need": 0.8, "verdict": "MET" if jac_min >= 0.8 else "MISSED"},
                     "P10": {"min_eight_kept": sl_min, "need": 6, "verdict": "MET" if sl_min >= 6 else "MISSED"}}
    out["C3_initialisation"] = c3

    # ---- C4 ridge-Tobit
    reg_em = json.loads((lib.RUNS / "pki_censored_em_2026-10-05" / "pki_censored_em.json").read_text("utf-8"))["conventions"]
    tob_by_group, c4 = {}, {"per_split": {}, "per_convention": {}}
    lower_all, passes_all = True, True
    for gi, g in enumerate(gs):
        a, m = ld(f"c4_g{gi}")
        tob_by_group[g[0]] = a["pool_pred"]
        c4["per_split"][g[0]] = {k: m[k] for k in ("lambda", "sigma", "converged", "n_iter", "n_features", "heldout_point_spearman",
                                                  "heldout_point_r2", "heldout_point_rmse", "upper_n", "upper_pass_frac",
                                                  "lower_n", "lower_pass_frac", "pool_min", "pool_max")}
        for nm in g:
            r = reg_em[nm]
            tob_sp, rf_sp = m["heldout_point_spearman"], r["spearman"]
            tob_up, rf_up = m["upper_pass_frac"], r["censored_heldout_test"]["upper_pass_frac"]
            c4["per_convention"][nm] = {"tobit_spearman": tob_sp, "forest_spearman_registered": rf_sp,
                                        "tobit_upper_pass": tob_up, "forest_upper_pass_registered": rf_up}
            lower_all &= tob_sp < rf_sp
            passes_all &= tob_up > rf_up
    rk_t, sl_t = lib.select(lib.expand(tob_by_group), uid, murcko)
    forest_mean = np.vstack(lib.expand(reg_by_group)).mean(axis=0)
    tob_mean = np.vstack(lib.expand(tob_by_group)).mean(axis=0)
    sp_pool = lib.sp(forest_mean, tob_mean)
    in8 = len(set(sl_t) & sets_prev8)
    c4["tobit_survivors"], c4["tobit_shortlist"] = rk_t, sl_t
    c4["pool_spearman_forest_vs_tobit_mean"] = sp_pool
    c4["tobit_shortlist_in_the_eight"] = in8
    c4["tobit_survivors_in_the_fifteen"] = len(set(rk_t) & sets_prev15)
    c4["grading"] = {"P11": {"tobit_lower_in_all_8": bool(lower_all), "verdict": "MET" if lower_all else "MISSED"},
                     "P12": {"spearman": sp_pool, "need": 0.5, "verdict": "MET" if sp_pool >= 0.5 else "MISSED"},
                     "P13": {"tobit_shortlist_in_eight": in8, "need_at_most": 4, "verdict": "MET" if in8 <= 4 else "MISSED"},
                     "P14": {"tobit_higher_in_all_8": bool(passes_all), "verdict": "MET" if passes_all else "MISSED"}}
    out["C4_ridge_tobit"] = c4

    # ---- C5 determinism
    a1, _ = ld("c5_a")
    a2, _ = ld("c5_b")
    a3, _ = ld("c5_c")
    same = hashlib.sha256(a1["final_pred"].tobytes()).hexdigest() == hashlib.sha256(a2["final_pred"].tobytes()).hexdigest()
    dif3 = float(np.max(np.abs(a3["final_pred"] - a1["final_pred"])))
    out["C5_determinism"] = {"one_thread_runs_byte_identical": bool(same), "max_abs_diff_three_threads": dif3,
                             "n_pool_predictions_differing_three_threads": int((a3["final_pred"] != a1["final_pred"]).sum()),
                             "grading": {"P15": {"verdict": "MET" if same else "MISSED"},
                                         "P16": {"verdict": "MET" if dif3 > 0 else "MISSED"}}}

    # ---- C1 simulation
    c1 = {"per_run": {}, "summary": {}}
    methods = ("drop", "bound", "em", "tobit")
    tab = {}
    for truth in ("lin", "rf"):
        for seed in range(1, 6):
            _, m = ld(f"c1_{truth}_{seed}")
            c1["per_run"][f"{truth}_{seed}"] = {k: m[k] for k in methods} | {"lo_bound": m["lo_bound"], "n_te_below": m["n_te_below"],
                                                                              "tobit_lambda": m["tobit_lambda"]}
            tab[(truth, seed)] = m
    def wins(truth, a, b, key="rmse"):
        return sum(tab[(truth, s)][a][key] < tab[(truth, s)][b][key] for s in range(1, 6))
    P1 = {t: wins(t, "em", "bound") for t in ("lin", "rf")}
    P2 = {t: wins(t, "em", "drop") for t in ("lin", "rf")}
    P3a, P3b = wins("lin", "tobit", "em"), wins("rf", "em", "tobit")
    P4 = {}
    for t in ("lin", "rf"):
        P4[t] = sum(tab[(t, s)]["em"]["spearman"] >= max(tab[(t, s)][k]["spearman"] for k in methods) - 0.03 for s in range(1, 6))
    mean_rmse = {t: {k: float(np.mean([tab[(t, s)][k]["rmse"] for s in range(1, 6)])) for k in methods} for t in ("lin", "rf")}
    mean_sp = {t: {k: float(np.mean([tab[(t, s)][k]["spearman"] for s in range(1, 6)])) for k in methods} for t in ("lin", "rf")}
    reach = {t: {k: float(np.mean([tab[(t, s)][k]["share_below_bound_predicted_below"] for s in range(1, 6)])) for k in methods} for t in ("lin", "rf")}
    c1["summary"] = {"mean_rmse": mean_rmse, "mean_spearman": mean_sp, "mean_share_below_bound_predicted_below": reach}
    c1["grading"] = {
        "P1": {"em_beats_bound_replicates": P1, "need_each": 4, "verdict": "MET" if all(v >= 4 for v in P1.values()) else "MISSED"},
        "P2": {"em_beats_drop_replicates": P2, "need_each": 4, "verdict": "MET" if all(v >= 4 for v in P2.values()) else "MISSED"},
        "P3": {"tobit_beats_em_under_lin": P3a, "em_beats_tobit_under_rf": P3b, "need": 4,
               "verdict": "MET" if (P3a >= 4 and P3b >= 4) else "MISSED"},
        "P4": {"em_spearman_within_0.03_of_best": P4, "need_each": 4, "verdict": "MET" if all(v >= 4 for v in P4.values()) else "MISSED"}}
    out["C1_simulation"] = c1

    grades = {}
    for sec in ("C1_simulation", "C2_convergence", "C3_initialisation", "C4_ridge_tobit"):
        for k, v in out[sec]["grading"].items():
            grades[k] = v["verdict"]
    for k, v in out["C5_determinism"]["grading"].items():
        grades[k] = v["verdict"]
    out["grades"] = dict(sorted(grades.items(), key=lambda kv: int(kv[0][1:])))
    out["n_met"] = sum(v == "MET" for v in grades.values())
    out["n_predictions"] = len(grades)
    (lib.OUT / "results.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("anchor:", ca["PASS"], ca["max_abs_diff_to_stored_4dp"])
    print("grades:", out["grades"], f"-> {out['n_met']} of {out['n_predictions']} met")


if __name__ == "__main__":
    main()
