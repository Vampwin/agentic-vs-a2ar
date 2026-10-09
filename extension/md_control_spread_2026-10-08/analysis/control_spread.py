"""Descriptive, post hoc: how much do the three 50 ns control runs differ, and how much would
a duplicate (n = 2) report depend on which two are used? Reads the series already written
by the unchanged ADR-AIIM-324/325 pipeline; re-derives nothing from trajectories.
Persistence = fraction of frames with the N6.55 minimum distance <= 3.5 A (0.35 nm)."""
import itertools, json, hashlib
from pathlib import Path
import numpy as np

R = Path(r"D:\Agentic VS Workflow\20_new_runs")
SRC = {"r1": R / "md_control_anchor_2026-10-07/analysis_prod",
       "r2": R / "md_control_replicates_2026-10-06/analysis_prod_r2",
       "r3": R / "md_control_replicates_2026-10-06/analysis_prod_r3"}
OUT = R / "md_control_spread_2026-10-08"


def xvg(p):
    a = np.array([[float(x) for x in l.split()[:2]] for l in p.read_text().splitlines()
                  if l.strip() and l[0] not in "#@"])
    return a[a[:, 0] <= 50.0 + 1e-9]


res, inputs = {}, {}
for k, d in SRC.items():
    n, r = d / "N6_55.xvg", d / "ligand_rmsd.xvg"
    inputs[str(n)] = hashlib.sha256(n.read_bytes()).hexdigest()
    inputs[str(r)] = hashlib.sha256(r.read_bytes()).hexdigest()
    N, Rm = xvg(n), xvg(r)
    ok = N[:, 1] * 10 <= 3.5
    blocks = [float(ok[(N[:, 0] > a) & (N[:, 0] <= a + 10)].mean()) for a in range(0, 50, 10)]
    lost = N[~ok, 0]
    # longest run of consecutive frames without the contact
    run = best = 0
    for v in ok:
        run = 0 if v else run + 1
        best = max(best, run)
    res[k] = {"frames": int(len(N)), "persistence": round(float(ok.mean()), 4),
              "persistence_by_10ns_block": [round(b, 3) for b in blocks],
              "first_frame_without_contact_ns": round(float(lost[0]), 1) if len(lost) else None,
              "longest_gap_ns": round(best * 0.1, 1),
              "rmsd_mean_A": round(float(Rm[:, 1].mean() * 10), 2),
              "rmsd_last_10ns_mean_A": round(float(Rm[Rm[:, 0] > 40, 1].mean() * 10), 2)}

p = np.array([res[k]["persistence"] for k in res]); m = np.array([res[k]["rmsd_mean_A"] for k in res])
spread = {"persistence": {"mean": round(float(p.mean()), 3), "sd": round(float(p.std(ddof=1)), 3),
                          "range": [float(p.min()), float(p.max())]},
          "rmsd_mean_A": {"mean": round(float(m.mean()), 2), "sd": round(float(m.std(ddof=1)), 2),
                          "range": [float(m.min()), float(m.max())]}}
pairs = {}
for a, b in itertools.combinations(res, 2):
    pairs[f"{a}+{b}"] = {"persistence_mean": round((res[a]["persistence"] + res[b]["persistence"]) / 2, 3),
                         "rmsd_mean_A": round((res[a]["rmsd_mean_A"] + res[b]["rmsd_mean_A"]) / 2, 2),
                         "persistence_difference": round(abs(res[a]["persistence"] - res[b]["persistence"]), 3)}
out = {"note": "descriptive and post hoc; no prediction graded", "runs": res, "spread_n3": spread,
       "duplicate_pairs": pairs, "inputs_sha256": inputs}
(OUT / "control_spread.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
print(json.dumps({k: v for k, v in out.items() if k != "inputs_sha256"}, indent=1))
