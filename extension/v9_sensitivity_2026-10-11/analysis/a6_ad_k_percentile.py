"""A6 (ADR-AIIM-399): the applicability rule of the verification study at k = 3, 5, 10 and percentile = 5, 10, 20.

The frozen 10_new_system/scripts/ad_descriptor_sweep.py is imported unchanged; only its module constants K and PCTL are set
here before main() is called, and its --out is pointed at a folder of this run. Anchor: k = 5, percentile 10 reproduce the
stored admissions of ad_descriptor_sweep_2026-08-15. Stop rule: 40 minutes; settings completed by then are reported as such.
"""
import importlib.util
import json
import sys
import time
from pathlib import Path

import pandas as pd

ROOT = Path(r"D:\Agentic VS Workflow")
OUT = ROOT / "20_new_runs" / "v9_sensitivity_2026-10-11" / "a6_ad"
STORED = ROOT / "20_new_runs" / "ad_descriptor_sweep_2026-08-15" / "ad_coverage_by_descriptor.csv"
SETTINGS = [(5, 10)] + [(k, p) for k in (3, 5, 10) for p in (5, 10, 20) if (k, p) != (5, 10)]   # the anchor first
STOP = 40 * 60

sys.path.insert(0, str(ROOT / "10_new_system" / "scripts"))
spec = importlib.util.spec_from_file_location("_sweep", ROOT / "10_new_system" / "scripts" / "ad_descriptor_sweep.py")
sweep = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sweep)

OUT.mkdir(parents=True, exist_ok=False)
t0 = time.time()
res = {"adr": "ADR-AIIM-399", "settings": {}, "completed": [], "not_run": []}
stored = pd.read_csv(STORED).set_index("descriptor_space")
for k, p in SETTINGS:
    if time.time() - t0 > STOP:
        res["not_run"].append([k, p]); continue
    sweep.K, sweep.PCTL = k, p
    d = OUT / f"k{k}_p{p}"
    sys.argv = ["ad_descriptor_sweep.py", "--out", str(d)]
    sweep.main()
    tab = pd.read_csv(d / "ad_coverage_by_descriptor.csv").set_index("descriptor_space")
    res["settings"][f"k{k}_p{p}"] = {s: {"n_inside": int(r["n_inside"]), "n_total": int(r["n_total"]), "threshold": float(r["threshold"])}
                                    for s, r in tab.iterrows()}
    res["completed"].append([k, p])
    if (k, p) == (5, 10):
        res["anchor_reproduces_stored"] = bool((tab["n_inside"] == stored["n_inside"]).all() and (tab["threshold"] == stored["threshold"]).all())
        print("ANCHOR", res["anchor_reproduces_stored"], flush=True)
    print(f"done k={k} p={p} [{time.time()-t0:.0f}s]", flush=True)
    res["seconds"] = round(time.time() - t0, 1)
    (OUT.parent / "a6_ad_k_percentile.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
S = res["settings"]
four = ("ECFP4", "FCFP4", "TopTorsion", "Gobbi_Pharm2D")
res["P13_anchor"] = res.get("anchor_reproduces_stored")
res["P14_four_spaces_at_most_5_everywhere"] = all(S[c][s]["n_inside"] <= 5 for c in S for s in four)
res["P15_physchem_between_50_and_80pct_everywhere"] = all(0.50 <= S[c]["PhysChem"]["n_inside"] / S[c]["PhysChem"]["n_total"] <= 0.80 for c in S)
res["max_inside_four_spaces"] = max(S[c][s]["n_inside"] for c in S for s in four)
res["physchem_range"] = [min(S[c]["PhysChem"]["n_inside"] for c in S), max(S[c]["PhysChem"]["n_inside"] for c in S)]
res["seconds"] = round(time.time() - t0, 1)
(OUT.parent / "a6_ad_k_percentile.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
print(json.dumps({k: v for k, v in res.items() if k != "settings"}, indent=1))
