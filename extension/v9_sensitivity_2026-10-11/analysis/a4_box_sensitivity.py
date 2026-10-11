"""A4 (ADR-AIIM-399, amendment 1): the four box-constrained expansion compounds that completed the full protocol, redocked with
the frozen g4_dock_library.py in the production box (22 A, the anchor) and in a 30 A box. Nothing but the box size differs.

Usage:  python a4_box_sensitivity.py launch     -> writes the inputs, starts the two WSL processes and waits for them (stop rule: 120 min)
        python a4_box_sensitivity.py collect    -> reads the two score files and writes a4_box_sensitivity.json
"""
import json
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

ROOT = Path(r"D:\Agentic VS Workflow")
OUT = ROOT / "20_new_runs" / "v9_sensitivity_2026-10-11" / "a4_box"
FUNNEL = ROOT / "20_new_runs" / "np_library_funnel_2026-08-15"
G2 = ROOT / "20_new_runs" / "snp_screen_2026-08-04" / "G2_preparation"
SCRIPT = ROOT / "10_new_system" / "scripts" / "g4_dock_library.py"
IDS = ["NPLIB01997", "NPLIB03133", "NPLIB04561", "NPLIB03343"]
WSL_DISTRO, CPU = "Ubuntu-22.04", 6


def wslpath(p: Path) -> str:
    s = str(p).replace("\\", "/")
    return "/mnt/" + s[0].lower() + s[2:]


def launch():
    OUT.mkdir(parents=True, exist_ok=False)          # refuses a non-empty earlier attempt
    led = pd.read_csv(FUNNEL / "np_library_dock_ledger.csv")
    sub = led[led["compound_id"].isin(IDS)].set_index("compound_id").loc[IDS].reset_index()
    assert len(sub) == 4
    sub.to_csv(OUT / "ledger_four.csv", index=False)
    box = json.loads((G2 / "docking_box.json").read_text(encoding="utf-8"))
    assert box["size"] == [22, 22, 22]
    box30 = dict(box); box30["size"] = [30, 30, 30]
    (OUT / "docking_box_30.json").write_text(json.dumps(box30, indent=2), encoding="utf-8")
    runs = {"box22": G2 / "docking_box.json", "box30": OUT / "docking_box_30.json"}
    procs = []
    meta = {"started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "ids": IDS, "cpu": CPU, "pids": {}}
    for tag, bx in runs.items():
        cmd = (f"source ~/miniconda3/etc/profile.d/conda.sh && conda activate diffdock && "
               f"python '{wslpath(SCRIPT)}' --ledger '{wslpath(OUT / 'ledger_four.csv')}' "
               f"--receptor '{wslpath(G2 / '4EIY_prepared.pdbqt')}' --box '{wslpath(bx)}' "
               f"--out '{wslpath(OUT / f'scores_{tag}.csv')}' --poses-dir '{wslpath(OUT / f'poses_{tag}')}' "
               f"--id-col compound_id --smiles-col SMILES --n-conformers 8 --exhaustiveness 32 --n-poses 10 --cpu {CPU}")
        (OUT / f"poses_{tag}").mkdir()
        fh = open(OUT / f"dock_{tag}.log", "w", encoding="utf-8")
        p = subprocess.Popen(["wsl", "-d", WSL_DISTRO, "--", "bash", "-lc", cmd], stdout=fh, stderr=subprocess.STDOUT,
                             creationflags=getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0))
        meta["pids"][tag] = p.pid
        procs.append(p)
        (OUT / f"cmd_{tag}.txt").write_text(cmd, encoding="utf-8")
    (OUT / "launch.json").write_text(json.dumps(meta, indent=1), encoding="utf-8")
    print("launched", meta, flush=True)
    # stop rule of amendment 1: 120 minutes of wall time for the two processes together
    t0, STOP = time.time(), 120 * 60
    while any(p.poll() is None for p in procs) and time.time() - t0 < STOP:
        time.sleep(20)
    stopped = [p.pid for p in procs if p.poll() is None]
    for p in procs:
        if p.poll() is None:
            p.kill()
    meta["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    meta["wall_seconds"] = round(time.time() - t0, 1)
    meta["stopped_by_rule"] = stopped
    (OUT / "launch.json").write_text(json.dumps(meta, indent=1), encoding="utf-8")
    print("finished", meta, flush=True)


def collect():
    stored = pd.concat([pd.read_csv(f) for f in sorted(FUNNEL.glob("full_scores_chunk*.csv"))])
    stored = stored[stored["compound_id"].isin(IDS)].set_index("compound_id")
    res = {"adr": "ADR-AIIM-399 (amendment 1)", "ids": IDS, "per_compound": {}}
    tabs = {}
    for tag in ("box22", "box30"):
        f = OUT / f"scores_{tag}.csv"
        tabs[tag] = pd.read_csv(f).set_index("compound_id") if f.exists() else None
    done = all(t is not None and set(IDS) <= set(t.index) for t in tabs.values())
    res["completed"] = bool(done)
    if done:
        for c in IDS:
            s, a, b = stored.loc[c], tabs["box22"].loc[c], tabs["box30"].loc[c]
            res["per_compound"][c] = {
                "stored_score": float(s["vina_score"]), "stored_dist": float(s["centroid_dist_A"]),
                "box22_status": a["status"], "box22_score": float(a["vina_score"]), "box22_dist": float(a["centroid_dist_A"]),
                "box30_status": b["status"], "box30_score": float(b["vina_score"]), "box30_dist": float(b["centroid_dist_A"]),
                "box30_on_pocket": bool(b["on_pocket"]), "elapsed_s_box22": float(a["elapsed_s"]), "elapsed_s_box30": float(b["elapsed_s"]),
                "delta_30_minus_22": round(float(b["vina_score"]) - float(a["vina_score"]), 3)}
        pc = res["per_compound"]
        res["P8_anchor_22_reproduces_stored"] = all(abs(v["box22_score"] - v["stored_score"]) <= 0.001 for v in pc.values())
        res["P9_at_least_two_better_by_0.3"] = sum(v["delta_30_minus_22"] < -0.3 for v in pc.values()) >= 2
        res["P10_all_on_pocket_in_30"] = all(v["box30_on_pocket"] for v in pc.values())
        res["n_better_by_0.3"] = int(sum(v["delta_30_minus_22"] < -0.3 for v in pc.values()))
        res["n_changed_by_more_than_0.1"] = int(sum(abs(v["delta_30_minus_22"]) > 0.1 for v in pc.values()))
    (OUT.parent / "a4_box_sensitivity.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    {"launch": launch, "collect": collect}[sys.argv[1]]()
