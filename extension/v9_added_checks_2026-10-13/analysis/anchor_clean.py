"""ADR-AIIM-427 anchor: every clean workspace under the frozen evaluators, with the added checks (direct) and deployed with the added checks.
Run before the pre-registration is written; no faulted workspace is built here."""
import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import harness as H


def main():
    H.write_plus_contracts()
    snap = H.frozen_tree_hashes()
    rows = []
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        for name, rec, cut in H.ext.G2_RECEPTORS:
            ws = td / f"g2_{name}"
            inp = H.build_g2_plus(ws, name, rec, cut)
            rows.append({"gate": "G2", "case": name, "frozen": H.run_frozen("G2", inp), "plus": H.run_plus("G2", inp),
                         "deployed_plus": H.run_deployed_plus("G2", inp, ws), "n_his": len(inp["protonation_decision"])})
        for cname, sid, reps in H.ext.clean_g7_cases():
            ws = td / f"g7_{cname}"
            inp = H.build_g7_plus(ws, sid, reps)
            last = {}
            df = H.pd.read_csv(inp["equilibration_check_csv"])
            for (c, r), g in df.groupby(["compound_id", "replicate"]):
                last[f"{c}/{r}"] = float(g["time_ns"].max())
            rows.append({"gate": "G7", "case": cname, "frozen": H.run_frozen("G7", inp), "plus": H.run_plus("G7", inp),
                         "deployed_plus": H.run_deployed_plus("G7", inp, ws), "last_time_ns": last})
        unrelated = [f"UNRELATED_{i:04d}" for i in range(40)]
        for name, *_ in H.ext2.G5_VARIANTS:
            for tag, ctx in (("", None), ("+context", unrelated if name in ("vina_only", "combined_zscore", "ifp_zm_only", "z_vina", "z_ifp") else None)):
                if tag and ctx is None:
                    continue
                ws = td / f"g5_{name}{tag}"
                g5, reg = H.build_g5_plus(ws, name, context_ids=ctx)
                fr = H.run_frozen_g5(g5, reg)
                g5, reg = H.build_g5_plus(td / f"g5p_{name}{tag}", name, context_ids=ctx)
                pl = H.run_plus_g5(g5, reg)
                g5d, _ = H.build_g5_plus(td / f"g5d_{name}{tag}", name, context_ids=ctx)
                rows.append({"gate": "G5", "case": name + tag, "frozen": fr, "plus": pl, "deployed_plus": H.run_deployed_plus_g5(td / f"g5d_{name}{tag}", g5d)})
    after = H.frozen_tree_hashes()
    out = {"adr": "ADR-AIIM-427 anchor (clean cases only)", "n_clean": len(rows), "frozen_tree_unchanged": snap == after, "n_frozen_files": len(snap), "rows": rows}
    (H.OUT / "anchor_clean.json").write_text(json.dumps(out, indent=1, default=str), encoding="utf-8")

    def state(r, k):
        v = r[k]
        return (v["G5a"]["state"] if v.get("G5a") else None) if r["gate"] == "G5" else v["state"]

    bad = [r["case"] for r in rows if not (state(r, "frozen") == "PASS" and state(r, "plus") == "PASS" and state(r, "deployed_plus") == "PASS")]
    print("clean cases:", len(rows), "| not PASS in all three:", bad, "| frozen tree unchanged:", snap == after, f"({len(snap)} files)")
    for r in rows:
        if r["gate"] == "G7":
            print("  ", r["case"], "min last time", round(min(r["last_time_ns"].values()), 2))


if __name__ == "__main__":
    main()
