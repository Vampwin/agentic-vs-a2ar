"""ADR-AIIM-416 anchor: every clean workspace must pass under both paths before any faulted case is run or read."""
import json, sys, tempfile, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import ext

def main():
    rows = []
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        for name, rec, cut in ext.G2_RECEPTORS:
            ws = td / f"g2_{name}"; inp = ext.build_g2(ws, name, rec, cut)
            rows.append({"gate": "G2", "case": name, "direct": ext.run_direct("G2", inp), "deployed": ext.run_deployed("G2", inp, ws)})
        for cname, sid, reps in ext.clean_g7_cases():
            ws = td / f"g7_{cname}"; inp = ext.build_g7(ws, sid, reps)
            rows.append({"gate": "G7", "case": cname, "direct": ext.run_direct("G7", inp), "deployed": ext.run_deployed("G7", inp, ws)})
    ok = all(r["direct"]["state"] == "PASS" and r["deployed"]["state"] == "PASS" and r["deployed"]["contained"] is False for r in rows)
    out = {"adr": "ADR-AIIM-416 anchor (clean cases only)", "n_clean": len(rows), "all_pass_both_paths": ok, "rows": rows}
    (ext.OUT / "anchor_clean.json").write_text(json.dumps(out, indent=1, default=str), encoding="utf-8")
    print("clean cases:", len(rows), "| all PASS under both paths and not contained:", ok)

if __name__ == "__main__":
    main()
