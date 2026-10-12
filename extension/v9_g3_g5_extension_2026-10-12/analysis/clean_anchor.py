"""ADR-AIIM-420 anchor: the clean workspaces of G3 and G5 under both paths, before any faulted case is run or read."""
import json, sys, tempfile
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import ext2

BASE_G3_FAILED = ["control_representativeness"]           # the production run's own outcome (ESCALATE)

def g3_clean_sources():
    r1 = ext2.RUNS / "snp_screen_repro_run1_2026-08-15" / "G3_redocking"
    r2 = ext2.RUNS / "snp_screen_repro_run2_2026-08-15" / "G3_redocking"
    src = [("production", ext2.G3D / "ZMA_redocked_4EIY.sdf", ext2.G3D / "atom_map_4eiy.json"),
           ("repro_run1", r1 / "ZMA_redocked_4EIY.sdf", r1 / "atom_map_4eiy.json"),
           ("repro_run2", r2 / "ZMA_redocked_4EIY.sdf", r2 / "atom_map_4eiy.json")]
    for k in range(1, 9):
        src.append((f"redock_seed{ext2.SEED + k}", ext2.OUT / "poses_clean" / f"redock_seed{ext2.SEED + k}.sdf", ext2.G3D / "atom_map_4eiy.json"))
    return src

def main():
    rows = []
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        for name, poses, amap in g3_clean_sources():
            ws = td / f"g3_{name}"
            inp, reg = ext2.build_g3(ws, poses, atom_map_src=amap)
            rows.append({"gate": "G3", "case": name, "direct": ext2.run_direct_g3(inp, reg), "deployed": ext2.run_deployed_g3(ws, inp)})
        for name, *_ in ext2.G5_VARIANTS:
            ws = td / f"g5_{name}"
            g5, reg = ext2.build_g5(ws, name)
            d = ext2.run_direct_g5(g5, reg)
            ws2 = td / f"g5d_{name}"
            g5b, _ = ext2.build_g5(ws2, name)
            rows.append({"gate": "G5", "case": name, "metrics": g5["metrics"], "direct": d, "deployed": ext2.run_deployed_g5(ws2, g5b)})
    ok_g3 = all(r["direct"]["checks_failed"] == BASE_G3_FAILED and r["deployed"]["checks_failed"] == BASE_G3_FAILED for r in rows if r["gate"] == "G3")
    ok_g5 = all(r["direct"]["G5a"]["state"] == "PASS" and r["deployed"]["G5a"]["state"] == "PASS" for r in rows if r["gate"] == "G5")
    out = {"adr": "ADR-AIIM-420 anchor (clean cases only)", "n_clean_g3": 11, "n_clean_g5": 11,
           "g3_failed_set_equals_production_all": ok_g3, "g5a_pass_all": ok_g5, "rows": rows}
    (ext2.OUT / "anchor_clean.json").write_text(json.dumps(out, indent=1, default=str), encoding="utf-8")
    print("G3: failed set == production's [control_representativeness] in all 11, both paths:", ok_g3)
    print("G5a PASS in all 11, both paths:", ok_g5)
    for r in rows:
        if r["gate"] == "G3": print(f"  {r['case']:24s} rmsd={r['direct']['rank1_rmsd']} state={r['direct']['state']} failed={r['direct']['checks_failed']}")

if __name__ == "__main__":
    main()
