"""ADR-AIIM-417: every injected instance must differ from its clean base in the file the fault targets (no vacuous injection), and the
file the fault does not target must be byte-identical to the clean one."""
import hashlib, json, sys, tempfile
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import ext

def h(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    rows = []
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        n, rc, cut = ext.G2_RECEPTORS[0]
        clean = td / "g2_clean"; ext.build_g2(clean, n, rc, cut)
        for i, p in enumerate(ext.F06_GRID):
            ws = td / f"f06_{i}"; ext.build_g2(ws, n, rc, cut); ext.inject_f06(ws, p, 7)
            same_other = all(h(ws / f) == h(clean / f) for f in ("ZMA_reference_4EIY.sdf", "docking_box.json"))
            rows.append({"family": "F06", "instance": i, "target_file_changed": h(ws / "4EIY_prepared.pdbqt") != h(clean / "4EIY_prepared.pdbqt"), "other_inputs_identical": same_other})
        clean7 = td / "g7_clean"; ext.build_g7(clean7, "SNP00157", ["r1", "r2", "r3"])
        for fam, grid in (("F15", ext.F15_GRID), ("F15L", ext.F15_GRID), ("F16", ext.F16_GRID)):
            for i, p in enumerate(grid):
                ws = td / f"{fam}_{i}"; ext.build_g7(ws, "SNP00157", ["r1", "r2", "r3"])
                if fam == "F16": ext.inject_f16(ws, p)
                else: ext.inject_f15(ws, p, log_updated=(fam == "F15L"))
                tgt = {"F15": "equilibration_check.csv", "F15L": "equilibration_check.csv", "F16": "md_results.csv"}[fam]
                others = [f for f in ("md_completion_log.csv", "equilibration_check.csv", "md_results.csv", "selected_candidates.csv") if f != tgt]
                if fam == "F15L": others = [f for f in others if f != "md_completion_log.csv"]
                rows.append({"family": fam, "instance": i, "target_file_changed": h(ws / tgt) != h(clean7 / tgt),
                             "other_inputs_identical": all(h(ws / f) == h(clean7 / f) for f in others),
                             "log_changed": h(ws / "md_completion_log.csv") != h(clean7 / "md_completion_log.csv")})
    out = {"n": len(rows), "all_target_changed": all(r["target_file_changed"] for r in rows),
           "all_other_identical": all(r["other_inputs_identical"] for r in rows),
           "F15_log_unchanged_F15L_log_changed": all((not r["log_changed"]) for r in rows if r["family"] in ("F15", "F16")) and all(r["log_changed"] for r in rows if r["family"] == "F15L"),
           "rows": rows}
    (ext.OUT / "verify_injections.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print({k: v for k, v in out.items() if k != "rows"})

if __name__ == "__main__":
    main()
