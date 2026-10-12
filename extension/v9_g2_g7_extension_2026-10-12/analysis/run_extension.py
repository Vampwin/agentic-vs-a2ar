"""ADR-AIIM-416 (pre-registration) / ADR-AIIM-417 (results): evaluate the faulted workspaces under both paths, grade P1 to P7."""
from __future__ import annotations

import json
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import ext

SEED = 7


def both(gate, inputs, ws):
    return {"direct": ext.run_direct(gate, inputs), "deployed": ext.run_deployed(gate, inputs, ws)}


def main():
    t0 = time.time()
    rows = []
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        base_g2 = ext.G2_RECEPTORS[0]            # the registered run's receptor is the base of the F06 instances
        for i, p in enumerate(ext.F06_GRID):
            ws = td / f"F06_{i:02d}"
            inp = ext.build_g2(ws, base_g2[0], base_g2[1], base_g2[2])
            rec = ext.inject_f06(ws, p, SEED)
            rows.append({"family": "F06", "gate": "G2", "instance": i, "params": p, "injection": rec, **both("G2", inp, ws)})
        base_g7 = ("SNP00157", ["r1", "r2", "r3"])  # the compound with three replicates simulated in this project
        for fam, grid, upd in (("F15", ext.F15_GRID, False), ("F15L", ext.F15_GRID, True), ("F16", ext.F16_GRID, None)):
            for i, p in enumerate(grid):
                ws = td / f"{fam}_{i:02d}"
                inp = ext.build_g7(ws, *base_g7)
                rec = ext.inject_f15(ws, p, log_updated=upd) if fam.startswith("F15") else ext.inject_f16(ws, p)
                rows.append({"family": fam, "gate": "G7", "instance": i, "params": p, "injection": rec, **both("G7", inp, ws)})

    anchor = json.loads((ext.OUT / "anchor_clean.json").read_text("utf-8"))
    out = {"adr": "ADR-AIIM-416 (pre-registration) / ADR-AIIM-417 (results)",
           "prereg_sha256": "90e0820e14bf1c5ab82817f38aafc296916b555f13b96aaa2632bff77638647b",
           "n_clean": anchor["n_clean"], "clean_all_pass_both_paths": anchor["all_pass_both_paths"], "cases": rows}

    def fam(f):
        rs = [r for r in rows if r["family"] == f and not r["injection"].get("vacuous")]
        n = len(rs)
        det_b = sum(r["direct"]["state"] != "PASS" for r in rs)
        det_c = sum(r["deployed"]["state"] != "PASS" for r in rs)
        cont = sum(bool(r["deployed"].get("contained")) for r in rs)
        fired = sorted({c for r in rs for c in r["direct"].get("checks_failed", [])})
        return {"n": n, "vacuous_excluded": sum(bool(r["injection"].get("vacuous")) for r in rows if r["family"] == f),
                "detected_direct": det_b, "detected_deployed": det_c, "contained_deployed": cont, "checks_that_fired_direct": fired}

    out["by_family"] = {f: fam(f) for f in ("F06", "F15", "F15L", "F16")}
    bf = out["by_family"]
    state_eq = all(r["direct"]["state"] == r["deployed"]["state"] for r in rows) and all(
        c["direct"]["state"] == c["deployed"]["state"] for c in anchor["rows"])
    cont_eq = all(bool(r["deployed"].get("contained")) == (r["direct"]["state"] != "PASS") for r in rows)
    out["grading"] = {
        "P1_F06_zero_detected": {"direct": bf["F06"]["detected_direct"], "deployed": bf["F06"]["detected_deployed"], "n": bf["F06"]["n"],
                                 "verdict": "MET" if bf["F06"]["detected_direct"] == 0 and bf["F06"]["detected_deployed"] == 0 else "MISSED"},
        "P2_F15_silent_zero_detected": {"direct": bf["F15"]["detected_direct"], "deployed": bf["F15"]["detected_deployed"], "n": bf["F15"]["n"],
                                        "verdict": "MET" if bf["F15"]["detected_direct"] == 0 and bf["F15"]["detected_deployed"] == 0 else "MISSED"},
        "P3_F15L_all_detected": {"direct": bf["F15L"]["detected_direct"], "deployed": bf["F15L"]["detected_deployed"], "n": bf["F15L"]["n"],
                                 "checks": bf["F15L"]["checks_that_fired_direct"],
                                 "verdict": "MET" if (bf["F15L"]["detected_direct"] == bf["F15L"]["n"] == bf["F15L"]["detected_deployed"]
                                                      and bf["F15L"]["checks_that_fired_direct"] == ["trajectory_completion"]) else "MISSED"},
        "P4_F16_zero_detected": {"direct": bf["F16"]["detected_direct"], "deployed": bf["F16"]["detected_deployed"], "n": bf["F16"]["n"],
                                 "verdict": "MET" if bf["F16"]["detected_direct"] == 0 and bf["F16"]["detected_deployed"] == 0 else "MISSED"},
        "P5_state_equal_both_paths": {"all_cases_equal": bool(state_eq), "verdict": "MET" if state_eq else "MISSED"},
        "P6_containment_equals_detection": {"all_cases": bool(cont_eq), "verdict": "MET" if cont_eq else "MISSED"},
        "P7_clean_anchor": {"n_clean": anchor["n_clean"], "all_pass": anchor["all_pass_both_paths"],
                            "verdict": "MET" if anchor["all_pass_both_paths"] else "MISSED"},
    }
    out["n_met"] = sum(g["verdict"] == "MET" for g in out["grading"].values())
    out["seconds"] = round(time.time() - t0, 1)
    (ext.OUT / "results.json").write_text(json.dumps(out, indent=1, default=str), encoding="utf-8")
    print(json.dumps(out["by_family"], indent=1))
    for k, g in out["grading"].items():
        print(f"{k:38s} {g['verdict']}")
    print("met", out["n_met"], "of 7;", out["seconds"], "s")


if __name__ == "__main__":
    main()
