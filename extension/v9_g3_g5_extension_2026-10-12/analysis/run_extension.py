"""ADR-AIIM-420 (pre-registration) / ADR-AIIM-421 (results): evaluate the faulted workspaces of G3 and G5 under both paths, grade P1 to P9."""
from __future__ import annotations

import json
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import ext2

BASE_G3 = ["control_representativeness"]
BASE_G5A = []
BASE_G5B = ["ef_1pct", "bedroc20"]            # the Vina ranking of the real benchmark


def new_g3(failed, state="OK"):
    return sorted(set(failed) - set(BASE_G3)) or ([] if state != "ERROR" else ["ERROR"])


def g3_case(family, i, params, ws, inputs, reg, extra):
    d = ext2.run_direct_g3(inputs, reg)
    c = ext2.run_deployed_g3(ws, inputs)
    nd = new_g3(d["checks_failed"], d["state"])
    nc = new_g3(c["checks_failed"], c["state"])
    return {"family": family, "gate": "G3", "instance": i, "params": params, "injection": extra, "direct": d, "deployed": c,
            "new_failed_direct": nd, "new_failed_deployed": nc, "detected_direct": bool(nd), "detected_deployed": bool(nc)}


def new_g5(res):
    new = []
    if res["G5a"]["state"] != "PASS":
        new += [f"G5a:{x}" for x in res["G5a"]["checks_failed"]] or ["G5a:ERROR"]
    elif res["G5b"] is not None:
        new += [f"G5b:{x}" for x in res["G5b"]["checks_failed"] if x not in BASE_G5B]
        if res["G5b"]["state"] == "ERROR":
            new.append("G5b:ERROR")
    return new


def g5_case(family, i, params, td, name, extra, **kw):
    ws = td / f"{family}_{i:02d}"
    g5, reg = ext2.build_g5(ws, "vina_only", **kw)
    d = ext2.run_direct_g5(g5, reg)
    ws2 = td / f"{family}_{i:02d}_dep"
    g5b, _ = ext2.build_g5(ws2, "vina_only", **kw)
    c = ext2.run_deployed_g5(ws2, g5b)
    nd, nc = new_g5(d), new_g5(c)
    return {"family": family, "gate": "G5", "instance": i, "params": params, "injection": extra, "metrics": g5["metrics"], "direct": d,
            "deployed": c, "new_failed_direct": nd, "new_failed_deployed": nc, "detected_direct": bool(nd), "detected_deployed": bool(nc)}


def main():
    t0 = time.time()
    rows = []
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        # ---- F07: the redocked poses of the shifted boxes
        for i, p in enumerate(ext2.F07_GRID):
            poses = ext2.OUT / "poses_f07" / f"F07_{i:02d}_{p['axis']}{p['offset_angstrom']}.sdf"
            ws = td / f"F07_{i:02d}"
            inputs, reg = ext2.build_g3(ws, poses)
            rows.append(g3_case("F07", i, p, ws, inputs, reg, {"box_center": ext2.f07_center(p), "poses": poses.name}))
        # ---- F08: the reference ligand rotated, the clean poses
        for i, p in enumerate(ext2.F08_GRID):
            ws = td / f"F08_{i:02d}"
            inputs, reg = ext2.build_g3(ws, ext2.G3D / "ZMA_redocked_4EIY.sdf")
            inj = ext2.rotate_reference(ext2.CRYSTAL, ws / "ZMA_crystal_4EIY.sdf", p["rotation_degrees"], p["rotation_axis"])
            rows.append(g3_case("F08", i, p, ws, inputs, reg, inj))
        # ---- F11 and F11L: known actives exposed to the ranking step
        for fam, from_bench in (("F11", False), ("F11L", True)):
            for i, p in enumerate(ext2.F11_GRID):
                ids = ext2.leaked_ids(p["n_leaked_ids"])
                rows.append(g5_case(fam, i, p, td, "vina_only", {"n_leaked": len(ids)}, leak_ids=ids, leak_from_benchmark=from_bench))
        # ---- F12: the labels swapped
        for i, p in enumerate(ext2.F12_GRID):
            lab = ext2.swapped_labels(p["swap_mode"])
            base = ext2.bench_frame()["label"].to_numpy()
            rows.append(g5_case("F12", i, p, td, "vina_only", {"n_labels_changed": int((lab != base).sum())}, label_swap={"labels": lab}))

    anchor = json.loads((ext2.OUT / "anchor_clean.json").read_text("utf-8"))
    out = {"adr": "ADR-AIIM-420 (pre-registration) / ADR-AIIM-421 (results)", "prereg_sha256": "7db40ba84d6818ecdbc3ba739d7de280e82a20f105b657b0b1e8e08868c9c48a",
           "anchor": {k: anchor[k] for k in ("n_clean_g3", "n_clean_g5", "g3_failed_set_equals_production_all", "g5a_pass_all")}, "cases": rows}

    def fam(f):
        rs = [r for r in rows if r["family"] == f]
        return {"n": len(rs), "detected_direct": sum(r["detected_direct"] for r in rs), "detected_deployed": sum(r["detected_deployed"] for r in rs),
                "new_failed_checks": sorted({c for r in rs for c in r["new_failed_direct"]})}
    out["by_family"] = {f: fam(f) for f in ("F07", "F08", "F11", "F11L", "F12")}

    def sel(f, cond):
        return [r for r in rows if r["family"] == f and cond(r["params"])]
    f08_big = sel("F08", lambda p: p["rotation_degrees"] >= 90)
    f07_far = sel("F07", lambda p: p["offset_angstrom"] >= 20)
    f07_near = sel("F07", lambda p: p["offset_angstrom"] == 5)
    f07_by_offset = {str(o): {"detected": sum(r["detected_direct"] for r in sel("F07", lambda p, o=o: p["offset_angstrom"] == o)),
                              "n": len(sel("F07", lambda p, o=o: p["offset_angstrom"] == o))} for o in (5, 10, 15, 20, 25)}
    f08_by_angle = {str(a): {"detected": sum(r["detected_direct"] for r in sel("F08", lambda p, a=a: p["rotation_degrees"] == a)),
                             "n": len(sel("F08", lambda p, a=a: p["rotation_degrees"] == a))} for a in (45, 90, 120, 135, 150)}
    both_all = lambda rs: all(r["detected_direct"] and r["detected_deployed"] for r in rs)
    none_any = lambda rs: not any(r["detected_direct"] or r["detected_deployed"] for r in rs)
    g3_eq = all(r["new_failed_direct"] == r["new_failed_deployed"] and sorted(r["direct"]["checks_failed"]) == sorted(r["deployed"]["checks_failed"])
                for r in rows if r["gate"] == "G3")
    g5_eq = all(r["direct"]["G5a"]["checks_failed"] == r["deployed"]["G5a"]["checks_failed"] and
                ((r["direct"]["G5b"] is None) == (r["deployed"]["G5b"] is None)) and
                (r["direct"]["G5b"] is None or r["direct"]["G5b"]["checks_failed"] == r["deployed"]["G5b"]["checks_failed"])
                for r in rows if r["gate"] == "G5")
    f11l = [r for r in rows if r["family"] == "F11L"]
    f12 = [r for r in rows if r["family"] == "F12"]
    f12_via_auc = all(r["direct"]["G5a"]["state"] == "PASS" and r["direct"]["G5b"] is not None and "auc_roc" in r["direct"]["G5b"]["checks_failed"] for r in f12)
    f11l_via_iso = all("G5a:label_isolation" in r["new_failed_direct"] for r in f11l)
    mm = ext2.metrics  # noqa: F841 (the metrics of the clean Vina ranking are in the anchor rows)
    anc_vina = next(r for r in anchor["rows"] if r["gate"] == "G5" and r["case"] == "vina_only")["metrics"]
    p1 = bool(anchor["g3_failed_set_equals_production_all"] and anchor["g5a_pass_all"] and
              all((r["direct"]["rank1_rmsd"] or 9) <= 2.0 for r in anchor["rows"] if r["gate"] == "G3") and
              (anc_vina["auc_roc"], anc_vina["bedroc20"], anc_vina["ef_1pct"]) == (0.6014, 0.2752, 1.6333))
    G = {
        "P1_anchor": {"verdict": "MET" if p1 else "MISSED"},
        "P2_F08_ge90_all_detected": {"detected_direct": sum(r["detected_direct"] for r in f08_big), "detected_deployed": sum(r["detected_deployed"] for r in f08_big),
                                     "n": len(f08_big), "verdict": "MET" if both_all(f08_big) else "MISSED"},
        "P3_F07_far_all_detected": {"detected_direct": sum(r["detected_direct"] for r in f07_far), "detected_deployed": sum(r["detected_deployed"] for r in f07_far),
                                    "n": len(f07_far), "verdict": "MET" if both_all(f07_far) else "MISSED"},
        "P4_F07_5A_none_detected": {"detected_direct": sum(r["detected_direct"] for r in f07_near), "detected_deployed": sum(r["detected_deployed"] for r in f07_near),
                                    "n": len(f07_near), "verdict": "MET" if none_any(f07_near) else "MISSED"},
        "P5_G3_deployed_equals_direct": {"verdict": "MET" if g3_eq else "MISSED"},
        "P6_F11_silent_none_detected": {"detected_direct": out["by_family"]["F11"]["detected_direct"], "detected_deployed": out["by_family"]["F11"]["detected_deployed"],
                                        "n": 10, "verdict": "MET" if none_any([r for r in rows if r["family"] == "F11"]) else "MISSED"},
        "P7_F11L_all_detected_by_label_isolation": {"detected_direct": out["by_family"]["F11L"]["detected_direct"], "detected_deployed": out["by_family"]["F11L"]["detected_deployed"],
                                                    "n": 10, "via_label_isolation": bool(f11l_via_iso),
                                                    "verdict": "MET" if both_all(f11l) and f11l_via_iso else "MISSED"},
        "P8_F12_all_detected_by_auc_roc": {"detected_direct": out["by_family"]["F12"]["detected_direct"], "detected_deployed": out["by_family"]["F12"]["detected_deployed"],
                                           "n": 3, "via_auc_roc_with_G5a_passing": bool(f12_via_auc),
                                           "verdict": "MET" if both_all(f12) and f12_via_auc else "MISSED"},
        "P9_G5_deployed_equals_direct": {"verdict": "MET" if g5_eq else "MISSED"},
    }
    out["grading"] = G
    out["F07_by_offset"], out["F08_by_angle"] = f07_by_offset, f08_by_angle
    out["F07_rank1_rmsd"] = {f"{r['params']['axis']}{r['params']['offset_angstrom']}": r["direct"].get("rank1_rmsd") for r in rows if r["family"] == "F07"}
    out["F08_rank1_rmsd"] = {f"{r['params']['rotation_axis']}{r['params']['rotation_degrees']}": r["direct"].get("rank1_rmsd") for r in rows if r["family"] == "F08"}
    out["n_met"] = sum(g["verdict"] == "MET" for g in G.values())
    out["seconds"] = round(time.time() - t0, 1)
    (ext2.OUT / "results.json").write_text(json.dumps(out, indent=1, default=str), encoding="utf-8")
    print(json.dumps(out["by_family"], indent=1))
    print("F07 by offset:", f07_by_offset)
    print("F08 by angle:", f08_by_angle)
    for k, g in G.items():
        print(f"{k:42s} {g['verdict']}")
    print("met", out["n_met"], "of 9;", out["seconds"], "s")


if __name__ == "__main__":
    main()
