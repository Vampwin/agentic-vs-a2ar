"""ADR-AIIM-427 (pre-registration) / ADR-AIIM-428 (results): the held-out and the seen faulted workspaces under the frozen evaluators (F), the
frozen evaluators with the four added checks called directly (P) and the deployed path on the copied contracts (D). Grades P1 to P12."""
from __future__ import annotations

import json
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import harness as H

SEED_HELD = 11          # F06 random choice of histidines
SEED_LEAK = 13          # held-out leaks
SEED_SEEN = 7           # the seeds of ADR-AIIM-416 and 420
BASES_G7 = [("SNP00163", ["r1", "r2", "r3"]), ("NPLIB01185", ["r1", "r2", "r3"]), ("ZM241385", ["r1", "r2", "r3"])]
HELD_FRACTIONS = (0.02, 0.30, 0.55, 0.75, 0.90, 0.94, 0.96, 0.99)
HELD_OFFSETS = (-100, -3, -2, -1, 1, 2, 3, 100)
HELD_VARIANTS = ("max_ifp_tanimoto", "asn253_filtered", "combined_zscore", "ifp_zm_only", "ifp_max_of_two")
PREREG_SHA256 = "46ca6a7cf054d8101771c424a9a66f7e8395957e6dcb437bc334d05207f773f2"       # replaced by the SHA-256 of PREREG_ADR427.md (CRLF form as stored) before the run is read


def leak_pool(kind: str, n: int, seed: int) -> list[str]:
    b = H.ext2.bench_frame()
    pool = b.loc[b["label"] == (1 if kind == "active" else 0), "compound_id"].tolist()
    return list(H.np.random.default_rng(seed).permutation(pool)[:n])


def g2g7_case(fam, gate, inp, ws, params, rec):
    return {"family": fam, "gate": gate, "params": params, "injection": rec, "frozen": H.run_frozen(gate, inp), "plus": H.run_plus(gate, inp),
            "deployed_plus": H.run_deployed_plus(gate, inp, ws)}


def g5_case(fam, td, tag, variant, params, rec, **kw):
    g5, reg = H.build_g5_plus(td / f"{tag}_f", variant, **kw)
    fr = H.run_frozen_g5(g5, reg)
    g5, reg = H.build_g5_plus(td / f"{tag}_p", variant, **kw)
    pl = H.run_plus_g5(g5, reg)
    g5d, _ = H.build_g5_plus(td / f"{tag}_d", variant, **kw)
    return {"family": fam, "gate": "G5", "variant": variant, "params": params, "injection": rec, "frozen": fr, "plus": pl,
            "deployed_plus": H.run_deployed_plus_g5(td / f"{tag}_d", g5d)}


def main():
    t0 = time.time()
    H.write_plus_contracts()
    snap = H.frozen_tree_hashes()
    rows = []
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        # ---- A: held-out F06, nine receptors other than the one the check was read from, one and three histidines at random
        for name, rec, cut in H.ext.G2_RECEPTORS[1:]:
            for n in (1, 3):
                ws = td / f"A_{name}_{n}"
                inp = H.build_g2_plus(ws, name, rec, cut)
                r = H.ext.inject_f06(ws, {"n_his_flipped": n, "residue_selection": "random"}, SEED_HELD)
                rows.append({"set": "A_heldout_F06", **g2g7_case("F06", "G2", inp, ws, {"receptor": name, "n_his_flipped": n}, r)})
        # ---- B and C: held-out F15 (silent truncation) and F16 (residue offsets)
        for sid, reps in BASES_G7:
            for f in HELD_FRACTIONS:
                ws = td / f"B_{sid}_{f}"
                inp = H.build_g7_plus(ws, sid, reps)
                r = H.ext.inject_f15(ws, {"fraction_kept": f}, log_updated=False)
                rows.append({"set": "B_heldout_F15", **g2g7_case("F15", "G7", inp, ws, {"base": sid, "fraction_kept": f}, r)})
            for o in HELD_OFFSETS:
                ws = td / f"C_{sid}_{o}"
                inp = H.build_g7_plus(ws, sid, reps)
                r = H.ext.inject_f16(ws, {"residue_offset": o})
                rows.append({"set": "C_heldout_F16", **g2g7_case("F16", "G7", inp, ws, {"base": sid, "residue_offset": o}, r)})
        # ---- D, E, F: held-out leaks of the benchmark's compound IDs to the ranking step
        for v in HELD_VARIANTS:
            for n in (1, 3, 10, 30):
                ids = leak_pool("active", n, SEED_LEAK)
                rows.append({"set": "D_heldout_F11_logged", **g5_case("F11", td, f"D_{v}_{n}", v, {"n_leaked": n}, {"n_leaked": len(ids)}, leak_ids=ids)})
            for n in (3, 30):
                ids = leak_pool("decoy", n, SEED_LEAK)
                rows.append({"set": "E_heldout_F11_decoy_ids", **g5_case("F11D", td, f"E_{v}_{n}", v, {"n_leaked": n}, {"n_leaked": len(ids)}, leak_ids=ids)})
            ids = leak_pool("active", 10, SEED_LEAK)
            rows.append({"set": "F_heldout_F11_unlogged", **g5_case("F11U", td, f"F_{v}", v, {"n_leaked": 10}, {"n_leaked": len(ids)}, leak_ids=ids, leak_mode="unlogged")})
        # ---- G: the faults of ADR-AIIM-416 and 420 (seen, which is why the checks exist), same bases and seeds
        base = H.ext.G2_RECEPTORS[0]
        for i, p in enumerate(H.ext.F06_GRID):
            ws = td / f"G_F06_{i}"
            inp = H.build_g2_plus(ws, base[0], base[1], base[2])
            r = H.ext.inject_f06(ws, p, SEED_SEEN)
            rows.append({"set": "G_seen", **g2g7_case("F06", "G2", inp, ws, p, r)})
        for fam, grid, upd in (("F15", H.ext.F15_GRID, False), ("F15L", H.ext.F15_GRID, True), ("F16", H.ext.F16_GRID, None)):
            for i, p in enumerate(grid):
                ws = td / f"G_{fam}_{i}"
                inp = H.build_g7_plus(ws, "SNP00157", ["r1", "r2", "r3"])
                r = H.ext.inject_f15(ws, p, log_updated=upd) if fam.startswith("F15") else H.ext.inject_f16(ws, p)
                rows.append({"set": "G_seen", **g2g7_case(fam, "G7", inp, ws, p, r)})
        for fam, from_bench in (("F11", False), ("F11L", True)):
            for i, p in enumerate(H.ext2.F11_GRID):
                ids = H.ext2.leaked_ids(p["n_leaked_ids"])
                rows.append({"set": "G_seen", **g5_case(fam, td, f"G_{fam}_{i}", "vina_only", p, {"n_leaked": len(ids)}, leak_ids=ids,
                                                       leak_mode="benchmark" if from_bench else "logged")})
    after = H.frozen_tree_hashes()
    anchor = json.loads((H.OUT / "anchor_clean.json").read_text("utf-8"))

    def st(r, k):
        v = r[k]
        return v["G5a"]["state"] if r["gate"] == "G5" else v["state"]

    def failed(r, k):
        v = r[k]
        return list(v["G5a"]["checks_failed"]) if r["gate"] == "G5" else list(v["checks_failed"])

    def added_of(r):
        return H.ADDED["G5a" if r["gate"] == "G5" else r["gate"]]

    for r in rows:
        r["vacuous"] = bool(r["injection"].get("vacuous"))
        r["detected_frozen"] = st(r, "frozen") != "PASS"
        r["detected_plus"] = st(r, "plus") != "PASS"
        r["detected_deployed"] = st(r, "deployed_plus") != "PASS"
        r["added_fired_plus"] = sorted(set(failed(r, "plus")) & set(added_of(r)))
        r["added_fired_deployed"] = sorted(set(failed(r, "deployed_plus")) & set(added_of(r)))
        dep = r["deployed_plus"]
        r["halted_deployed"] = bool((dep["G5a"] if r["gate"] == "G5" else dep).get("halted"))

    def summ(sel):
        rs = [r for r in rows if sel(r) and not r["vacuous"]]
        return {"n": len(rs), "n_vacuous_excluded": sum(1 for r in rows if sel(r) and r["vacuous"]),
                "detected_frozen": sum(r["detected_frozen"] for r in rs), "detected_plus": sum(r["detected_plus"] for r in rs),
                "detected_deployed": sum(r["detected_deployed"] for r in rs),
                "by_added_check_plus": sum(bool(r["added_fired_plus"]) for r in rs), "by_added_check_deployed": sum(bool(r["added_fired_deployed"]) for r in rs),
                "added_checks_that_fired": sorted({c for r in rs for c in r["added_fired_plus"]})}

    s = {
        "A": summ(lambda r: r["set"] == "A_heldout_F06"),
        "B_le_0.94": summ(lambda r: r["set"] == "B_heldout_F15" and r["params"]["fraction_kept"] <= 0.94),
        "B_ge_0.96": summ(lambda r: r["set"] == "B_heldout_F15" and r["params"]["fraction_kept"] >= 0.96),
        "C": summ(lambda r: r["set"] == "C_heldout_F16"),
        "D": summ(lambda r: r["set"] == "D_heldout_F11_logged"),
        "E": summ(lambda r: r["set"] == "E_heldout_F11_decoy_ids"),
        "F": summ(lambda r: r["set"] == "F_heldout_F11_unlogged"),
        "G_F06": summ(lambda r: r["set"] == "G_seen" and r["family"] == "F06"),
        "G_F15": summ(lambda r: r["set"] == "G_seen" and r["family"] == "F15"),
        "G_F15L": summ(lambda r: r["set"] == "G_seen" and r["family"] == "F15L"),
        "G_F16": summ(lambda r: r["set"] == "G_seen" and r["family"] == "F16"),
        "G_F11": summ(lambda r: r["set"] == "G_seen" and r["family"] == "F11"),
        "G_F11L": summ(lambda r: r["set"] == "G_seen" and r["family"] == "F11L"),
    }
    nv = [r for r in rows if not r["vacuous"]]
    det = [r for r in nv if r["detected_plus"]]

    def only(rs, cid):
        return all(r["added_fired_plus"] == [cid] for r in rs)

    def verdict(c):
        return "MET" if c else "MISSED"

    grading = {
        "P1_anchor": {"n_clean": anchor["n_clean"], "verdict": verdict(all(
            st(r, "frozen") == "PASS" and st(r, "plus") == "PASS" and st(r, "deployed_plus") == "PASS" for r in anchor["rows"]))},
        "P2_A": {**s["A"], "verdict": verdict(s["A"]["n"] == 16 and s["A"]["n_vacuous_excluded"] == 2 and s["A"]["detected_frozen"] == 0
                                              and s["A"]["detected_plus"] == 16 == s["A"]["detected_deployed"]
                                              and s["A"]["added_checks_that_fired"] == ["protonation_matches_record"])},
        "P3_B": {"le_0.94": s["B_le_0.94"], "ge_0.96": s["B_ge_0.96"],
                 "verdict": verdict(s["B_le_0.94"]["n"] == 18 == s["B_le_0.94"]["detected_plus"] == s["B_le_0.94"]["detected_deployed"]
                                    and s["B_le_0.94"]["detected_frozen"] == 0 and s["B_ge_0.96"]["n"] == 6 and s["B_ge_0.96"]["detected_plus"] == 0
                                    == s["B_ge_0.96"]["detected_deployed"] == s["B_ge_0.96"]["detected_frozen"]
                                    and s["B_le_0.94"]["added_checks_that_fired"] == ["trajectory_covers_plan"])},
        "P4_C": {**s["C"], "verdict": verdict(s["C"]["n"] == 24 == s["C"]["detected_plus"] == s["C"]["detected_deployed"] and s["C"]["detected_frozen"] == 0
                                              and s["C"]["added_checks_that_fired"] == ["binding_site_matches_receptor"])},
        "P5_D": {**s["D"], "verdict": verdict(s["D"]["n"] == 20 == s["D"]["detected_plus"] == s["D"]["detected_deployed"] and s["D"]["detected_frozen"] == 0
                                              and s["D"]["added_checks_that_fired"] == ["benchmark_ids_not_exposed"])},
        "P6_E": {**s["E"], "verdict": verdict(s["E"]["n"] == 10 == s["E"]["detected_plus"] == s["E"]["detected_deployed"] and s["E"]["detected_frozen"] == 0)},
        "P7_F": {**s["F"], "verdict": verdict(s["F"]["n"] == 5 and s["F"]["detected_plus"] == 0 == s["F"]["detected_deployed"] == s["F"]["detected_frozen"])},
        "P8_G": {k: v for k, v in s.items() if k.startswith("G_")},
        "P9_state_equal_P_and_D": {"verdict": verdict(all(st(r, "plus") == st(r, "deployed_plus") for r in rows))},
        "P10_detected_are_halted_under_D": {"verdict": verdict(all(r["halted_deployed"] for r in nv if r["detected_deployed"] and r["added_fired_deployed"]))},
        "P11_frozen_tree_unchanged": {"n_files": len(snap), "verdict": verdict(snap == after)},
    }
    g = s
    grading["P8_G"]["verdict"] = verdict(
        g["G_F06"]["n"] == 9 == g["G_F06"]["by_added_check_plus"] == g["G_F06"]["by_added_check_deployed"]
        and g["G_F15"]["n"] == 10 == g["G_F15"]["by_added_check_plus"] == g["G_F15"]["by_added_check_deployed"]
        and g["G_F15L"]["n"] == 10 == g["G_F15L"]["by_added_check_plus"] == g["G_F15L"]["by_added_check_deployed"]
        and g["G_F16"]["n"] == 10 == g["G_F16"]["by_added_check_plus"] == g["G_F16"]["by_added_check_deployed"]
        and g["G_F11"]["n"] == 10 == g["G_F11"]["by_added_check_plus"] == g["G_F11"]["by_added_check_deployed"]
        and g["G_F11L"]["n"] == 10 == g["G_F11L"]["detected_plus"] == g["G_F11L"]["detected_deployed"] == g["G_F11L"]["by_added_check_plus"])
    out = {"adr": "ADR-AIIM-427 (pre-registration) / ADR-AIIM-428 (results)", "prereg_sha256": PREREG_SHA256, "n_cases": len(rows),
           "n_clean": anchor["n_clean"], "summary": s, "grading": grading, "cases": rows}
    out["n_met"] = sum(v["verdict"] == "MET" for v in grading.values())
    out["seconds"] = round(time.time() - t0, 1)
    (H.OUT / "results.json").write_text(json.dumps(out, indent=1, default=str), encoding="utf-8")
    print(json.dumps(s, indent=1))
    for k, v in grading.items():
        print(f"{k:34s} {v['verdict']}")
    print("met", out["n_met"], "of", len(grading), ";", out["seconds"], "s")


if __name__ == "__main__":
    main()
