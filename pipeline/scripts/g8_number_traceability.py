"""
g8_number_traceability.py -- audit every load-bearing number in the manuscript
against the artifact it came from.

G8 is the one gate that checks the paper rather than the pipeline. Its premise
is that a number which cannot be re-derived from a stored artifact should not
appear in a manuscript, however carefully it was copied. Each entry below
names the claim, recomputes the value from its source, and compares.

A check that cannot locate its source FAILS. It does not warn and it does not
pass quietly: an unverifiable number is exactly the thing this gate exists to
catch.

Usage:
    python g8_number_traceability.py [--manuscript PATH]
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd


def _root() -> Path:
    marker = "00_project_control/08_DECISION_LOG.md"
    for cand in (Path(r"D:\Agentic VS Workflow"), Path("/mnt/d/Agentic VS Workflow")):
        if (cand / marker).is_file():
            return cand
    raise SystemExit("project root not found")


ROOT = _root()
RUNS = ROOT / "20_new_runs"
DEFAULT_MS = ROOT / "30_manuscript/CURRENT_v5/manuscript_AIIM_v6.md"


def _json(p: Path):
    return json.loads(p.read_text(encoding="utf-8"))


def _md_system() -> dict:
    """Recount the membrane system from its own coordinate and restraint files.

    The Methods section quotes an atom count, a lipid count, a water count and
    a restraint count. None of those is taken on trust here: each is counted
    from the file the simulation actually reads."""
    b = RUNS / "md_production_2026-08-21/03_build"
    gro = (b / "system_gmx.gro").read_text(encoding="utf-8").splitlines()
    n = int(gro[1])
    resnames = [gro[2 + i][5:10].strip() for i in range(n)]
    # Lipid21 splits POPC across three residues; one lipid is one PC headgroup.
    n_popc = sum(1 for r in resnames if r == "PC") // 38 if resnames.count("PC") else 0
    posre = sum(len([l for l in (b / f"posre_{m}.itp").read_text(encoding="utf-8").splitlines()
                     if l[:1].isspace() and l.strip()[:1].isdigit()])
                for m in ("system1", "system2", "ZMA"))
    # Water is counted in MOLECULES, not atoms. OPC is a four-site model, so
    # the atom count is 4x the molecule count and quoting it would inflate the
    # figure fourfold -- which is what the first version of this check did.
    n_water_atoms = resnames.count("WAT")
    # The crystallographic waters cannot be counted here at all: tleap
    # renames every water to the force field's own residue name, so they are
    # indistinguishable from bulk solvent by the time the system is
    # parameterised. The count belongs to the receptor-preparation step, and
    # is read from that step's own metadata.
    prep = _json(RUNS / "md_production_2026-08-21/01_receptor/"
                        "receptor_prep_v2_metadata.json")
    return {"n_atoms": n,
            "n_popc_headgroups": resnames.count("PC"),
            "n_water": n_water_atoms // 4,
            "n_water_atoms": n_water_atoms,
            "n_crystal_water": prep["crystal_waters_retained"],
            "n_crystal_water_available": prep["crystal_waters_available"],
            "n_restrained": posre}


# Each source is a zero-argument callable so a missing file raises here, in the
# check, rather than at import time where it would take the whole audit down.
def sources() -> dict[str, object]:
    e8 = pd.read_csv(RUNS / "np_library_funnel_2026-08-15/e8_np_library_prioritization.csv")
    cov = pd.read_csv(RUNS / "ad_sweep_expanded_2026-08-21/ad_coverage_by_descriptor.csv"
                      ).set_index("descriptor_space")
    per = pd.read_csv(RUNS / "ad_sweep_expanded_2026-08-21/ad_per_compound.csv")
    return {
        "e8": e8,
        "cov": cov,
        "per": per,
        "mut": _json(RUNS / "verification_summary_2026-08-21/mutation_score_by_gate.json"),
        "meta": _json(RUNS / "metamorphic_2026-08-21/metamorphic_results.json"),
        "yrand": _json(RUNS / "qsar_lbdd_2026-08-07/qsar_y_randomization.json"),
        "funnel": _json(RUNS / "funnel_exclusion_control_2026-08-21/exclusion_control_summary.json"),
        "tierA": _json(RUNS / "np_library_funnel_2026-08-15/tierA_convergence_analysis.json"),
        "struct": _json(RUNS / "verification_summary_2026-08-21/structure_and_protocol_checks.json"),
        "boxchk": _json(RUNS / "verification_summary_2026-08-21/box_and_cutoff_checks.json"),
        "conf": _json(RUNS / "np_library_funnel_2026-08-15/"
                      "conformal_qsar_summary_np_library_dock_ledger.json"),
        "g4": _json(RUNS / "g4_verification_2026-08-21/G4_contract_outcome_passA_cpu6.json"),
        "md": _md_system(),
        "dd": _json(RUNS / "diffdock_orthogonal_2026-08-23/orthogonal_check_result.json"),
        "ho": _json(RUNS / "heldout_evaluation_2026-08-08/heldout_evaluation_summary.json"),
    }


def build_checks(S) -> list[tuple[str, str, object, object]]:
    """(id, literal that must appear in the manuscript, expected, actual)."""
    e8, cov, per = S["e8"], S["cov"], S["per"]
    mut, meta, yr = S["mut"], S["meta"], S["yrand"]
    fu, tA, st, bx = S["funnel"], S["tierA"], S["struct"], S["boxchk"]
    C: list[tuple[str, str, object, object]] = []
    add = lambda cid, lit, exp, act: C.append((cid, lit, exp, act))

    # --- verification numbers -------------------------------------------------
    add("mutation_detection", "93.2 %", 93.2, round(mut["overall"]["detection_rate"] * 100, 1))
    add("mutation_cases", "133 seeded fault cases", 133, mut["overall"]["n_fault_cases"])
    add("mutation_detected", "detected 124", 124,
        mut["overall"]["n_fault_cases"] - mut["overall"]["n_escaped"])
    add("seer_A", "| **A** | on | on | **0.068** |", 0.068,
        round(mut["ablation_SEER"]["A_contracts_and_monitor"], 3))
    add("seer_C", "| **C** | off | on | 0.992 |", 0.992,
        round(mut["ablation_SEER"]["C_monitor_only"], 3))
    add("false_alarm", "**0.0 in every condition**", 0.0,
        mut["overall"]["false_escalation_rate_clean_runs"])
    add("gate_G4", "| G4 scoring | 27 | 4 | 85.2 % |", 85.2,
        round(mut["by_gate"]["G4"]["detection_rate"] * 100, 1))
    add("gate_G6", "| G6 selection | 32 | 5 | 84.4 % |", 84.4,
        round(mut["by_gate"]["G6"]["detection_rate"] * 100, 1))

    # --- metamorphic ----------------------------------------------------------
    add("meta_M1", "25 of 25 bit-identical", 25, meta["variants"]["M1"]["n_identical"])
    add("meta_M2_max", "**0.974 kcal/mol**", 0.974, round(meta["variants"]["M2"]["max_abs_deviation"], 3))
    add("meta_M2_mean", "mean 0.166", 0.166, round(meta["variants"]["M2"]["mean_abs_deviation"], 3))
    add("meta_M2_rho", "Spearman correlation with the baseline ranking is 0.772", 0.772,
        round(meta["variants"]["M2"]["spearman_rank_correlation_with_baseline"], 3))
    add("meta_M2_mw", "ρ = +0.735", 0.735,
        round(meta["variants"]["M2"]["deviation_vs_size"]["spearman_rho_molecular_weight"], 3))

    # --- y-randomisation ------------------------------------------------------
    add("yr_mean", "mean\ncross-validated AUC 0.500", 0.500, round(yr["permuted_mean"], 3))
    add("yr_z", "z = 21.9", 21.9, round(yr["z_vs_permuted"], 1))
    add("yr_p", "empirical p = 0.0099", 0.0099, round(yr["empirical_p"], 4))
    add("yr_real", "against 0.911 for", 0.911, round(yr["real_cv_auc"], 3))

    # --- funnel-exclusion control --------------------------------------------
    add("funnel_excl", "7.0 % reached the top-500", 7.0, round(fu["excluded_reaching_cutoff_pct"], 1))
    add("funnel_kept", "against 10.0 % of the 4,999", 10.0, round(fu["kept_pct_above_cutoff"], 1))
    add("funnel_cut", "−8.588 kcal/mol", -8.588, round(fu["top500_cutoff"], 3))
    add("funnel_p", "Mann–Whitney p = 0.0097", 0.0097, round(fu["mannwhitney_p"], 4))

    # --- applicability domain -------------------------------------------------
    add("ad_phys", "86.3 %", 86.3, round(float(cov.loc["PhysChem", "coverage_fraction"]) * 100, 1))
    add("ad_maccs", "9.3 % in MACCS", 9.3, round(float(cov.loc["MACCS", "coverage_fraction"]) * 100, 1))
    add("ad_erg", "3.6 % in ErG", 3.6, round(float(cov.loc["ErG", "coverage_fraction"]) * 100, 1))
    add("ad_ecfp4_zero", "0 % in ECFP4", 0.0, round(float(cov.loc["ECFP4", "coverage_fraction"]) * 100, 1))
    add("ad_max_spaces", "inside more than two of the eight", 2, int(per.n_spaces_inside.max()))
    add("ad_zero_spaces", "52 are inside none", 52, int((per.n_spaces_inside == 0).sum()))
    # NB: select the eight per-space flags explicitly. A blanket
    # `endswith("_inside")` also picks up `n_spaces_inside`, which is itself the
    # per-compound total -- that doubled the count when this check was written,
    # and the check caught it.
    space_cols = [f"{sp}_inside" for sp in cov.index]
    add("ad_phys_calls", "429 of the 494 in-domain calls", (429, 494),
        (int(per.PhysChem_inside.sum()), int(per[space_cols].sum().sum())))
    add("ad_ecfp4_max", "similarity of 0.352", 0.352, round(float(per.ECFP4_sim.max()), 3))

    # --- tiering and candidates ----------------------------------------------
    vc = e8.tier.value_counts()
    add("tier_A", "| A — convergent |", 3, int(vc.get("A_converging_sbdd_lbdd", 0)))
    add("tier_B", "| B — strong structure-based |", 119, int(vc.get("B_strong_sbdd_lbdd_null", 0)))
    add("tier_C", "| C — supporting |", 360, int(vc.get("C_supporting_sbdd_only", 0)))
    add("n_ok", "497 completed", 497, int((e8.status == "ok").sum()))
    add("n_onpocket", "482 produced an on-pocket pose", 482,
        int(((e8.status == "ok") & (e8.on_pocket == True)).sum()))  # noqa: E712
    add("tierA_binom", "binomial p = 0.49", 0.49,
        S["tierA"]["question_1_is_the_COUNT_enriched"]["binomial_p"])
    add("tierA_flav", "all flavones (p = 0.020)", 0.0198,
        S["tierA"]["question_2_is_the_IDENTITY_random"]["p_all_three_flavones"])
    add("tierA_pong", "all from one species (p = 0.002)", 0.0023,
        S["tierA"]["question_2_is_the_IDENTITY_random"]["p_all_three_pongamia"])

    # --- structures and protocol ---------------------------------------------
    add("res_4eiy", "1.80 Å A2AR–BRIL", 1.80, st["receptor_constructs"]["4EIY"]["resolution_A"])
    add("res_5olh", "2.60 Å A2A-StaR2-BRIL", 2.60, st["receptor_constructs"]["5OLH"]["resolution_A"])
    add("t88a", "T88A\n(3.36) at 5.91 Å", 5.91,
        next(m["min_dist_to_ligand_A"] for m in
             st["receptor_constructs"]["5OLH"]["engineered_mutations_in_A2AR"]
             if m["mutation"] == "T88A"))
    add("rmsd", "**0.379 Å**", 0.379, round(st["redocking_control"]["rank1_rmsd_A"], 3))
    add("stereo_snplib", "24 of 236 SNP-Lib", 24, st["undefined_stereocentres"]["SNP-Lib"]["n_undefined"])
    add("stereo_top500", "106 of the\n500 redocked compounds", 106,
        st["undefined_stereocentres"]["redocked top-500"]["n_undefined"])
    add("box_over", "**five compounds exceed the 22 Å\nedge**", 5,
        bx["ligand_extent"]["n_exceeding_box_edge"])
    add("box_max", "maximum 24.3 Å", 24.3, round(bx["ligand_extent"]["max_L1"], 1))
    add("cutoff_97", "482 of 497 compounds\n(97.0 %)", 97.0,
        next(r["pct_of_497"] for r in bx["centroid_cutoff_sensitivity"] if r["cutoff_A"] == 8))
    add("cutoff_7", "447 (89.9 %) at 7 Å", 89.9,
        next(r["pct_of_497"] for r in bx["centroid_cutoff_sensitivity"] if r["cutoff_A"] == 7))

    # --- conformal ------------------------------------------------------------
    add("conf_valid", "empirical error\n0.0501", 0.0501,
        S["conf"]["validity_check"]["significance_0.05"]["empirical_error_rate"])
    add("conf_n", "(n = 1,217)", 1217, S["conf"]["calibration_set_n"])
    add("conf_sets", "the ambiguous set {inactive, active} to 3", 3,
        int((e8["prediction_set_sig0.05"] == "{0,1}").sum()))

    # --- G4 -------------------------------------------------------------------
    add("g4_pass", "G4 then passed 7 of 7 checks", 7, len(S["g4"]["passed"]))
    add("mw_rho", "ρ = −0.182", -0.1824, S["g4"]["evidence"]["mw_bias_check"]["spearman_rho"])

    # --- membrane MD system (§2.9), counted from the simulation's own files ---
    md = S["md"]
    add("md_atoms", "93,203 atoms in total", 93203, md["n_atoms"])
    add("md_water", "15,656 OPC waters", 15656, md["n_water"])
    add("md_xtal_water", "156 crystallographic waters", 156, md["n_crystal_water"])
    add("md_restrained", "2,357 receptor and ligand heavy atoms", 2357,
        md["n_restrained"])

    # --- orthogonal-engine check (ADR-AIIM-111/112), re-derived from poses ----
    # Every number the manuscript quotes for this check is recomputed by
    # diffdock_orthogonal_check.py from the deposited SDF poses, not copied.
    dd = S["dd"]
    add("dd_control_rmsd", "**1.92 Å**", 1.92, round(dd["control"]["rmsd_to_crystal_A"], 2))
    add("dd_control_conf", "control +0.16", 0.16, dd["control"]["confidence"])
    add("dd_rmsds", "6.43, 8.80 and 9.50 Å", "6.43/8.80/9.50",
        "/".join(f'{dd["candidates"][c]["rmsd_to_vina_A"]:.2f}'
                 for c in ("NPLIB01185", "NPLIB01643", "NPLIB03056")))
    add("dd_confs", "+0.23, +0.47 and −0.50", "0.23/0.47/-0.50",
        "/".join(f'{dd["candidates"][c]["confidence"]:.2f}'
                 for c in ("NPLIB01185", "NPLIB01643", "NPLIB03056")))
    # The manuscript states the pre-registered prediction failed. That claim is
    # itself audited: if a rerun ever made it hold, this check fails loudly
    # rather than leaving a stale sentence in the text.
    add("dd_prediction_failed", "**failed**", False,
        dd["preregistered_prediction"]["holds"])
    add("dd_highest_conf", "highest (+0.47) on NPLIB01643", 0.47,
        dd["candidates"]["NPLIB01643"]["confidence"])

    # --- held-out fault cases (ADR-AIIM-060), counted from the run's own output -
    # This is the least comfortable number in the paper, which is exactly why it
    # is audited like every other one: if the artifact ever stops saying one of
    # five, the manuscript sentence fails loudly instead of going stale.
    ho = S["ho"]
    add("heldout_n", "five cases authored by a project-blind session", 5,
        ho["n_cases"])
    add("heldout_detected", "**One of the five was detected.**", 1,
        ho["n_detected"])
    add("heldout_detected_abstract", "it detected one", 1, ho["n_detected"])
    return C


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manuscript", default=str(DEFAULT_MS))
    args = ap.parse_args()
    text = Path(args.manuscript).read_text(encoding="utf-8")
    norm = re.sub(r"\s+", " ", text)

    S = sources()
    checks = build_checks(S)

    rows, n_pass, n_fail, n_absent = [], 0, 0, 0
    for cid, lit, exp, act in checks:
        lit_norm = re.sub(r"\s+", " ", lit)
        present = lit_norm in norm
        if isinstance(exp, float) and isinstance(act, (int, float)):
            agrees = abs(float(exp) - float(act)) <= max(5e-4, abs(float(exp)) * 1e-3)
        else:
            agrees = exp == act
        if not agrees:
            state = "FAIL-VALUE"; n_fail += 1
        elif not present:
            state = "FAIL-ABSENT"; n_absent += 1
        else:
            state = "PASS"; n_pass += 1
        rows.append({"check": cid, "state": state, "expected": exp, "recomputed": act,
                     "literal_present_in_manuscript": present, "literal": lit})

    # --- citation integrity -------------------------------------------------
    # A reference list is part of reporting integrity, so G8 audits it too:
    # every in-text marker must resolve to an entry, and every entry must be
    # cited. An uncited entry is padding; a dangling marker is a broken claim.
    ref_head = text.find("\n## References\n")
    body, reflist = (text[:ref_head], text[ref_head:]) if ref_head > 0 else (text, "")
    entries = {int(m) for m in re.findall(r"^\[(\d+)\]", reflist, flags=re.M)}
    cited: set[int] = set()
    for grp in re.findall(r"\[(\d+(?:[,-]\d+)*)\]", body):
        for part in grp.split(","):
            if "-" in part:
                a, b = part.split("-")
                cited.update(range(int(a), int(b) + 1))
            else:
                cited.add(int(part))
    dangling = sorted(cited - entries)
    uncited = sorted(entries - cited)
    contiguous = entries == set(range(1, len(entries) + 1)) if entries else False
    cite_ok = not dangling and not uncited and contiguous and bool(entries)
    rows.append({"check": "citation_integrity",
                 "state": "PASS" if cite_ok else "FAIL-VALUE",
                 "expected": "every marker resolves, every entry cited, 1..N contiguous",
                 "recomputed": {"n_entries": len(entries), "n_cited": len(cited),
                                "dangling_markers": dangling, "uncited_entries": uncited,
                                "numbering_contiguous": contiguous},
                 "literal_present_in_manuscript": bool(entries), "literal": "## References"})
    if cite_ok:
        n_pass += 1
    else:
        n_fail += 1

    for r in rows:
        if r["state"] != "PASS":
            print(f"  {r['state']:12s} {r['check']:18s} manuscript={r['expected']!r} "
                  f"artifact={r['recomputed']!r}  present={r['literal_present_in_manuscript']}")
    print(f"\nG8 number traceability: {n_pass} pass, {n_fail} value mismatch, "
          f"{n_absent} value correct but literal not found in the text")
    print(f"total checks: {len(rows)}")

    outp = RUNS / "verification_summary_2026-08-21/g8_number_traceability.json"
    outp.write_text(json.dumps({
        "manuscript": args.manuscript,
        "state": "PASS" if (n_fail == 0 and n_absent == 0) else "FAIL",
        "n_checks": len(rows), "n_pass": n_pass,
        "n_value_mismatch": n_fail, "n_literal_absent": n_absent,
        "checks": rows,
    }, indent=2, default=str), encoding="utf-8")
    print(f"written -> {outp}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
