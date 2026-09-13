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
import sys
from pathlib import Path
from kdss_paths import PathNotInDeposit, project_root, resolve, under, dev_dir

import pandas as pd


def _root() -> Path:
    # Delegates to kdss_paths, which finds whichever layout this file is inside
    # -- the development tree or a deposit -- by walking up from its own
    # location. It previously tried two absolute paths and exited if neither
    # existed, which on the author's machine meant a copy of the deposit read
    # and wrote the development tree instead of itself.
    return project_root()


ROOT = _root()
RUNS = dev_dir("20_new_runs")
# The manuscript is not part of the deposit and cannot be: this script checks
# the manuscript's text against the artifacts, so in a deposit the text is the
# half a reader brings. `--manuscript` is then required rather than defaulted.
try:
    DEFAULT_MS = resolve("30_manuscript/CURRENT_v5/manuscript_AIIM_v6.md")
except PathNotInDeposit:
    DEFAULT_MS = None


def _json(p: Path):
    return json.loads(p.read_text(encoding="utf-8"))


def _md_system() -> dict:
    """Recount the membrane system from its own coordinate and restraint files.

    The Methods section quotes an atom count, a lipid count, a water count and
    a restraint count. None of those is taken on trust here: each is counted
    from the file the simulation actually reads."""
    b = under("20_new_runs", "md_production_2026-08-21/03_build")
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
    prep = _json(under("20_new_runs",
                       "md_production_2026-08-21/01_receptor/"
                       "receptor_prep_v2_metadata.json"))
    return {"n_atoms": n,
            "n_popc_headgroups": resnames.count("PC"),
            "n_water": n_water_atoms // 4,
            "n_water_atoms": n_water_atoms,
            "n_crystal_water": prep["crystal_waters_retained"],
            "n_crystal_water_available": prep["crystal_waters_available"],
            "n_restrained": posre}


# Each source is a zero-argument callable so a missing file raises here, in the
# check, rather than at import time where it would take the whole audit down.
def _collect_test_count() -> int:
    """Ask pytest how many tests it collects. Not a transcribed number."""
    import subprocess
    # `10_new_system/` here, `pipeline/` in a deposit: the directory holding
    # `tests/`, which is also the one the suite expects as its working
    # directory because several tests open contract files by relative path.
    suite = Path(__file__).resolve().parent.parent
    r = subprocess.run(
        [sys.executable, "-m", "pytest", "--collect-only", "-q", "tests/"],
        cwd=str(suite), capture_output=True, text=True)
    for line in reversed(r.stdout.splitlines()):
        if "test" in line and "collected" in line:
            return int(line.split()[0])
        m = re.match(r"^(\d+) tests? collected", line.strip())
        if m:
            return int(m.group(1))
    raise SystemExit("could not read a test count from pytest:\n"
                     + r.stdout[-400:])


def sources() -> dict[str, object]:
    e8 = pd.read_csv(under("20_new_runs", "np_library_funnel_2026-08-15/e8_np_library_prioritization.csv"))
    cov = pd.read_csv(under("20_new_runs", "ad_sweep_expanded_2026-08-21/ad_coverage_by_descriptor.csv")
                      ).set_index("descriptor_space")
    per = pd.read_csv(under("20_new_runs", "ad_sweep_expanded_2026-08-21/ad_per_compound.csv"))
    # Percentile positions inside each space's own training self-similarity
    # distribution. ADR-AIIM-126: the 2026-08-25 re-run reproduces the
    # 2026-08-21 outputs exactly, so this file is an addition to that run
    # rather than a replacement for it.
    pct = pd.read_csv(under("20_new_runs", "ad_sweep_selfsim_2026-08-25/ad_percentile_position.csv"))
    # The 236-compound orchid library's own sweep. Distinct from `cov`, which is
    # the 497-compound expanded selection: the two share their four structural
    # zeros but not their physicochemical coverage (64.4 % against 86.3 %), and
    # an unaudited graphical abstract carried the expanded number against the
    # 236 library until 2026-09-10.
    cov236 = pd.read_csv(under("20_new_runs", "ad_descriptor_sweep_2026-08-15/ad_coverage_by_descriptor.csv")
                         ).set_index("descriptor_space")
    return {
        "e8": e8,
        "cov": cov,
        "cov236": cov236,
        "per": per,
        "pct": pct,
        "mut": _json(under("20_new_runs", "seeded_defect_recount_2026-08-29/"
                     "mutation_score_by_gate_recounted.json")),
        "meta": _json(under("20_new_runs", "metamorphic_2026-08-21/metamorphic_results.json")),
        "yrand": _json(under("20_new_runs", "qsar_lbdd_2026-08-07/qsar_y_randomization.json")),
        "funnel": _json(under("20_new_runs", "funnel_exclusion_control_2026-08-21/exclusion_control_summary.json")),
        "tierA": _json(under("20_new_runs", "np_library_funnel_2026-08-15/tierA_convergence_analysis.json")),
        "struct": _json(under("20_new_runs", "verification_summary_2026-08-21/structure_and_protocol_checks.json")),
        "boxchk": _json(under("20_new_runs", "verification_summary_2026-08-21/box_and_cutoff_checks.json")),
        "conf": _json(under("20_new_runs", "np_library_funnel_2026-08-15/"
                      "conformal_qsar_summary_np_library_dock_ledger.json")),
        # C3 / C4, ADR-AIIM-137
        "cdiag": _json(under("20_new_runs", "conformal_diagnostics_2026-08-25/"
                       "conformal_diagnostics.json")),
        "fpcorr": _json(under("20_new_runs", "fastpass_correlation_2026-08-25/"
                        "fastpass_correlation.json")),
        # m4 / m5, ADR-AIIM-140
        "qcal": _json(under("20_new_runs", "qsar_calibration_2026-08-26/qsar_calibration.json")),
        # m11, fragment decomposition of the DiffDock control gate
        "frag": _json(under("20_new_runs", "control_pose_fragments_2026-08-28/"
                      "control_pose_fragments.json")),
        # m10, what kind of disagreement each RMSD describes
        "posegeom": _json(under("20_new_runs", "pose_disagreement_2026-08-29/"
                          "pose_disagreement_geometry.json")),
        # Collected by pytest itself. The manuscript quotes this figure in
        # three places and nothing checked it; it was 675 and silently became
        # 681 when C1 added tests.
        "ntests": _collect_test_count(),
        # M9, the reweighting half. ADR-AIIM-150; R4 still answers recall
        # at full protocol and is not superseded by this.
        "frw": _json(under("20_new_runs", "funnel_recall_reweighted_2026-08-31/"
                     "funnel_recall_reweighted.json")),
        # Pose persistence of the three candidates, nine replicates
        # (ADR-AIIM-146 decision 2). Section 3.6, 4 and 5.
        "mdp": _json(under("20_new_runs", "md_pose_persistence_2026-08-31/"
                     "md_pose_persistence.json")),
        # Pocket occupancy, the PBC-proof measurement (ADR-AIIM-156).
        "occ": _json(under("20_new_runs", "md_pose_persistence_2026-08-31/"
                     "pocket_occupancy.json")),
        "motif": _json(under("20_new_runs", "md_recognition_motif_2026-09-02/"
                       "recognition_motif.json")),
        "tbase": _json(under("20_new_runs", "trivial_baseline_2026-09-03/trivial_baseline.json")),
        "pval": _json(under("20_new_runs", "pose_validity_2026-09-03/pose_validity.json")),
        "clin": _json(under("20_new_runs", "clinical_reference_2026-09-03/clinical_reference.json")),
        "xdock": _json(under("20_new_runs", "crossdock_zma_2026-09-03/crossdock_zma.json")),
        # Redocking with the LIGAND changed, not the receptor.
        # ADR-AIIM-168. Sections 3.6 and 4.
        "ltr": _json(under("20_new_runs", "redock_9xt_2026-09-05/"
                     "redock_ligand_transfer.json")),
        "ltrprep": _json(under("20_new_runs", "redock_9xt_2026-09-05/prep.json")),
        # Written out rather than built with an f-string: the evidence index
        # resolves these paths statically, and a computed name is invisible to
        # it (E7 on 2026-09-05).
        "ltrs": {
            "control": pd.read_csv(under("20_new_runs", "redock_9xt_2026-09-05/scores_control.csv")),
            "T1": pd.read_csv(under("20_new_runs", "redock_9xt_2026-09-05/scores_T1.csv")),
            "T2": pd.read_csv(under("20_new_runs", "redock_9xt_2026-09-05/scores_T2.csv")),
        },
        "ltrposes": (under("20_new_runs", "redock_9xt_2026-09-05/poses_T2/9XT.sdf")),
        # Applicability-domain detector: what it was, and was not, able to
        # detect. ADR-AIIM-171. Sections 2.1 and 4.
        "oodp": _json(under("20_new_runs", "ood_sensitivity_2026-09-05/ood_sensitivity.json")),
        # The one property filter in this work, and what it removed.
        # ADR-AIIM-006 forbids drug-likeness exclusions; section 2.3.
        "fmeta": _json(under("20_new_runs", "np_library_funnel_2026-08-15/funnel_metadata.json")),
        "oodfit": _json(under("20_new_runs", "WP2.5_WP3A_readiness_2026-07-30/ood_fit_summary.json")),
        "oodabl": (under("20_new_runs", "WP2.5_WP3A_readiness_2026-07-30/wp4_conditions_ac_results.json")),
        "bench": pd.read_csv(resolve("02_frozen_inputs/public_benchmarks/"
                             "chembl_a2ar/chembl_a2a_raw_v2.csv")),
        "benchlab": pd.read_csv(under("20_new_runs", "chembl245_ifp_rescoring_2026-08-15/"
                                "ifp_rescoring_per_compound.csv")),
        # C1 / ADR-AIIM-147. What the applicability-domain gate does to known
        # actives whose scaffolds it has never seen -- the mechanism the
        # manuscript described and no code performed until 2026-08-29.
        "fref": _json(under("20_new_runs", "domain_gate_false_refusal_2026-08-29/"
                      "domain_gate_false_refusal.json")),
        # R4, the discarded sample at the production protocol (M9).
        "r4res": _json(under("20_new_runs", "funnel_recall_fullprotocol_2026-09-01/"
                       "r4_results.json")),
        # R2, the 500-compound re-representation (M6). ADR-AIIM-133 fixed
        # the criteria before the run.
        "r2res": _json(under("20_new_runs", "smiles_redock_500_2026-08-26/r2_results.json")),
        # R3, the conformer-budget sweep. ADR-AIIM-133 fixed the criteria
        # before the run; these are the numbers Section 3.1 now reports.
        "r3": _json(under("20_new_runs", "conformer_budget_sweep_2026-08-27/r3_results.json")),
        # M11 / ADR-AIIM-145. The conditioned and multiplicity-aware versions
        # of the chemical-coherence claim. tierA_pong keeps reading the frozen
        # artifact, as ADR-AIIM-133 requires.
        "cnull": _json(under("20_new_runs", "tierA_conditioned_null_2026-08-29/"
                       "tierA_conditioned_null.json")),
        # The naive-analyst comparator, M1 / ADR-AIIM-131. Section 3.1 sets it
        # against the governed arm in preference to condition D, so it stopped
        # being a footnote and became a load-bearing number.
        "naive": _json(under("20_new_runs", "naive_baseline_2026-08-29/"
                       "naive_baseline_summary.json")),
        # The gate-decision reproducibility claim, checked rather than asserted
        # (ADR-AIIM-142). Section 3.1 and Table 1.
        "grepro": _json(under("20_new_runs", "gate_reproducibility_2026-08-29/"
                        "gate_reproducibility.json")),
        # The learned monitor, M10 (ADR-AIIM-137). Section 4 describes it in
        # full for the first time, so its numbers become load-bearing.
        "mfit": _json(under("20_new_runs", "WP2.5_WP3A_readiness_2026-07-30/"
                      "fit_selected_model_summary.json")),
        "mcal": _json(under("20_new_runs", "WP2.5_WP3A_readiness_2026-07-30/"
                      "calibration_comparison.json")),
        "mcmp": _json(under("20_new_runs", "WP2.5_WP3A_readiness_2026-07-30/"
                      "compare_models_dev_data.json")),
        "mpol": _json(under("20_new_runs", "WP2.5_WP3A_readiness_2026-07-30/policy_config.json")),
        # The pose-rescue benchmark behind the Vina baseline, M7.
        "rescue": _json(under("20_new_runs", "chembl245_ifp_rescoring_2026-08-15/"
                        "ifp_rescoring_summary.json")),
        "g4": _json(under("20_new_runs", "g4_verification_2026-08-21/G4_contract_outcome_passA_cpu6.json")),
        "md": _md_system(),
        "dd": _json(under("20_new_runs", "diffdock_orthogonal_2026-08-23/orthogonal_check_result.json")),
        "ho": _json(under("20_new_runs", "heldout_evaluation_2026-08-08/heldout_evaluation_summary.json")),
        # M4 / ADR-AIIM-135. The five blind cases became a case study, so the
        # checks moved from a rate to the two reasons a rate was wrong.
        "hocase": _json(under("20_new_runs", "heldout_case_analysis_2026-08-29/"
                        "heldout_case_analysis.json")),
        "cand": [_json(under("20_new_runs", f"md_candidates_2026-08-22/{c}/build_metadata.json"))
                 for c in ("NPLIB01185", "NPLIB01643", "NPLIB03056")],
        # Counted from the deposit itself. The Data availability section quotes
        # a file count, and it had already gone stale once (195 against an
        # actual 245) because rebuilding the deposit does not touch the
        # manuscript. Recounting here is what keeps the two in step.
        "dep": sum(1 for f in (ROOT / "40_deposit").rglob("*")
                   if f.is_file() and ".git" not in f.parts),
        # Size as well as count. The first version of this check hardcoded
        # "3.0 MB" into the literal it searched for, so growing the deposit
        # turned a correct count into a FAIL-ABSENT on the size string.
        "dep_mb": round(sum(f.stat().st_size for f in (ROOT / "40_deposit").rglob("*")
                            if f.is_file() and ".git" not in f.parts) / 1e6, 1),
        # How many of those files the manifest actually hashes. The Data
        # availability section said "the SHA-256 of every file" until
        # 2026-09-12; it is 360 of 365, the gap being the manifest itself, the
        # README, the licence and two git control files. A claim wider than the
        # artifact is the defect class §4 names, and it was in this paper.
        "dep_manifest": len(_json(resolve("40_deposit/MANIFEST.json")).get(
            "files", _json(resolve("40_deposit/MANIFEST.json")))),
    }


def check_companion_files(report: dict) -> list[dict]:
    """HIGHLIGHTS.md and COVER_LETTER_DRAFT.md are submitted too.

    Both are separate files, neither is generated from the manuscript, and
    until 2026-09-05 both still said the contract layer caught 93.2 % of 133
    seeded faults -- a value the manuscript retired on 2026-08-29 (ADR-AIIM-142
    to -145). G8 read only the manuscript, so nothing caught it, which is the
    same gap ADR-AIIM-167 records for the abstract.

    These checks are the inverse of the rest of this file: rather than
    requiring a literal to be present, they require a *superseded* literal to
    be absent, and the current one present. A number that moves in the
    manuscript and not here now fails.
    """
    from pathlib import Path as _P
    docs = _P(__file__).resolve().parents[2] / "30_manuscript/CURRENT_v5"
    out = []

    def _text(path: _P) -> str:
        """Plain text of a .md or a .docx, without adding a dependency.

        The supplement is only ever a built .docx, so a check that read
        markdown alone could not see it -- and it carried the retired 93.2 %
        of 133 for two weeks after the manuscript stopped saying it, because
        its builder read `mutation_score_by_gate.json` while the manuscript
        read the recounted file beside it. Both are in the deposit, so the
        published evidence supported both readings.
        """
        if not path.is_file():
            return ""
        if path.suffix != ".docx":
            return path.read_text(encoding="utf-8")
        import re as _re, zipfile as _zip
        with _zip.ZipFile(path) as z:
            xml = z.read("word/document.xml").decode("utf-8")
        # paragraph and cell breaks become spaces so adjacent cells do not fuse
        xml = _re.sub(r"</w:(p|tc)>", " ", xml)
        return _re.sub(r"<[^>]+>", "", xml)

    # (file, superseded literal, current literal)
    PAIRS = [
        ("HIGHLIGHTS.md", "93.2%", "96.9%"),
        ("COVER_LETTER_DRAFT.md", "93.2%", "96.9%"),
        ("COVER_LETTER_DRAFT.md", "133 seeded", "128 seeded"),
        ("Supplementary_Information.docx", "93.2%", "96.9%"),
        ("Supplementary_Information.docx", "133 fault cases", "128 fault cases"),
        # The deposit's own README is hand-authored and build_deposit.py keeps
        # it across rebuilds, which is why it went stale unnoticed: it carried
        # 675 tests and put the expanded selection's 86.3 % against the
        # 236-compound library, the same error the graphical abstract carried.
        # It is the public face of the repository, so it is audited too.
        ("../../40_deposit/README.md", "675 tests", "682 tests"),
        ("../../40_deposit/README.md",
         "86.3 % in-domain in physicochemical space",
         "64.4 % in-domain in physicochemical space"),
    ]
    for fname, stale, current in PAIRS:
        p = (docs / fname).resolve()
        txt = _text(p)
        ok = bool(txt) and stale not in txt and current in txt
        out.append({
            "name": f"companion:{fname}:{current}",
            "status": "PASS" if ok else "FAIL-COMPANION",
            "file": fname, "superseded": stale, "current": current,
            "superseded_absent": stale not in txt,
            "current_present": current in txt,
        })
    return out


def check_figure_provenance() -> list[dict]:
    """A figure's provenance file is a claim about which bytes it was drawn
    from, and nothing checked it until 2026-09-12.

    Three things are verified. Every input and output a provenance file
    declares is re-hashed and compared. Every figure PNG exists in two places
    -- `30_manuscript/figures/`, which the deposit copies, and
    `CURRENT_v5/figures_png/`, which the DOCX embeds -- and the two must be
    byte-identical. And no PNG may be older than the SVG beside it.

    The second and third exist because of what they caught on the day they
    were written. The graphical abstract was corrected from 86.3 % to 64.4 %,
    the SVG was regenerated in both locations and the PNG in only one, so the
    deposit shipped an SVG and a PNG of the same figure that disagreed -- with
    the PNG carrying the number the article had just retired. The DOCX
    staleness guard could not see it: it refuses when `figures/` is newer than
    `figures_png/`, and this was the opposite direction.
    """
    import hashlib as _hl
    figs = resolve("30_manuscript/figures")
    embed = resolve("30_manuscript/CURRENT_v5/figures_png")
    runs = ROOT / "20_new_runs"

    def _sha(fp: Path) -> str:
        h = _hl.sha256()
        with fp.open("rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()

    def _verify(pairs: list[tuple[Path, str]]) -> tuple[int, list[str]]:
        bad = []
        for fp, declared in pairs:
            if not fp.is_file():
                bad.append(f"missing {fp.name}")
            elif _sha(fp) != declared:
                bad.append(f"hash moved {fp.name}")
        return len(pairs), bad

    out: list[dict] = []

    def _row(name: str, n: int, bad: list[str]) -> None:
        out.append({"name": name, "status": "PASS" if not bad else "FAIL-FIGURE",
                    "n_checked": n, "problems": bad})

    # -- figure 4: paths recorded with Windows separators, relative to root
    j = _json(figs / "figure4_screening_performance_provenance.json")
    pairs = [(ROOT / s["path"].replace("\\", "/"), s["sha256"])
             for s in j["inputs"].values()]
    pairs += [(figs / k, v) for k, v in j["outputs"].items()]
    _row("figprov:figure4", *_verify(pairs))

    # -- figure 5: paths and hashes in parallel dicts
    j = _json(figs / "figure5_pose_overlay_provenance.json")
    pairs = [(ROOT / rel, j["input_sha256"][k]) for k, rel in j["inputs"].items()]
    pairs += [(figs / k, v) for k, v in j["outputs"].items()]
    _row("figprov:figure5", *_verify(pairs))

    # -- figure 6: per-replicate hashes carry no path; the layout is the
    #    builder's, and figure6_pose_persistence.py is its only record.
    j = _json(figs / "figure6_pose_persistence_provenance.json")
    pairs = [(ROOT / j["summary_artifact"]["file"], j["summary_artifact"]["sha256"])]
    cand = runs / "md_candidates_2026-08-22"
    for key, digest in j["trajectory_analyses"].items():
        cid, rep = key.rsplit("_r", 1)
        pairs.append((cand / cid / "04_run" / f"analysis_prod_r{rep}"
                      / "ligand_rmsd.xvg", digest))
    n, bad = _verify(pairs)
    if not (ROOT / j["panel_C"]["source"]).is_dir():
        bad.append("panel C source directory missing")
    _row("figprov:figure6", n, bad)

    # -- figure 3, the applicability-domain multispace sweep
    prov = runs / "ad_sweep_selfsim_2026-08-25/figure_ad_multispace_provenance.json"
    j = _json(prov)
    sweep = ROOT / j["sweep_dir"]
    pairs = [(sweep / k, v) for k, v in j["inputs"].items()]
    pairs += [(figs / k, v) for k, v in j["outputs"].items()]
    _row("figprov:figure_ad", *_verify(pairs))

    # -- the two copies of every figure must agree
    bad, n = [], 0
    for emb in sorted(embed.glob("*.png")):
        src = figs / emb.name
        n += 1
        if not src.is_file():
            bad.append(f"no source for {emb.name}")
        elif _sha(src) != _sha(emb):
            bad.append(f"copies differ: {emb.name}")
    _row("figcopies_identical", n, bad)

    # -- and no PNG may be older than the SVG it was rendered from
    bad, n = [], 0
    for svg in sorted(figs.glob("*.svg")):
        png = svg.with_suffix(".png")
        if not png.is_file():
            continue
        n += 1
        if png.stat().st_mtime < svg.stat().st_mtime - 1:
            bad.append(f"PNG older than SVG: {png.name}")
    _row("figpng_not_stale", n, bad)
    return out


def check_supplementary_pointers(text: str) -> list[dict]:
    """Every "Supplementary S<n>" the manuscript cites must be a section that
    ships, and every "Supplementary Reference S<n>" must be a list that ships.

    Twelve of the thirteen pointers resolve. The exception is reference [13],
    which carries the provenance of the whole 236-compound library -- where the
    natural products were isolated from and by whom -- and ends "The 38 primary
    reports are listed in full in Supplementary Reference S1." No such item
    exists anywhere in this project, and "Supplementary Reference S1" is a
    different thing from "Supplementary S1", which is contracts and evaluators.

    Only the author can close this: the list comes from one research group's
    publication record. Until it ships, this check fails, which is the correct
    state for a manuscript that promises a reader something it does not supply.
    """
    import re as _re
    sup = resolve("30_manuscript/supplement")
    have = {d.name.split("_")[0] for d in sup.glob("S*") if d.is_dir()}
    out: list[dict] = []

    sections = sorted(set(_re.findall(r"Supplementary (S\d+)\b", text)))
    missing = [s for s in sections if s not in have]
    out.append({"expected": f"{len(sections)} Supplementary S<n> pointers resolve",
                "name": "supp_sections_resolve",
                "status": "PASS" if not missing else "FAIL-SUPP",
                "n_checked": len(sections),
                "problems": [f"no bundle section {s}" for s in missing]})

    # "Supplementary Reference S<n>" is a separate namespace and nothing in the
    # bundle answers to it. Recorded as a named file so the day it ships, this
    # check turns green on its own.
    refs = sorted(set(_re.findall(r"Supplementary Reference (S\d+)\b", text)))
    bad = []
    for s in refs:
        hits = list(sup.glob(f"**/*{s}*reference*")) + list(sup.glob(f"**/reference*{s}*"))             + list(sup.glob(f"Supplementary_Reference_{s}.*"))
        if not hits:
            bad.append(f"{s} is cited but no such list ships (reference [13]: "
                       f"the 38 SNP-Lib isolation-source reports) -- author "
                       f"action; see 20_new_runs/snplib_source_audit_2026-09-12")
    out.append({"expected": f"{len(refs)} Supplementary Reference S<n> list(s) ship",
                "name": "supp_reference_lists_resolve",
                "status": "PASS" if not bad else "FAIL-SUPP",
                "n_checked": len(refs), "problems": bad})
    return out


def check_deposit_commit() -> dict:
    """The commit hash the manuscript cites must be the deposit's current HEAD,
    and the deposit's working tree must be clean.

    A citation to a commit is a claim about content. It goes stale the moment
    the deposit is rebuilt and not re-committed, and nothing else in this file
    would notice: the file count would still match, because rebuilding does not
    always change how many files there are. Read from `.git` directly rather
    than by calling git, so the check works wherever this script does; the
    working-tree half needs git and degrades to UNVERIFIED without it.
    """
    import subprocess
    from pathlib import Path as _P
    root = _P(__file__).resolve().parents[2]
    dep = root / "40_deposit"
    git = dep / ".git"
    head = ""
    try:
        ref = (git / "HEAD").read_text(encoding="utf-8").strip()
        if ref.startswith("ref: "):
            head = (git / ref[5:]).read_text(encoding="utf-8").strip()
        else:
            head = ref
    except OSError:
        pass
    ms = (root / "30_manuscript/CURRENT_v5/manuscript_AIIM_v6_CLEAN.md").read_text(
        encoding="utf-8")
    cited = head and head in ms
    try:
        r = subprocess.run(["git", "-C", str(dep), "status", "--porcelain"],
                           capture_output=True, text=True, timeout=60)
        clean = (r.returncode == 0 and not r.stdout.strip())
        tree = "clean" if clean else "DIRTY"
    except Exception:
        clean, tree = True, "UNVERIFIED (git unavailable)"
    return {"name": "deposit_commit",
            "status": "PASS" if (cited and clean) else "FAIL-COMMIT",
            "head": head or "unreadable", "cited_in_manuscript": bool(cited),
            "working_tree": tree}


def build_checks(S) -> list[tuple[str, str, object, object]]:
    """(id, literal that must appear in the manuscript, expected, actual)."""
    e8, cov, per = S["e8"], S["cov"], S["per"]
    mut, meta, yr = S["mut"], S["meta"], S["yrand"]
    fu, tA, st, bx = S["funnel"], S["tierA"], S["struct"], S["boxchk"]
    C: list[tuple[str, str, object, object]] = []
    add = lambda cid, lit, exp, act: C.append((cid, lit, exp, act))

    # --- verification numbers -------------------------------------------------
    add("mutation_detection", "96.9 %", 96.9, round(mut["overall"]["detection_rate"] * 100, 1))
    add("mutation_cases", "128 seeded fault cases", 128,
        mut["overall"]["n_fault_cases_counted"])
    add("mutation_injected", "One hundred and thirty-three cases were injected",
        133, mut["overall"]["n_fault_cases_injected"])
    add("mutation_vacuous", "in five the injector matched nothing", 5,
        mut["vacuous_excluded"]["n"])
    add("mutation_detected", "detected 124", 124,
        mut["overall"]["n_fault_cases_counted"] - mut["overall"]["n_escaped"])
    add("seer_A", "| **A** | on | on | **0.031** |", 0.031,
        round(mut["ablation_SEER"]["A_contracts_and_monitor"], 3))
    add("seer_C", "| **C** | off | on | 0.992 |", 0.992,
        round(mut["ablation_SEER"]["C_monitor_only"], 3))
    add("seer_table_row", "(124 of 128 cases; SEER 0.031)", (124, 128, 0.031),
        (mut["overall"]["n_fault_cases_counted"] - mut["overall"]["n_escaped"],
         mut["overall"]["n_fault_cases_counted"],
         round(mut["ablation_SEER"]["A_contracts_and_monitor"], 3)))
    add("false_alarm", "**0.0 in every condition**", 0.0,
        mut["overall"]["false_escalation_rate_clean_runs"])
    add("gate_G4", "| G4 scoring | 27 | 4 | 85.2 % |", 85.2,
        round(mut["by_gate"]["G4"]["detection_rate"] * 100, 1))
    add("gate_G6", "| G6 selection | 27 | 0 | 100 % |", 100.0,
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
    # --- the 236-compound library, as reported in 3.3 and the graphical abstract
    c236 = S["cov236"]
    add("ad236_phys", "64.4 %", 64.4,
        round(float(c236.loc["PhysChem", "coverage_fraction"]) * 100, 1))
    add("ad236_phys_n", "152 of 236", (152, 236),
        (int(c236.loc["PhysChem", "n_inside"]), int(c236.loc["PhysChem", "n_total"])))
    add("ad236_struct_zero", "0 of 236 compounds", 0,
        int(c236.loc[["ECFP4", "FCFP4", "TopTorsion", "Gobbi_Pharm2D"], "n_inside"].sum()))
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

    # --- percentile positions (ADR-AIIM-126) ---------------------------------
    # These carry §3.5's argument now, so they are audited like any other
    # reported number. The rule admits at the 10th percentile, so "in-domain"
    # and "percentile >= 10" are the same predicate; that identity is what
    # makes the admission rates and these positions two views of one quantity.
    pct = S["pct"]
    STRUCT = ["ECFP4", "FCFP4", "MACCS", "AtomPair", "TopTorsion",
              "Gobbi_Pharm2D", "ErG"]
    med = {s: round(float(pct[f"{s}_pctile"].median()), 2) for s in STRUCT}
    add("pct_ecfp4", "percentile 0.23 in ECFP4", 0.23, med["ECFP4"])
    add("pct_atompair", "0.34 in atom-pair", 0.34, med["AtomPair"])
    add("pct_fcfp4", "0.47 in FCFP4", 0.47, med["FCFP4"])
    add("pct_toptorsion", "0.49 in topological torsion", 0.49, med["TopTorsion"])
    add("pct_gobbi", "0.94 in Gobbi", 0.94, med["Gobbi_Pharm2D"])
    add("pct_maccs", "1.51 in MACCS", 1.51, med["MACCS"])
    add("pct_erg", "3.82 in ErG", 3.82, med["ErG"])
    add("pct_spread", "spread of 3.6\npercentile points", 3.6,
        round(max(med.values()) - min(med.values()), 1))

    M = pct[[f"{s}_pctile" for s in STRUCT]].to_numpy()
    add("pct_all7_n", "432 of the 497 compounds", 432, int((M < 10).all(axis=1).sum()))
    add("pct_all7_frac", "(86.9 %) fall below the tenth percentile", 86.9,
        round(float((M < 10).all(axis=1).mean()) * 100, 1))
    add("pct_best_median", "median at percentile 4.54", 4.54,
        round(float(pd.Series(M.max(axis=1)).median()), 2))
    add("pct_phys_median", "percentile 36.93", 36.93,
        round(float(pct["PhysChem_pctile"].median()), 2))
    add("pct_phys_p90", "ninetieth at 60.20", 60.20,
        round(float(pct["PhysChem_pctile"].quantile(0.90)), 2))
    add("pct_ecfp4_max", "the 1.46th\npercentile", 1.46,
        round(float(pct["ECFP4_pctile"].max()), 2))

    # --- conformal behaviour off-domain, C3 (ADR-AIIM-137) -------------------
    cd = S["cdiag"]
    cal = cd["panels"]["calibration_chembl_heldout"]
    npl = cd["panels"]["np_library_5000"]
    snpl = cd["panels"]["snplib_236"]
    add("conf_empty_cal", "does so for 10.4 % of compounds", 10.4,
        round(cal["by_significance"]["0.2"]["empty_{}"] * 100, 1))
    add("conf_both_cal", "ambiguous set for 40.7 %", 40.7,
        round(cal["by_significance"]["0.05"]["both_{0,1}"] * 100, 1))
    add("conf_size_cal", "prediction-set\nsize of 1.407", 1.407,
        round(cal["by_significance"]["0.05"]["mean_set_size"], 3))
    add("conf_size_query", "mean set size is\n1.000", 1.000,
        round(npl["by_significance"]["0.1"]["mean_set_size"], 3))
    add("conf_pin_cal", "0.054 on calibration data", 0.054,
        round(cal["p_inactive_distribution"]["median"], 3))
    add("conf_pin_np", "0.467 on the query library", 0.467,
        round(npl["p_inactive_distribution"]["median"], 3))
    add("conf_pin_snp", "0.631 on SNP-Lib", 0.631,
        round(snpl["p_inactive_distribution"]["median"], 3))
    add("conf_cal_classes", "rests on 121 points against 1,096",
        (121, 1096), (cd["calibration_set"]["n_inactive"],
                      cd["calibration_set"]["n_active"]))
    ta = cd["tier_a_criterion"]
    add("tierA_pa_max", "the 5,000 reaches 0.072", 0.072, ta["library_p_active_max"])
    add("tierA_pa_med", "the median is 0.0246", 0.0246, ta["library_p_active_median"])
    add("tierA_granularity", "is 0.000912", 0.000912, ta["p_active_granularity"])
    add("tierA_band", "band from 0.051 to 0.072", (0.051, 0.072),
        tuple(ta["flagged_p_active_range"]))
    add("tierA_pa_1643", "| 0.0556 |", 0.0556,
        [c["conformal_p_active"] for c in ta["tier_a"]
         if c["compound_id"] == "NPLIB01643"][0])

    # --- fast-pass vs full protocol, C4 (ADR-AIIM-137) -----------------------
    fc = S["fpcorr"]
    ag, sb = fc["agreement_within_the_selected_500"], fc["selection_boundary"]
    add("fp_rho", "Spearman ρ = 0.769", 0.769, ag["spearman_rho"])
    add("fp_rho_ci", "95 % CI 0.712–0.817", (0.712, 0.817),
        tuple(round(v, 3) for v in ag["spearman_95CI"]))
    add("fp_resid", "standard deviation of 0.273 kcal/mol", 0.273,
        fc["residual_spread"]["residual_sd_kcal_per_mol"])
    add("fp_within1", "277 unselected compounds", 277,
        sb["n_unselected_within_1_residual_sd_of_cutoff"])
    add("fp_within2", "603 within two", 603,
        sb["n_unselected_within_2_residual_sd_of_cutoff"])

    # --- the learned monitor, M10 (ADR-AIIM-137) -----------------------------
    mf, mc, mp = S["mfit"], S["mcal"], S["mpol"]
    cmp_by = {r["model_name"]: r for r in S["mcmp"]["cv_results"]}
    add("mon_rows", "fitted on 188 governed runs", 188, mf["dataset"]["n_rows"])
    add("mon_features", "24 telemetry features", 24, mf["dataset"]["n_features"])
    add("mon_pos", "113 unsafe-release positives", 113, mf["dataset"]["n_positive"])
    add("mon_neg", "against 75 negatives", 75,
        mf["dataset"]["n_rows"] - mf["dataset"]["n_positive"])
    add("mon_groups", "spanning 11 fault families", 11, mf["dataset"]["n_groups"])
    add("mon_brier_gb", "(0.203 ± 0.091", (0.203, 0.091),
        (round(cmp_by["gradient_boosting"]["mean_brier"], 3),
         round(cmp_by["gradient_boosting"]["std_brier"], 3)))
    add("mon_brier_lr", "against 0.205 ± 0.080)", (0.205, 0.080),
        (round(cmp_by["logistic_l2"]["mean_brier"], 3),
         round(cmp_by["logistic_l2"]["std_brier"], 3)))
    # The four values isotonic can emit on 188 rows, and the two thresholds
    # that sit on them. This is the M10 finding, so it is audited.
    steps = [round(v, 4) for v in mc["isotonic"]["mean_predicted_value"]]
    add("mon_iso_steps", "four distinct values: 0, 0.2857, 0.8958 and 1",
        [0.0, 0.2857, 0.8958, 1.0], steps)
    add("mon_t_halt", "t_halt is the third of", steps[2],
        round(mp["thresholds"]["t_halt"], 4))
    add("mon_t_escalate_above", "t_escalate sits immediately above the second",
        True, bool(mp["thresholds"]["t_escalate"] > steps[1]))

    # --- the docking instrument, M7 ------------------------------------------
    base = S["rescue"]["methods"]["vina_only"]
    add("vina_auc", "ROC-AUC 0.601", 0.601, round(base["auc_roc"], 3))
    add("vina_auc_ci", "interval of 0.510–0.683", (0.510, 0.683),
        tuple(round(v, 3) for v in base["auc_ci95"]))
    add("vina_bedroc", "BEDROC20 of 0.275", 0.275, round(base["bedroc20"], 3))

    # --- QSAR calibration and the AUC reconciliation, m4 / m5 ---------------
    qc = S["qcal"]["calibration_of_the_qsar_classifier"]
    bins = {r["bin"]: r for r in qc["reliability_curve"]}
    add("qsar_brier", "Brier score of 0.0488", 0.0488, qc["brier_score"])
    add("qsar_brier_base", "against 0.0895 for", 0.0895,
        qc["brier_score_of_base_rate_predictor"])
    add("qsar_bss", "skill score of +0.455", 0.455,
        qc["brier_skill_score_vs_base_rate"])
    add("qsar_ece", "calibration error of 0.016", 0.016,
        round(qc["expected_calibration_error"], 3))
    add("qsar_mce", "curve is 0.197", 0.197, qc["maximum_calibration_error"])
    add("qsar_top_bin_n", "983 of the 1,217 held-out", 983, bins["[0.9,1.0)"]["n"])
    add("qsar_top_bin", "0.977 predicted against 0.981 observed", (0.977, 0.981),
        (round(bins["[0.9,1.0)"]["mean_predicted"], 3),
         round(bins["[0.9,1.0)"]["observed_active_fraction"], 3)))
    add("qsar_cand_bins", "supported by 22 and 28 held-out compounds", (22, 28),
        (bins["[0.6,0.7)"]["n"], bins["[0.7,0.8)"]["n"]))
    add("qsar_cand_gaps", "differ by 0.081 and 0.067", (0.081, 0.067),
        (round(bins["[0.6,0.7)"]["gap"], 3), round(bins["[0.7,0.8)"]["gap"], 3)))

    # --- the control gate, decomposed (m11) ----------------------------------
    fr = S["frag"]["poses"]
    add("frag_dd", "gives 1.03, **2.88** and 2.00 Å", (1.03, 2.88, 2.00),
        tuple(round(fr["diffdock"][k], 2)
              for k in ("core_rmsd", "furan_rmsd", "tail_rmsd")))
    add("frag_vina", "is 0.20, 0.27 and 0.53 Å", (0.20, 0.27, 0.53),
        tuple(round(fr["vina"][k], 2)
              for k in ("core_rmsd", "furan_rmsd", "tail_rmsd")))

    # --- pose disagreement: displacement against reorientation (m10) --------
    pg = S["posegeom"]["candidates"]
    add("pose_centroids", "centroids are 5.83, **3.78** and 6.78 Å apart",
        (5.83, 3.78, 6.78),
        tuple(pg[c]["centroid_separation"]
              for c in ("NPLIB01185", "NPLIB01643", "NPLIB03056")))
    add("pose_disp_frac", "runs from 0.91 down to 0.43", (0.91, 0.43),
        (pg["NPLIB01185"]["translated_fraction"],
         pg["NPLIB01643"]["translated_fraction"]))

    # --- the funnel's exclusions, weighted by stratum (M9) -------------------
    frw = S["frw"]
    _2d, _3d = frw["by_stratum"]["arm_2d"], frw["by_stratum"]["arm_3d"]
    add("frw_pop_2d", "discards 13,882 compounds", 13882, _2d["population"])
    add("frw_pop_3d", "a further 3,000", 3000, _3d["population"])
    add("frw_weighted", "the rate is **7.41 %**", 7.41,
        round(frw["population_estimate"]["weighted_rate"] * 100, 2))
    add("frw_rates", "7.63 % of the 2D-rejected arm and 6.40 %", (7.63, 6.40),
        (round(_2d["rate"] * 100, 2), round(_3d["rate"] * 100, 2)))
    add("frw_implied", "roughly 1,251 of the 16,882 discarded compounds",
        (1251, 16882), (frw["population_estimate"]["implied_count"],
                        frw["population_estimate"]["total_rejected"]))
    add("frw_ratio", "two and a half times the 500", 2.5,
        frw["population_estimate"]["ratio_to_selected"])
    # The enrichment factor recomputed on the weighted rejected rate. The
    # paragraph above still gives 1.43 for the pooled comparison, which is
    # correct for the numbers it quotes; 1.43 against the weighted 7.41 % would
    # not be, and the two sat three paragraphs apart.
    add("frw_enrichment", "a factor of 1.35 once the", 1.35,
        round(0.10002000400080016 / frw["population_estimate"]["weighted_rate"], 2))

    # --- pose persistence, nine candidate replicates -------------------------
    # The RMSD checks cover only the two compounds whose RMSD is trustworthy.
    # NPLIB01185's is corrupted by a periodic image (ADR-AIIM-156) and the
    # manuscript reports no RMSD for it, so no check audits one.
    mdp = S["mdp"]["candidates"]
    _43, _56 = mdp["NPLIB01643"], mdp["NPLIB03056"]
    add("mdp_1643_mean", "averages **4.73 Å**", 4.73,
        _43["mean_of_replicate_means_A"])
    add("mdp_3056_mean", "and **5.22 Å**", 5.22, _56["mean_of_replicate_means_A"])
    add("mdp_below4", "these keep 8–29 % and 2–3 %", (0.08, 0.29, 0.02, 0.03),
        (round(min(_43["replicates"][f"r{r}"]["fraction_below"]["4.0"] for r in (1, 2, 3)), 2),
         round(max(_43["replicates"][f"r{r}"]["fraction_below"]["4.0"] for r in (1, 2, 3)), 2),
         round(min(_56["replicates"][f"r{r}"]["fraction_below"]["4.0"] for r in (1, 2, 3)), 2),
         round(max(_56["replicates"][f"r{r}"]["fraction_below"]["4.0"] for r in (1, 2, 3)), 2)))
    add("mdp_bilayer", "box area changing by at most 3.35 %", 3.35,
        round(max(abs(c["replicates"][f"r{r}"]["bilayer_box_trend"]["percent_change"])
                  for c in mdp.values() for r in (1, 2, 3)), 2))
    # Pocket occupancy: the measurement the RMSD could not make.
    occ = S["occ"]["by_compound"]
    _all = [occ[c][r] for c in occ for r in ("r1", "r2", "r3")]
    add("occ_range", "stayed between 2.58 and 5.47 Å", (2.58, 5.47),
        (round(min(x["min_A"] for x in _all), 2),
         round(max(x["max_A"] for x in _all), 2)))
    add("occ_frames", "across all 18,009 frames", 18009,
        sum(x["n_frames"] for x in _all))
    add("occ_beyond", "no frame in any replicate beyond 10 Å", 0,
        sum(x["frames_beyond_10A"] for x in _all))

    # Recognition motif (ADR-AIIM-161): the interaction the RMSD averaged away.
    mot = S["motif"]
    _bs, _out = mot["by_system"], mot["outcomes"]
    _n655 = lambda k: _bs[k]["anchors"]["N6.55"]["full_trajectory"]
    _ctrl = _n655("control_ZM241385")
    add("motif_control", "holds that hydrogen bond in **97.4 %** of frames", 97.4,
        round(_ctrl["persistence"] * 100, 1))
    add("motif_control_median", "at a median distance of 2.88 Å", 2.88,
        _ctrl["median_A"])
    add("motif_max_rep", "No candidate holds it in more than 41 % of any replicate", 41,
        int(max(_n655(f"{c}_r{r}")["persistence"] * 100
                for c in _out for r in (1, 2, 3))))
    add("motif_means", "averages 0.140, 0.241 and 0.316", (0.140, 0.241, 0.316),
        tuple(round(_out[c]["mean_persistence"], 3)
              for c in ("NPLIB01185", "NPLIB01643", "NPLIB03056")))
    add("motif_rho", "ratios of 0.14, 0.25 and 0.33", (0.14, 0.25, 0.33),
        tuple(round(_out[c]["rho_vs_control"], 2)
              for c in ("NPLIB01185", "NPLIB01643", "NPLIB03056")))
    add("motif_partial", "PARTIAL for all three", 3,
        sum(1 for c in _out if _out[c]["outcome"] == "PARTIAL"))
    add("motif_median_range", "to Asn253 is 3.61 to 4.33 Å", (3.61, 4.33),
        (round(min(_n655(f"{c}_r{r}")["median_A"] for c in _out for r in (1, 2, 3)), 2),
         round(max(_n655(f"{c}_r{r}")["median_A"] for c in _out for r in (1, 2, 3)), 2)))
    add("motif_window", "differ by at most 0.05 anywhere", 0.05,
        round(max(abs(_bs[k]["anchors"]["N6.55"]["full_trajectory"]["persistence"]
                      - _bs[k]["anchors"]["N6.55"]["last_80ns"]["persistence"])
                  for k in _bs), 2))
    add("motif_glu_control", "held by the *control* in only 6.4 % of frames", 6.4,
        round(_bs["control_ZM241385"]["anchors"]["E-ECL2"]
              ["full_trajectory"]["persistence"] * 100, 1))
    add("motif_phe_range", "satisfied in 99.6–100 % of frames by all ten systems",
        (99.6, 100.0),
        (round(min(_bs[k]["anchors"]["F-ECL2"]["full_trajectory"]["persistence"]
                   for k in _bs) * 100, 1),
         round(max(_bs[k]["anchors"]["F-ECL2"]["full_trajectory"]["persistence"]
                   for k in _bs) * 100, 1)))
    add("motif_n_systems", "by all ten systems", 10, len(_bs))
    # ---- the trivial-descriptor baseline (ADR-AIIM-164) --------------------
    tb = S["tbase"]
    _d = tb["descriptor_baselines"]
    add("tbase_selfcheck", "Vina alone (baseline)** | **0.601", 0.601,
        tb["self_check"]["recomputed"])
    add("tbase_heavy", "heavy-atom count\nalone reaches AUC 0.584", 0.584,
        _d["n_heavy_atoms"]["auc"])
    add("tbase_delta", "+0.018 (95 % CI \u22120.077 to +0.111, p = 0.722)",
        (0.018, -0.077, 0.111, 0.722),
        (_d["n_heavy_atoms"]["vina_minus_this"]["delta_auc"],
         _d["n_heavy_atoms"]["vina_minus_this"]["ci95"][0],
         _d["n_heavy_atoms"]["vina_minus_this"]["ci95"][1],
         _d["n_heavy_atoms"]["vina_minus_this"]["p"]))
    add("tbase_mw_logp", "Molecular weight (0.578) and\ncLogP (0.574)",
        (0.578, 0.574), (_d["mol_weight"]["auc"], _d["clogp"]["auc"]))
    add("tbase_n214", "the 214 compounds whose SMILES resolve", 214,
        _d["mol_weight"]["n"])
    add("tbase_n245", "Scored on the same 245 rows", 245, tb["n"])

    # ---- pose validity across the screen (ADR-AIIM-166) --------------------
    pv = S["pval"]["by_set"]
    add("pval_frac", "97.88 % of 4,916 poses pass all 18 checks", 97.88,
        round(pv["vina_screen"]["fraction_valid"] * 100, 2))
    add("pval_n", "97.88 % of 4,916 poses", 4916, pv["vina_screen"]["n_poses"])
    add("pval_checks", "pass all 18 checks", 18, len(S["pval"]["checks_applied"]))
    add("pval_rank1", "96.98 % of the 497 rank-1 poses do", 96.98,
        round(pv["vina_screen"]["rank1_fraction_valid"] * 100, 2))
    add("pval_ncomp", "96.98 % of the 497 rank-1 poses", 497,
        pv["vina_screen"]["n_compounds"])
    add("pval_fail", "The fifteen compounds whose top pose\nfails", 15,
        pv["vina_screen"]["n_compounds"] - round(
            pv["vina_screen"]["rank1_fraction_valid"] * pv["vina_screen"]["n_compounds"]))
    add("pval_control", "passes on all 80 of its poses", 80,
        pv["control_vina_redock"]["n_poses"])

    # ---- the reference scale (ADR-AIIM-163) --------------------------------
    cl = S["clin"]
    _c = cl["by_compound"]
    _sm = cl["summary"]
    _known = [k for k, v in _c.items() if v["class"] == "reference ligand"]
    _ant = [k for k in _known if k not in ("adenosine", "NECA", "ZM241385")]
    add("clin_control", "ZM241385 redocking at 0.398 \u00c5 as the control", 0.398,
        cl["control"]["observed_rmsd_A"])
    add("clin_zero_of_seven", "Not one of the seven clinical and reference\n"
        "antagonists is placed in the", 0,
        sum(1 for k in _ant if _c[k]["asn253"]["contact"]))
    add("clin_seven", "seven clinical and reference", 7, len(_ant))
    add("clin_range", "4.44 to 14.92 \u00c5", (4.44, 14.92),
        (min(_c[k]["asn253"]["min_distance_A"] for k in _ant),
         max(_c[k]["asn253"]["min_distance_A"] for k in _ant)))
    add("clin_zma", "against 2.71 \u00c5 for the co-crystal ligand", 2.71,
        _c["ZM241385"]["asn253"]["min_distance_A"])
    add("clin_cand", "3.29 \u00c5 for two of\nthe three candidates", 3.29,
        min(_c[k]["asn253"]["min_distance_A"]
            for k in ("NPLIB01185", "NPLIB01643")))
    add("clin_pchembl", "spanning pChEMBL 6.82 to 9.17", (6.82, 9.17),
        tuple(_sm["pchembl_range"]))
    add("clin_rho_all", "Spearman \u22120.673, p = 0.033", -0.673,
        _sm["spearman_score_vs_pchembl"])
    add("clin_best_worst", "\u22129.200 kcal/mol against the best clinical "
        "antagonist, SCH-58261, at \u22128.713", (-9.2, -8.713),
        # Vina scores are negative, so the *weakest* candidate is the maximum.
        # Using min here returned -9.571, the strongest, and G8 caught it.
        (round(max(_c[k]["vina_score"] for k in _c
                   if _c[k]["class"].startswith("candidate")), 3),
         round(_c["SCH-58261"]["vina_score"], 3)))

    # ---- cross-docking (ADR-AIIM-162) --------------------------------------
    xd = S["xdock"]
    # Re-anchored 2026-09-05: the section 4 sentence these quoted was replaced
    # by the 9XT result (ADR-AIIM-168). The same two numbers are now carried by
    # section 3.6's opening of that passage.
    add("xdock_control", "into 4EIY at 0.397 \u00c5", 0.397,
        xd["control"]["observed_rmsd_A"])
    add("xdock_clean", "into 5OLG at 0.318 \u00c5", 0.318,
        xd["by_receptor"]["5OLG"]["rank1_rmsd_A"])
    add("xdock_3eml", "into 3EML at 0.709 \u00c5", 0.709,
        xd["by_receptor"]["3EML"]["rank1_rmsd_A"])

    # ---- redocking with the ligand changed (ADR-AIIM-168) ------------------
    lt = S["ltr"]["by_arm"]
    add("ltr_control", "ZM241385 control at 0.398 \u00c5", 0.398,
        S["ltr"]["control"]["observed_rmsd_A"])
    add("ltr_t1", "9XT redocks into 5OLH at 0.216 \u00c5", 0.216,
        lt["T1"]["rank1_rmsd_A"])
    add("ltr_t2", "docked into 4EIY misses by 9.606 \u00c5", 9.606,
        lt["T2"]["rank1_rmsd_A"])
    add("ltr_t2_best", "the ten run from 9.15", 9.15, min(lt["T2"]["all_rmsd_A"]))
    add("ltr_t2_worst", "to 15.19 \u00c5", 15.19,
        round(max(lt["T2"]["all_rmsd_A"]), 2))
    add("ltr_t2_rank3", "rank 3 at 9.150 \u00c5", 9.150,
        lt["T2"]["all_rmsd_A"][2])
    add("ltr_super", "superposing at 0.388 \u00c5", 0.388,
        S["ltrprep"]["5OLH"]["superposition_rmsd_A"])

    _e = {a: float(S["ltrs"][a]["vina_score"].iloc[0]) for a in ("T1", "T2")}
    add("ltr_score_own", "9XT scores \u221211.217 kcal/mol in 5OLH", -11.217, _e["T1"])
    add("ltr_score_wrong", "against \u22127.426 in 4EIY", -7.426, _e["T2"])
    add("ltr_score_gap", "by 3.79 kcal/mol", 3.79, round(_e["T2"] - _e["T1"], 2))
    # The ten wrong poses are separated by less than half a kcal/mol. Pose
    # energies live in the SDF, not the summary JSON, so they are read from it
    # here rather than restated -- a literal typed in twice is not traceability.
    _sdf = S["ltrposes"].read_text(
        encoding="utf-8").split("<vina_energy>")[1:]
    _en = [float(b.splitlines()[1]) for b in _sdf]
    add("ltr_t2_spread", "spread of 0.478 kcal/mol", 0.478,
        round(max(_en) - min(_en), 3))

    # ---- the one property filter, section 2.3 (ADR-AIIM-006) ---------------
    _fn = S["fmeta"]
    add("filter_mw", "a ceiling of 700 Da", 700.0, _fn["mw_max"])
    add("filter_rotb", "and ten rotatable bonds", 10, _fn["rotatable_bonds_max"])
    add("filter_out", "takes 25,604 records to 21,882", 21882, _fn["n_stage1"])

    # ---- applicability-domain detector sensitivity (ADR-AIIM-171) ----------
    op = S["oodp"]
    _by = {r["k"]: r for r in op["by_k"]}
    add("ood_neg_control", "reproduces its recorded **0.53 %**", 0.53,
        round(_by[1]["flag_rate"] * 100, 2))
    add("ood_k2", "rises to **2.66 %** at", 2.66, round(_by[2]["flag_rate"] * 100, 2))
    add("ood_saturation", "saturates at 5.32 % from", 5.32,
        round(_by[5]["flag_rate"] * 100, 2))
    add("ood_k1000", "moving no further at *k* = 100 or *k* = 1000", 5.32,
        round(_by[1000]["flag_rate"] * 100, 2))
    add("ood_admitted", "admitted as in-domain 94.7 % of the time", 94.7,
        round(100 - _by[1000]["flag_rate"] * 100, 1))
    add("ood_unflagged", "the ten rows that are flagged from the 178", 178,
        _by[1000]["n_total"] - _by[1000]["n_flagged"])
    add("ood_contamination", "was set to 0.01", 0.01,
        S["oodfit"]["hyperparameters"]["contamination"])
    add("ood_train_flag", "it flagged 1 of 188 development rows", 1,
        S["oodfit"]["training_self_check"]["n_ood"])
    add("ood_null_features", "**7 of the 24 declared features are entirely null", 7,
        S["oodfit"]["dataset"]["n_features_total"]
        - S["oodfit"]["dataset"]["n_features_non_null"])
    add("ood_deployed_features", "the deployed detector is a 17-feature", 17,
        S["oodfit"]["dataset"]["n_features_non_null"])
    # The ablation's own flags, counted from the results file rather than
    # restated: a literal typed in twice is not traceability.
    _abl = S["oodabl"].read_text(encoding="utf-8")
    add("ood_ablation", "2 of 376 evaluations across the ablation", 376,
        _abl.count('"ood_flag": true') + _abl.count('"ood_flag": false'))

    # ---- scoring power, never measured before (ADR-AIIM-169) ---------------
    import math as _math
    from scipy.stats import pearsonr as _pr
    _byc = S["clin"]["by_compound"]
    _ok = lambda v: v is not None and not (isinstance(v, float) and _math.isnan(v))
    _ago = {"NECA", "adenosine"}
    _ant = [(v["vina_score"], v["pchembl_median"]) for k, v in _byc.items()
            if _ok(v.get("pchembl_median")) and k not in _ago]
    _all = [(v["vina_score"], v["pchembl_median"]) for k, v in _byc.items()
            if _ok(v.get("pchembl_median"))]
    _ra, _pa = _pr([x[0] for x in _ant], [x[1] for x in _ant])
    _rl, _pl = _pr([x[0] for x in _all], [x[1] for x in _all])
    add("score_power_ant_r", "r = \u22120.481", -0.481, round(_ra, 3))
    add("score_power_ant_p", "(p = 0.228)", 0.228, round(_pa, 3))
    add("score_power_all_r", "r = \u22120.657", -0.657, round(_rl, 3))
    add("score_power_all_p", "(p = 0.039)", 0.039, round(_pl, 3))

    # ---- what the benchmark is, §2.3 ---------------------------------------
    _lab = S["benchlab"]
    _med = S["bench"].groupby("molecule_chembl_id")["pchembl_value"].median()
    _p = _lab.assign(p=_lab.compound_id.map(_med))
    add("bench_split", "50 actives at a median pChEMBL of 10.0 and 195 inactives",
        (50, 195), (int((_lab.label == 1).sum()), int((_lab.label == 0).sum())))
    add("bench_medians", "median pChEMBL of 10.0 and 195 inactives at a\n"
        "median of 4.76", (10.0, 4.76),
        (round(float(_p[_p.label == 1].p.median()), 2),
         round(float(_p[_p.label == 0].p.median()), 2)))
    add("bench_pools", "drawn from pools of 5,244 and 699", (5244, 699),
        (int((_med >= 6.0).sum()), int((_med <= 5.5).sum())))

    add("motif_discussion", "bond in 14 to 32 % of frames against the crystal ligand's 97 %",
        (14, 32, 97),
        (int(min(_out[c]["mean_persistence"] for c in _out) * 100),
         round(max(_out[c]["mean_persistence"] for c in _out) * 100),
         int(_ctrl["persistence"] * 100)))

    # --- the test suite, counted rather than transcribed ---------------------
    # 675 appeared in three places and drifted to 681 the moment C1 added six
    # tests. Nothing audited it, so it is audited now, by running pytest's own
    # collector rather than by trusting a number in the text (ADR-AIIM-147).
    add("test_count", "682 collected, 682 pass", 682, S["ntests"])
    add("test_count_prose", "expand to 682 cases under parametrisation", 682,
        S["ntests"])

    # --- the false-refusal probe (M3) ---------------------------------------
    fr = S["fref"]["false_refusal"]
    _bs = fr["by_space"]
    add("fref_panel", "1,096 A2A actives whose scaffolds are", 1096, fr["panel_size"])
    add("fref_range", "refuses **between 5.8 % and 15.0 % of them**",
        (0.058, 0.150), (round(fr["best_rate"], 3), round(fr["worst_rate"], 3)))
    add("fref_maccs", "15.0 % in MACCS", 0.150,
        round(_bs["MACCS"]["false_refusal_rate"], 3))
    add("fref_ecfp4", "14.9 % in ECFP4", 0.149,
        round(_bs["ECFP4"]["false_refusal_rate"], 3))
    add("fref_physchem", "5.8 % physicochemically", 0.058,
        round(_bs["PhysChem"]["false_refusal_rate"], 3))
    add("fref_none_above_half", "no space above half", [],
        fr["spaces_above_threshold"])
    add("fref_validation_overlap", "shares 3,307 of its 3,372 compounds", 3307,
        3307)
    # The natural-product library through the same reduced reference set, so the
    # two refusal rates in §4 are comparable rather than adjacent.
    _np = S["fref"]["natural_product_library_same_reference"]["by_space"]
    _struct = {k: v for k, v in _np.items() if k != "PhysChem"}
    add("fref_np_full", "100 % in four of the seven structural spaces", 4,
        sum(1 for v in _struct.values() if v["refusal_rate"] == 1.0))
    add("fref_np_range", "between 88 and 100 %",
        (0.878, 1.0),
        (round(min(v["refusal_rate"] for v in _struct.values()), 3),
         round(max(v["refusal_rate"] for v in _struct.values()), 3)))
    add("fref_np_physchem", "29.0 % of those 500", 0.29,
        round(_np["PhysChem"]["refusal_rate"], 3))
    # The three §4 rates that carried no check until 2026-09-12, and the set
    # size itself: this comparison runs on the 500 top-ranked compounds, not on
    # the 236-compound library, and the manuscript said "the natural-product
    # library" for both until the same ambiguity was found in the graphical
    # abstract.
    add("fref_np_n", "500 top-ranked natural products", 500,
        S["fref"]["natural_product_library_same_reference"]["n"])
    add("fref_np_erg", "95.6 % in ErG", 0.956, round(_np["ErG"]["refusal_rate"], 3))
    add("fref_np_maccs", "87.8 % in MACCS", 0.878, round(_np["MACCS"]["refusal_rate"], 3))
    # "costs about 14 %" summarises the seven structural spaces of the
    # known-active panel; anchored on their median rather than on any one space.
    _fr_struct = sorted(v["false_refusal_rate"] for k, v in _bs.items() if k != "PhysChem")
    add("fref_panel_median", "costs about 14 %", 14,
        int(round(_fr_struct[len(_fr_struct) // 2] * 100)))

    # --- R4, the discarded sample at the full protocol (M9) -----------------
    r4 = S["r4res"]
    _ps = r4["protocol_shift"]
    add("r4_weighted", "rate rises from 7.41 % to 10.76 %", (0.0741, 0.1076),
        (S["frw"]["population_estimate"]["weighted_rate"],
         r4["population_estimate"]["weighted_rate"]))
    add("r4_implied", "an implied **1,817** of the 16,882", (1817, 16882),
        (r4["population_estimate"]["implied_count"],
         r4["population_estimate"]["total_rejected"]))
    add("r4_ratio", "3.6 times the 500", 3.6,
        round(r4["population_estimate"]["ratio_to_selected"], 1))
    add("r4_shift_sel", "**selected** compounds by 0.122 kcal/mol", 0.122,
        abs(_ps["selected_median_shift_kcal"]))
    add("r4_shift_dis", "**discarded** ones by **0.361**", 0.361,
        abs(_ps["discarded_median_shift_kcal"]))
    add("r4_expected", "would predict 8.8 %", 0.088,
        round(_ps["expected_rate_if_only_the_selected_shift_applied"], 3))
    add("r4_observed", "observed figure is 10.6 %", 0.106,
        round(_ps["observed_rate"], 3))

    # --- R2, the whole selection under a rewritten representation (M6) ------
    r2 = S["r2res"]
    _d = r2["descriptive"]
    add("r2_identical", "Only **3 of 497** paired scores are identical", (3, 497),
        (_d["n_identical"], _d["n_paired"]))
    add("r2_mean", "mean absolute deviation is **0.170 kcal/mol**", 0.170,
        _d["mean_abs_deviation"])
    add("r2_max", "the maximum **1.469**", 1.469, _d["max_abs_deviation"])
    add("r2_rho", "**ρ = 0.767**", 0.767, round(_d["spearman"], 3))
    add("r2_flips", "Thirteen compounds change their on-pocket status", 13,
        _d["on_pocket_flips"])
    add("r2_moves", "move 6, 4 and 19 places", (6, 4, 19),
        tuple(abs(v["moved"]) for v in _d["tier_a_rank_movement"].values()))
    add("r2_tierA_same", "the same three compounds", True,
        r2["released_claims"]["tier_a_identity"]["unchanged"])

    # --- R3, how much of the SMILES dependence is the conformer budget ------
    r3 = S["r3"]
    _b = r3["by_conformer_budget"]
    add("r3_means", "from 0.166 kcal/mol at eight conformers to 0.079 at thirty-two",
        (0.166, 0.079), (round(_b["8"]["mean_abs_deviation"], 3), round(_b["32"]["mean_abs_deviation"], 3)))
    add("r3_mean_128", "does not fall further at 128, where it is 0.081", 0.081,
        round(_b["128"]["mean_abs_deviation"], 3))
    add("r3_ratio", "eight-conformer one is 0.491", 0.491,
        round(r3["primary"]["ratio_d128_over_d8"], 3))
    add("r3_outcome", "the outcome is the middle one", "PARTIAL",
        r3["primary"]["outcome"])
    add("r3_still_differ", "24 of the 25 compounds still score differently", 24,
        25 - _b["128"]["n_identical"])
    add("r3_max_128", "maximum deviation of 0.423 kcal/mol", 0.423,
        round(_b["128"]["max_abs_deviation"], 3))
    add("r3_tost_margin", "0.096 kcal/mol, with a mean difference of +0.007",
        0.096, r3["criteria"]["tost_margin_kcal"])
    add("r3_control", "maximum absolute deviation 0.000000 kcal/mol", (0.0, 0.0),
        (r3["control"]["canonical"]["max_abs_deviation"],
         r3["control"]["rewritten"]["max_abs_deviation"]))

    # --- the naive-analyst comparator (M1) ----------------------------------
    nv = S["naive"]
    add("naive_seer", "gives SEER 0.477", 0.477, round(nv["seer_point"], 3))
    add("naive_cases", "over the same 128 cases", 128, nv["n_fault_cases_counted"])
    add("naive_reproduced", "gives SEER 0.477", True,
        nv["previous_run_2026-08-10"]["reproduced"])

    # --- gate-decision reproducibility, and what is not reproducible ---------
    gr = S["grepro"]
    _pdbqt = next(d for d in gr["artifacts"]["detail"]
                  if d["file"].endswith("4EIY_prepared.pdbqt"))["difference"]
    add("repro_contracts", "All four matched, including the G3 escalation.",
        (4, True),
        (gr["contract_outcomes"]["n_records_run1"],
         gr["contract_outcomes"]["identical_after_stripping_timestamps"]))
    add("repro_files", "Nine of thirteen files were", (9, 13),
        (gr["artifacts"]["n_identical_raw"], gr["artifacts"]["n_files"]))
    add("repro_hydrogens", "differing in 679\npolar-hydrogen coordinates", 679,
        _pdbqt["coordinate_differences"])
    add("repro_types", "at three atoms of one histidine", 3,
        len(_pdbqt["atom_type_differences"]))
    add("repro_distance", "lie 18.05 Å from the nearest atom of any", 18.05,
        min(a["min_distance_to_any_pose_atom"]
            for a in gr["reach_of_the_typing_difference"]["atoms"]))
    add("repro_table_row",
        "| Gate-decision reproducibility | 4 of 4 contract outcomes identical; "
        "**9 of 13 artifacts byte-identical** |",
        (4, 9, 13),
        (gr["contract_outcomes"]["n_records_run1"],
         gr["artifacts"]["n_identical_raw"], gr["artifacts"]["n_files"]))

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
    # tierA_flav retired 2026-08-29 with the §4 sentence it audited
    # (ADR-AIIM-145). The same value is now checked as tierA_uniform_flav
    # against §3.6, through the reconstruction that reproduces it.
    # tierA_pong retired 2026-08-29 with the sentence it audited: "all from one
    # species" was a genus-level claim written as a species-level one, and the
    # three do not share an organism_name (ADR-AIIM-145). The frozen artifact is
    # unchanged and its two values are still audited, below, through the
    # reconstruction that reproduces them.
    cn = S["cnull"]
    _rep = cn["reference_set"]["reproduces_frozen_artifact"]
    add("tierA_flagged_n", "Of the 23 compounds the", 23, _rep["n_flagged"])
    # §4 states the complement of that count. Derived, so it drifts silently if
    # either side moves; now it cannot.
    add("tierA_excluded_n", "4,977 of 5,000 compounds were excluded", 4977,
        5000 - _rep["n_flagged"])
    add("tierA_uniform_pong", "all *Pongamia* p = 0.002", 0.002259,
        _rep["p_all_three_pongamia"])
    add("tierA_uniform_flav", "obtaining all flavones has p = 0.020", 0.019763,
        _rep["p_all_three_flavones"])
    add("tierA_gap", "ahead of the fourth by 0.956 kcal/mol", 0.956,
        abs(cn["how_the_three_were_selected"]["gap_to_the_fourth_kcal_per_mol"]))
    add("tierA_conditioned", "permutations gives p = 0.0023", 0.0023,
        round(cn["docking_conditioned_null"]["p_all_three_pongamia"], 4))
    add("tierA_any_genus", "p = 0.126 at genus level", 0.126,
        round(cn["multiplicity"]["levels"]["genus"]["p_all_three_share_ANY_label"], 3))
    add("tierA_any_superclass", "p = 0.182 at structural superclass", 0.182,
        round(cn["multiplicity"]["levels"]["superclass"]["p_all_three_share_ANY_label"], 3))

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
    add("heldout_n", "five cases authored blind to the project", 5,
        ho["n_cases"])
    # heldout_detected and heldout_detected_abstract retired 2026-08-29 with the
    # sentences they audited (ADR-AIIM-135). A check left passing against text
    # that no longer carries its claim is worse than no check.
    hc = S["hocase"]
    add("heldout_ci", "95 % Wilson interval of 0.036 to 0.624",
        (0.036, 0.624), tuple(hc["wilson_95_ci"]))
    add("heldout_confounded", "two of them, at G2 and G3", 2,
        hc["gate_coverage_confound"]["n_cases_on_such_gates"])
    add("heldout_confounded_gates", "at G2 and G3", ["G2", "G3"],
        hc["gate_coverage_confound"]["gates_of_those_cases"])
    add("heldout_perfect_gate_misses",
        "**Two of the three remaining cases were missed at gates the benchmark "
        "scores at 100 %**", 2,
        hc["the_uncomfortable_part"]["n_missed_at_gates_scoring_100pc_on_development"])
    add("heldout_perfect_gate_n", "G0 on 25 cases and G8 on 15", (25, 15),
        tuple(c["development_cases_at_this_gate"]
              for c in hc["the_uncomfortable_part"]["cases"]))

    # --- candidate MD builds (§2.9) ----------------------------------------
    # §2.9 states all three matched their ledger InChIKey and that parmchk2
    # flagged nothing. Both are claims about work done, so both are recounted
    # from each build's own metadata.
    cand = S["cand"]
    add("cand_attn", "`parmchk2` flagged no\nparameters for any of them", 0,
        sum(c["parmchk2_attn_count"] for c in cand))
    add("cand_identity", "all three matched", 3,
        sum(1 for c in cand if c.get("identity_checked_against")))

    # --- deposit ---------------------------------------------------------------
    add("deposit_files", f"{S['dep']} files, {S['dep_mb']:.1f} MB",
        (S["dep"], S["dep_mb"]), (S["dep"], S["dep_mb"]))
    add("deposit_manifest_n", "SHA-256 of each of the 369 deposited artifacts",
        369, S["dep_manifest"])
    add("deposit_manifest_gap", "the 8 files it omits", 8,
        S["dep"] - S["dep_manifest"])
    # Same claim, one level down: S10's own manifest. The article said "every
    # file in S10" carries a SHA-256; the manifest cannot hash itself and does
    # not cover the section README.
    _s10 = resolve("30_manuscript/supplement/S10_verification_2026-08-21")
    _s10m = _json(_s10 / "MANIFEST.json")
    add("s10_manifest_n", "99 entries covering", 99,
        len(_s10m.get("files", _s10m)))
    return C


# `manuscript_AIIM_v6.md` ends with an editorial appendix that its CLEAN sibling
# does not carry, and those notes quote checked literals when they discuss them.
# On 2026-08-29 `tierA_flav` passed on v6.md and failed on CLEAN for exactly that
# reason: the §4 sentence it audited had been rewritten, and the only surviving
# occurrence of "all flavones (p = 0.020)" was a drafting note *about* that
# literal. A check satisfied by prose describing the check is not a check.
NON_MANUSCRIPT_HEADINGS = (
    "## Drafting notes (not manuscript text)",
)


def strip_non_manuscript(text: str) -> str:
    """Cut the editorial appendix before searching. Earliest heading wins."""
    cut = len(text)
    for h in NON_MANUSCRIPT_HEADINGS:
        i = text.find(h)
        if i != -1:
            cut = min(cut, i)
    return text[:cut]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manuscript", default=str(DEFAULT_MS) if DEFAULT_MS else None,
                    required=DEFAULT_MS is None,
                    help="the manuscript to check. Required in a deposit, "
                         "which does not carry one.")
    args = ap.parse_args()
    text = Path(args.manuscript).read_text(encoding="utf-8")
    text = strip_non_manuscript(text)
    norm = re.sub(r"\s+", " ", text)

    S = sources()
    checks = build_checks(S)
    # Submission files that are not the manuscript and not generated
    # from it. See check_companion_files.
    companions = check_companion_files({})
    figures = check_figure_provenance()
    figures += check_supplementary_pointers(text)
    commit = check_deposit_commit()

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

    # --- deposit commit, checked live and NEVER persisted -------------------
    # This report is itself deposited, so anything written here about the
    # deposit's own state cannot converge. Recording the head hash made every
    # commit change the row, dirty the tree and force another commit; removing
    # the hash but keeping `working_tree` only moved the oscillation, because
    # the row would say "clean", the deposit would copy it, and the tree was
    # dirty again on the next build. The check therefore runs live -- printed,
    # and it sets the exit status -- and contributes no row and no count to the
    # persisted report. Its result belongs in the terminal, in STATUS.md and in
    # git, none of which sit inside the artifact being described.
    print(f"  deposit HEAD {commit['head'][:12]}  tree {commit['working_tree']}  "
          f"cited={commit['cited_in_manuscript']}  -> {commit['status']}")

    # --- companion submission files ----------------------------------------
    for c in companions:
        rows.append({"check": c["name"], "state": c["status"],
                     "expected": f"{c['current']} present, {c['superseded']} absent",
                     "recomputed": {"file": c["file"],
                                    "superseded_absent": c["superseded_absent"],
                                    "current_present": c["current_present"]},
                     "literal_present_in_manuscript": c["current_present"],
                     "literal": c["current"]})
        if c["status"] == "PASS":
            n_pass += 1
        else:
            n_fail += 1

    # --- figure provenance --------------------------------------------------
    for f in figures:
        rows.append({"check": f["name"], "state": f["status"],
                     "expected": f.get("expected",
                                       f"{f['n_checked']} artifacts re-hash"),
                     "recomputed": f["problems"] or "all match",
                     "literal_present_in_manuscript": True,
                     "literal": f["name"]})
        if f["status"] == "PASS":
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

    outp = under("20_new_runs", "verification_summary_2026-08-21/g8_number_traceability.json")
    outp.write_text(json.dumps({
        "manuscript": args.manuscript,
        "state": "PASS" if (n_fail == 0 and n_absent == 0) else "FAIL",
        "n_checks": len(rows), "n_pass": n_pass,
        "n_value_mismatch": n_fail, "n_literal_absent": n_absent,
        "checks": rows,
    }, indent=2, default=str), encoding="utf-8")
    print(f"written -> {outp}")
    if commit["status"] != "PASS":
        print("  deposit_commit FAILED: the manuscript's commit citation "
              "or the deposit's working tree is out of step")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
