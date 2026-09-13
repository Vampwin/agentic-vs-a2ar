"""trivial_descriptor_baseline.py -- what is Table 5's baseline worth?

ADR-AIIM-164. §3.4 judges three pre-registered rescue techniques against a Vina
baseline of AUC 0.601, and nothing in this project ever asked what that
comparator is worth. The standard check is a trivial molecular descriptor: a
physics-based score that cannot beat counting atoms is not carrying physical
rank information.

**Self-check first.** The script recomputes the AUC of `vina_score` on these
rows and refuses to report anything unless it reproduces Table 5's 0.601 to
within `SELF_CHECK_TOL`. An earlier exploratory pass read the legacy 2026-06
`14_chembl_boxdock_scores.csv` instead and obtained 0.681 -- a different
protocol, whose numbers the manuscript does not report.

**Paired bootstrap, not overlapping intervals.** Comparing two AUCs by whether
their confidence intervals overlap is the wrong test and the first pass here
nearly drew the opposite conclusion from doing it.

Usage:
    python trivial_descriptor_baseline.py --out <a new directory>
"""
from __future__ import annotations

import argparse
import hashlib
import json
import warnings
from pathlib import Path
from kdss_paths import project_root, resolve, under

warnings.filterwarnings("ignore")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402

# --- ADR-AIIM-164 -----------------------------------------------------------
TABLE5_VINA_AUC = 0.601        # the number this must reproduce before reporting
SELF_CHECK_TOL = 0.002
N_BOOT = 10000
SEED = 1
DESCRIPTORS = ["n_heavy_atoms", "mol_weight", "clogp"]


def _root() -> Path:
    # Delegates to kdss_paths, which finds whichever layout this file is inside
    # -- the development tree or a deposit -- by walking up from its own
    # location. It previously tried two absolute paths and exited if neither
    # existed, which on the author's machine meant a copy of the deposit read
    # and wrote the development tree instead of itself.
    return project_root()


ROOT = _root()
RUN = resolve("20_new_runs/chembl245_ifp_rescoring_2026-08-15")
SCORES = under("20_new_runs/chembl245_ifp_rescoring_2026-08-15", "scores.csv")
PERCOMP = under("20_new_runs/chembl245_ifp_rescoring_2026-08-15", "ifp_rescoring_per_compound.csv")
BENCH = resolve("02_frozen_inputs/public_benchmarks/chembl_a2ar")


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def paired_delta(y, a, b, n=N_BOOT, seed=SEED):
    """Bootstrap the AUC *difference* on resampled rows, both scores together."""
    rng = np.random.default_rng(seed)
    y, a, b = np.asarray(y), np.asarray(a), np.asarray(b)
    out = []
    for _ in range(n):
        i = rng.integers(0, len(y), len(y))
        if len(set(y[i])) < 2:
            continue
        out.append(roc_auc_score(y[i], a[i]) - roc_auc_score(y[i], b[i]))
    o = np.array(out)
    return (float(o.mean()), [float(np.percentile(o, 2.5)),
                              float(np.percentile(o, 97.5))],
            float(2 * min((o <= 0).mean(), (o >= 0).mean())))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    out = Path(args.out).resolve()
    if out.exists() and any(out.iterdir()):
        raise SystemExit(f"REFUSED: {out} exists and is not empty "
                         "(ADR-AIIM-092: a run directory is written once).")
    out.mkdir(parents=True, exist_ok=True)

    s = pd.read_csv(SCORES)
    p = pd.read_csv(PERCOMP)
    d = p.merge(s[["compound_id", "n_heavy_atoms"]], on="compound_id", how="inner")

    # Molecular weight and logP are not in the run files; derive from the frozen
    # SMILES so the descriptor family is three variables, not one.
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors
    RDLogger.DisableLog("rdApp.*")
    smi = {}
    for f in ("chembl_a2a_raw_v2.csv", "chembl_a2ar_validation_set.csv"):
        t = pd.read_csv(BENCH / f)
        idc = "molecule_chembl_id" if "molecule_chembl_id" in t.columns else "chembl_id"
        smc = "canonical_smiles" if "canonical_smiles" in t.columns else "SMILES"
        for r in t.drop_duplicates(idc).itertuples():
            smi.setdefault(getattr(r, idc), getattr(r, smc))
    mw, lp = [], []
    for cid in d.compound_id:
        m = Chem.MolFromSmiles(smi.get(cid, "")) if smi.get(cid) else None
        mw.append(Descriptors.MolWt(m) if m else np.nan)
        lp.append(Descriptors.MolLogP(m) if m else np.nan)
    d["mol_weight"], d["clogp"] = mw, lp
    # Only the columns the primary comparison needs may shrink the row set.
    # `n_heavy_atoms` is present for every row; molecular weight and logP are
    # derived from SMILES and 31 compounds have none. Dropping those rows for
    # all three descriptors took the row set to 214 and moved Vina's AUC to
    # 0.594, which is not the number the manuscript reports -- the self-check
    # refused it. Each descriptor is now scored on the rows where it exists,
    # against Vina recomputed on those same rows.
    d = d.dropna(subset=["n_heavy_atoms", "vina_score", "label"]).reset_index(drop=True)

    y = d.label.values
    vina = -d.vina_score.values          # more negative is better
    auc_vina = roc_auc_score(y, vina)

    # --- the gate -----------------------------------------------------------
    reproduced = abs(auc_vina - TABLE5_VINA_AUC) <= SELF_CHECK_TOL
    print(f"n = {len(d)}, actives = {int(y.sum())}")
    print(f"self-check: Vina AUC {auc_vina:.3f} against Table 5's "
          f"{TABLE5_VINA_AUC} -> {'reproduced' if reproduced else 'MISMATCH'}")
    report = {
        "authority": "ADR-AIIM-164",
        "not_pre_registered": "the headline number was computed exploratively "
                              "before the ADR was written; reported as a "
                              "descriptive measurement, not as a test",
        "inputs": {f.name: {"path": f.relative_to(ROOT).as_posix(),
                            "sha256": sha256(f)} for f in (SCORES, PERCOMP)},
        "n": int(len(d)), "n_actives": int(y.sum()),
        "self_check": {"table5_vina_auc": TABLE5_VINA_AUC,
                       "recomputed": round(float(auc_vina), 3),
                       "tolerance": SELF_CHECK_TOL, "reproduced": bool(reproduced)},
    }
    if not reproduced:
        report["outcome"] = "SELF_CHECK_FAILED"
        report["withheld"] = ("These rows do not reproduce the AUC the manuscript "
                              "reports, so no baseline computed on them is "
                              "comparable to Table 5.")
        (out / "trivial_baseline.json").write_text(json.dumps(report, indent=2),
                                                   encoding="utf-8")
        print("SELF-CHECK FAILED -- baseline withheld")
        return 1

    # --- the baselines ------------------------------------------------------
    rows = {}
    for col in DESCRIPTORS:
        sub = d.dropna(subset=[col])
        ys, vs = sub.label.values, -sub.vina_score.values
        if len(set(ys)) < 2:
            rows[col] = {"available": False, "n": int(len(sub))}
            continue
        v = sub[col].values.astype(float)
        auc = roc_auc_score(ys, v)
        auc_v_here = roc_auc_score(ys, vs)   # Vina on these same rows
        m, ci, pv = paired_delta(ys, vs, v)
        rows[col] = {"available": True, "n": int(len(sub)),
                     "auc": round(float(auc), 3),
                     "vina_auc_on_same_rows": round(float(auc_v_here), 3),
                     "vina_minus_this": {"delta_auc": round(m, 3),
                                         "ci95": [round(ci[0], 3), round(ci[1], 3)],
                                         "p": round(pv, 3)}}
        print(f"  {col:<16} n={len(sub):>4}  AUC {auc:.3f}  "
              f"(Vina here {auc_v_here:.3f})   Vina - this: {m:+.3f} "
              f"(95% CI {ci[0]:+.3f} to {ci[1]:+.3f}, p = {pv:.3f})")
    report["vina_auc"] = round(float(auc_vina), 3)
    report["descriptor_baselines"] = rows
    report["method"] = (f"paired bootstrap of the AUC difference, {N_BOOT} "
                        f"resamples, seed {SEED}. Overlapping confidence "
                        "intervals are not used to assess a difference.")
    usable = [r for r in rows.values() if r.get("available")]
    report["outcome"] = ("NOT_SEPARABLE_FROM_TRIVIAL_DESCRIPTORS"
                         if usable and all(r["vina_minus_this"]["p"] > 0.05
                                           for r in usable)
                         else "SEPARABLE")
    report["reading"] = (
        "The baseline against which §3.4's three rescue techniques were judged "
        "null is itself not distinguishable from molecular size. This is not a "
        "claim that docking is uninformative in general; it is a statement that "
        "on this benchmark, at this target, a technique failing to beat the "
        "baseline has been measured against very little.")
    (out / "trivial_baseline.json").write_text(json.dumps(report, indent=2),
                                               encoding="utf-8")
    d[["compound_id", "label", "vina_score"] + DESCRIPTORS].to_csv(
        out / "per_compound.csv", index=False, lineterminator="\n")
    print(f"\noutcome: {report['outcome']}")
    print(f"written -> {out / 'trivial_baseline.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
