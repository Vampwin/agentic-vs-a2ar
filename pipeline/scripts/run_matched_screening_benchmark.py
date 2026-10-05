"""run_matched_screening_benchmark.py -- screening power with the size confound removed.

Executes ADR-AIIM-204, frozen before any metric here was computed.

§3.4 reports Vina at ROC-AUC 0.601 against 0.584 for heavy-atom count over the
245-compound ChEMBL benchmark, a paired difference of +0.018 at p = 0.722. That
comparison is uninterpretable for a stated reason: the classes are separated by
5.24 log units of median pChEMBL and differ in size, so the benchmark measures
discrimination at wide separation and partly measures molecular weight.

This removes the size confound by matching each active to a physicochemically
similar inactive, and recomputes the same metrics on the matched subset.

Nothing is docked. Every Vina score is read from the frozen artifact
`20_new_runs/chembl245_ifp_rescoring_2026-08-15/ifp_rescoring_per_compound.csv`.

Output -> 20_new_runs/matched_screening_2026-09-21/
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kdss_paths import under  # noqa: E402

from rdkit import Chem, RDLogger  # noqa: E402
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors  # noqa: E402

RDLogger.DisableLog("rdApp.*")

BOOTSTRAP_N = 1000
BOOTSTRAP_SEED = 42
BEDROC_ALPHA = 20.0
EF_FRACTIONS = (0.01, 0.05, 0.10)
CALIPER_SD = 1.0
MATCH_VARS = ("MW", "logP", "HBD", "HBA", "rotB", "charge")

BENCH = under("20_new_runs", "chembl245_ifp_rescoring_2026-08-15/"
                             "ifp_rescoring_per_compound.csv")
RAW = under("02_frozen_inputs/public_benchmarks/chembl_a2ar", "chembl_a2a_raw_v2.csv")


# --------------------------------------------------------------------------
# Metrics. Higher score must mean "more likely active" for every scorer.
# --------------------------------------------------------------------------

def roc_auc(y: np.ndarray, s: np.ndarray) -> float:
    pos, neg = s[y == 1], s[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    order = np.argsort(s, kind="mergesort")
    ranks = np.empty(len(s), float)
    ranks[order] = np.arange(1, len(s) + 1)
    # average ranks for ties
    df = pd.DataFrame({"s": s, "r": ranks})
    ranks = df.groupby("s")["r"].transform("mean").to_numpy()
    return (ranks[y == 1].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def bedroc(y: np.ndarray, s: np.ndarray, alpha: float = BEDROC_ALPHA) -> float:
    n = len(y)
    na = int(y.sum())
    if na == 0 or na == n:
        return float("nan")
    order = np.argsort(-s, kind="mergesort")
    ranks = np.where(y[order] == 1)[0] + 1
    ra = na / n
    rie = (np.exp(-alpha * ranks / n).sum() / na) * (
        (1 - np.exp(-alpha)) / (alpha / n * (1 - np.exp(-alpha / n) * 0 + 0) + 1e-300)
    ) if False else (np.exp(-alpha * ranks / n).sum()
                     / (na * (1 - np.exp(-alpha)) / (n * (np.exp(alpha / n) - 1))))
    rie_min = ((1 - np.exp(alpha * ra)) / (ra * (1 - np.exp(alpha)))) * 0 + \
              (1 - np.exp(alpha * ra)) / (ra * (1 - np.exp(alpha)))
    rie_max = (1 - np.exp(-alpha * ra)) / (ra * (1 - np.exp(-alpha)))
    if abs(rie_max - rie_min) < 1e-12:
        return float("nan")
    return float((rie - rie_min) / (rie_max - rie_min))


def enrichment(y: np.ndarray, s: np.ndarray, frac: float) -> float:
    n = len(y)
    # ceil, not round: this reproduces the published EF convention of the
    # 245-compound benchmark (EF1 = 1.6333 at k = ceil(2.45) = 3).
    k = max(1, int(np.ceil(n * frac)))
    order = np.argsort(-s, kind="mergesort")
    top = y[order][:k]
    ra = y.mean()
    return float(top.mean() / ra) if ra > 0 else float("nan")


def boot_ci(y, s, fn, n_boot=BOOTSTRAP_N, seed=BOOTSTRAP_SEED):
    rng = np.random.default_rng(seed)
    vals = []
    idx = np.arange(len(y))
    for _ in range(n_boot):
        b = rng.choice(idx, size=len(idx), replace=True)
        if len(np.unique(y[b])) < 2:
            continue
        vals.append(fn(y[b], s[b]))
    if not vals:
        return (float("nan"), float("nan"))
    return (float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5)))


def paired_delta_ci(y, s1, s2, n_boot=BOOTSTRAP_N, seed=BOOTSTRAP_SEED):
    """Bootstrap interval for AUC(s1) - AUC(s2) on the same rows."""
    rng = np.random.default_rng(seed)
    idx = np.arange(len(y))
    d = []
    for _ in range(n_boot):
        b = rng.choice(idx, size=len(idx), replace=True)
        if len(np.unique(y[b])) < 2:
            continue
        d.append(roc_auc(y[b], s1[b]) - roc_auc(y[b], s2[b]))
    d = np.array(d)
    return (float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5)),
            float((d <= 0).mean()))


# --------------------------------------------------------------------------

def load(bench: str | Path = BENCH) -> tuple[pd.DataFrame, dict]:
    """Read a 245-compound score file and attach matching descriptors.

    `bench` is a parameter so the identical matching rule can be applied to the
    single-conformer re-run (ADR-AIIM-227) as well as to the eight-pooled scores
    this analysis was first run on. The descriptors come from SMILES and are
    protocol-independent; only the `vina` column differs between the two.
    """
    b = pd.read_csv(bench)
    raw = pd.read_csv(RAW)
    smi = raw.drop_duplicates("molecule_chembl_id").set_index(
        "molecule_chembl_id")["canonical_smiles"]
    pch = raw.groupby("molecule_chembl_id")["pchembl_value"].median()
    b["smiles"] = b.compound_id.map(smi)
    b["pchembl"] = b.compound_id.map(pch)

    dropped = b[b.smiles.isna()]
    rows = []
    for _, r in b.dropna(subset=["smiles"]).iterrows():
        m = Chem.MolFromSmiles(r.smiles)
        if m is None:
            continue
        rows.append({
            "compound_id": r.compound_id, "label": int(r.label),
            "vina": float(r.vina_score), "pchembl": r.pchembl,
            "MW": Descriptors.MolWt(m), "logP": Crippen.MolLogP(m),
            "HBD": rdMolDescriptors.CalcNumHBD(m),
            "HBA": rdMolDescriptors.CalcNumHBA(m),
            "rotB": rdMolDescriptors.CalcNumRotatableBonds(m),
            "charge": Chem.GetFormalCharge(m),
            "heavy": m.GetNumHeavyAtoms(),
        })
    d = pd.DataFrame(rows)
    loss = {
        "n_benchmark": int(len(b)),
        "n_entering": int(len(d)),
        "n_excluded_no_smiles": int(len(dropped)),
        "excluded_label_counts": dropped.label.value_counts().to_dict(),
    }
    return d, loss


def match(d: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Greedy 1:1 nearest-neighbour matching, caliper 1.0 SD. Deterministic."""
    z = d[list(MATCH_VARS)].astype(float)
    z = (z - z.mean()) / z.std(ddof=0).replace(0, 1.0)
    A = d.index[d.label == 1].to_numpy()
    I = d.index[d.label == 0].to_numpy()
    za, zi = z.loc[A].to_numpy(), z.loc[I].to_numpy()

    dist = np.linalg.norm(za[:, None, :] - zi[None, :, :], axis=2)
    # tightest pairs first, so the rule does not depend on row order
    order = np.argsort(dist.min(axis=1), kind="mergesort")
    used, pairs = set(), []
    for ai in order:
        cand = np.argsort(dist[ai], kind="mergesort")
        for ci in cand:
            if ci in used:
                continue
            if dist[ai, ci] > CALIPER_SD:
                break
            used.add(ci)
            pairs.append((A[ai], I[ci], float(dist[ai, ci])))
            break
    keep = [p[0] for p in pairs] + [p[1] for p in pairs]
    m = d.loc[keep].copy()
    info = {
        "n_actives_available": int((d.label == 1).sum()),
        "n_inactives_available": int((d.label == 0).sum()),
        "n_pairs_formed": len(pairs),
        "n_actives_unmatched": int((d.label == 1).sum()) - len(pairs),
        "caliper_sd": CALIPER_SD,
        "mean_pair_distance": round(float(np.mean([p[2] for p in pairs])), 4) if pairs else None,
        "max_pair_distance": round(float(np.max([p[2] for p in pairs])), 4) if pairs else None,
    }
    return m, info


def score_block(y, s, name) -> dict:
    lo, hi = boot_ci(y, s, roc_auc)
    return {
        "scorer": name,
        "auc_roc": round(roc_auc(y, s), 4),
        "auc_ci95": [round(lo, 4), round(hi, 4)],
        "bedroc20": round(bedroc(y, s), 4),
        **{f"EF{int(f*100)}": round(enrichment(y, s, f), 3) for f in EF_FRACTIONS},
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(under("20_new_runs",
                                               "matched_screening_2026-09-21")))
    ap.add_argument("--bench", default=str(BENCH),
                    help="245-compound score file; defaults to the eight-pooled "
                         "scores this analysis was first run on")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    d, loss = load(args.bench)
    m, info = match(d)
    m.to_csv(out / "matched_pairs.csv", index=False)

    res = {}
    for tag, frame in (("unmatched_all", d), ("matched", m)):
        y = frame.label.to_numpy()
        # Vina energies are negative-is-better, so negate for "higher = active".
        res[tag] = {
            "n": int(len(frame)), "n_actives": int(y.sum()),
            "n_inactives": int((y == 0).sum()),
            "median_pchembl_actives": round(float(frame[frame.label == 1].pchembl.median()), 2),
            "median_pchembl_inactives": round(float(frame[frame.label == 0].pchembl.median()), 2),
            "mean_heavy_actives": round(float(frame[frame.label == 1].heavy.mean()), 2),
            "mean_heavy_inactives": round(float(frame[frame.label == 0].heavy.mean()), 2),
            "vina": score_block(y, -frame.vina.to_numpy(), "vina (negated energy)"),
            "heavy_atoms": score_block(y, frame.heavy.to_numpy().astype(float),
                                       "heavy-atom count"),
        }
        lo, hi, p_le0 = paired_delta_ci(y, -frame.vina.to_numpy(),
                                        frame.heavy.to_numpy().astype(float))
        res[tag]["paired_delta_auc_vina_minus_heavy"] = {
            "point": round(res[tag]["vina"]["auc_roc"]
                           - res[tag]["heavy_atoms"]["auc_roc"], 4),
            "ci95": [round(lo, 4), round(hi, 4)],
            "bootstrap_fraction_le_zero": round(p_le0, 4),
            "includes_zero": bool(lo <= 0 <= hi),
        }

    summary = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "pre_registration": "ADR-AIIM-204",
        "nothing_docked": True,
        "source_scores": "20_new_runs/chembl245_ifp_rescoring_2026-08-15/"
                         "ifp_rescoring_per_compound.csv",
        "design": "greedy 1:1 nearest-neighbour matching on the six DUD-E "
                  "properties, standardised, caliper 1.0 SD; heavy-atom count "
                  "deliberately excluded from the matching variables",
        "match_vars": list(MATCH_VARS),
        "bootstrap": {"n": BOOTSTRAP_N, "seed": BOOTSTRAP_SEED},
        "population": loss,
        "matching": info,
        "results": res,
        "prior_published": {"auc_vina": 0.601, "auc_heavy": 0.584,
                            "source": "§3.4, unmatched n=245"},
    }
    with open(out / "matched_screening_summary.json", "w", encoding="utf-8",
              newline="\n") as f:
        json.dump(summary, f, indent=1, ensure_ascii=False)
        f.write("\n")

    print("=" * 74)
    print("SCREENING POWER WITH THE SIZE CONFOUND REMOVED (ADR-AIIM-204)")
    print("=" * 74)
    print(f"entering {loss['n_entering']} of {loss['n_benchmark']}; "
          f"{loss['n_excluded_no_smiles']} lacked a SMILES "
          f"{loss['excluded_label_counts']}")
    print(f"pairs formed {info['n_pairs_formed']} of {info['n_actives_available']} "
          f"actives; mean distance {info['mean_pair_distance']} SD, "
          f"max {info['max_pair_distance']}")
    print()
    for tag in ("unmatched_all", "matched"):
        r = res[tag]
        print(f"--- {tag}  n={r['n']} ({r['n_actives']} active / "
              f"{r['n_inactives']} inactive) ---")
        print(f"    heavy atoms  actives {r['mean_heavy_actives']} vs "
              f"inactives {r['mean_heavy_inactives']}")
        for k in ("vina", "heavy_atoms"):
            b = r[k]
            print(f"    {b['scorer']:26} AUC {b['auc_roc']:.4f} "
                  f"[{b['auc_ci95'][0]:.3f}, {b['auc_ci95'][1]:.3f}]  "
                  f"BEDROC20 {b['bedroc20']:.3f}  EF1 {b['EF1']}")
        p = r["paired_delta_auc_vina_minus_heavy"]
        print(f"    paired dAUC (vina - heavy) {p['point']:+.4f} "
              f"[{p['ci95'][0]:+.3f}, {p['ci95'][1]:+.3f}]  "
              f"includes zero: {p['includes_zero']}")
        print()
    print(f"written -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
