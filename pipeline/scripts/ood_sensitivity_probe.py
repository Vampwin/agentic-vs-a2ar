"""ood_sensitivity_probe.py -- does the applicability-domain detector detect anything?

The frozen detector `ood_isoforest_dev_2026-08-02` (ADR-AIIM-031) is load-bearing:
`risk_monitor/inference.py` fails closed when `out_of_domain_flag` is True or
None, so an unresolved domain status halts inference exactly as a leakage-guard
rejection does.

Everything recorded about it is in-distribution. The training self-check flagged
1 of 188 rows (0.53 %) and the 2x2 ablation flagged 2 of 376 evaluations
(0.53 %) -- and `contamination` was set to 0.01, which *defines* the expected
in-training anomaly fraction. Both numbers are that setting reproducing itself.
Nothing establishes that the detector flags telemetry which really is out of
domain, which is the only thing it exists to do.

So: displace the development rows away from the training distribution by a
factor k and report the flag rate against k. k = 1 is the negative control and
must reproduce the recorded rate, or nothing else here is reported.

Criteria and the reading of each outcome are fixed in ADR-AIIM-171, written
before this script was run.

**No re-fitting.** The registered artifact is loaded and its SHA-256 checked
against `ood_fit_summary.json` before use.

**Telemetry only.** The leakage guard runs over the perturbed columns before
anything is scored, and the run aborts if it objects.

Usage:
    python ood_sensitivity_probe.py --out 20_new_runs/ood_sensitivity_2026-09-05
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from kdss_paths import project_root, resolve, under, dev_dir
from types import SimpleNamespace

import numpy as np

# --- pre-registered, ADR-AIIM-171 -------------------------------------------
K_FACTORS = [1, 2, 5, 10, 100, 1000]
RECORDED_TRAIN_RATE = 1 / 188          # the training self-check, 0.005319...
NEG_CONTROL_TOL = 0.01                 # absolute, on a rate
ARTIFACT_ID = "ood_isoforest_dev_2026-08-02"


def _root() -> Path:
    # Delegates to kdss_paths, which finds whichever layout this file is inside
    # -- the development tree or a deposit -- by walking up from its own
    # location. It previously tried two absolute paths and exited if neither
    # existed, which on the author's machine meant a copy of the deposit read
    # and wrote the development tree instead of itself.
    return project_root()


ROOT = _root()
SYS = dev_dir("10_new_system")
READY = resolve("20_new_runs/WP2.5_WP3A_readiness_2026-07-30")
PILOT_JSON = READY / "pilot_postfix_2026-08-04.json"
OOD_SUMMARY = READY / "ood_fit_summary.json"
MODEL = READY / "models" / f"{ARTIFACT_ID}.joblib"


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    out = Path(args.out).resolve()
    if out.exists() and any(out.iterdir()):
        raise SystemExit(f"REFUSED: {out} exists and is not empty (ADR-AIIM-092).")
    out.mkdir(parents=True, exist_ok=True)

    sys.path.insert(0, str(SYS))
    from fault_injection.pilot_runner import build_training_dataset
    from telemetry.ood_detector import OutOfDomainDetector
    from telemetry.leakage_guard import check_field_names, check_allowlist

    # --- the registered artifact, verified before use -----------------------
    summary = json.loads(OOD_SUMMARY.read_text(encoding="utf-8"))
    want = summary["artifact"]["sha256"].lower()
    got = sha256(MODEL).lower()
    if want != got:
        raise SystemExit(f"REFUSED: {ARTIFACT_ID} hash mismatch\n"
                         f"  registered {want}\n  on disk    {got}")
    det = OutOfDomainDetector.load(MODEL)

    # --- the same rows the model was fitted on ------------------------------
    raw = json.loads(PILOT_JSON.read_text(encoding="utf-8"))
    X, y, groups = build_training_dataset(SimpleNamespace(outcomes=raw["outcomes"]),
                                          reference_condition="B")
    check_field_names(X.columns, context="ood_sensitivity_probe X columns")
    check_allowlist(X.columns, context="ood_sensitivity_probe X columns")

    used = [c for c in det.non_null_features if c in X.columns]
    numeric = [c for c in used if np.issubdtype(X[c].dtype, np.number)]

    # --- the probe -----------------------------------------------------------
    rows = []
    for k in K_FACTORS:
        Xk = X.copy()
        for c in numeric:
            Xk[c] = Xk[c] * k
        check_allowlist(Xk.columns, context=f"perturbed columns k={k}")
        flags = [bool(det.predict_one(r.to_dict())) for _, r in Xk.iterrows()]
        n = sum(flags)
        rows.append({"k": k, "n_flagged": n, "n_total": len(flags),
                     "flag_rate": round(n / len(flags), 5)})
        print(f"  k={k:<5} flagged {n:>3}/{len(flags)}  = {n/len(flags):.3%}")

    neg = rows[0]
    ok = abs(neg["flag_rate"] - RECORDED_TRAIN_RATE) <= NEG_CONTROL_TOL
    crossing = next((r["k"] for r in rows if r["flag_rate"] >= 0.5), None)

    report = {
        "question": "Does the frozen applicability-domain detector flag telemetry "
                    "that is genuinely outside its training distribution? Every "
                    "number recorded for it until now was in-distribution.",
        "authority": "ADR-AIIM-171, criteria fixed before the run",
        "artifact": {"id": ARTIFACT_ID, "sha256_verified": True,
                     "refitted": False,
                     "contamination": summary["hyperparameters"]["contamination"],
                     "n_features_declared": summary["dataset"]["n_features_total"],
                     "n_features_non_null": summary["dataset"]["n_features_non_null"],
                     "n_features_perturbed": len(numeric)},
        "method": "Every non-null numeric telemetry feature multiplied by k. "
                  "k=1 is the negative control. No re-fitting; the registered "
                  "artifact is loaded and hash-checked. Leakage guard run over "
                  "the perturbed columns before scoring.",
        "negative_control": {"k": 1, "observed_rate": neg["flag_rate"],
                             "recorded_training_rate": round(RECORDED_TRAIN_RATE, 5),
                             "tolerance": NEG_CONTROL_TOL, "passed": bool(ok)},
        "by_k": rows,
        "k_at_50pc": crossing,
        "declared_reading": {
            "approaches_100pc_by_k10": "the detector does what its fail-closed "
                                       "role requires",
            "stays_at_contamination_floor": "it cannot distinguish telemetry a "
                                            "thousandfold outside its training "
                                            "range from telemetry inside it; "
                                            "section 2.1's description must be "
                                            "narrowed to what it demonstrably does",
        },
    }
    if not ok:
        report["outcome"] = "NEGATIVE_CONTROL_FAILED"
        report["results_withheld"] = (
            f"At k=1 the detector flagged {neg['flag_rate']:.5f} against a "
            f"recorded {RECORDED_TRAIN_RATE:.5f}. The artifact loaded does not "
            "behave as the one registered; no other number here is reported.")
        print(f"\nNEGATIVE CONTROL FAILED -- results withheld")
    else:
        last = rows[-1]
        report["outcome"] = ("DETECTS" if crossing is not None and crossing <= 10
                             else "PARTIAL" if crossing is not None
                             else "DOES_NOT_DETECT")
        print(f"\nnegative control {neg['flag_rate']:.3%} -- passed.  "
              f"at k={last['k']} flag rate {last['flag_rate']:.3%}  "
              f"-> {report['outcome']}")
    (out / "ood_sensitivity.json").write_text(json.dumps(report, indent=2),
                                              encoding="utf-8")
    print(f"\nwritten -> {out / 'ood_sensitivity.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
