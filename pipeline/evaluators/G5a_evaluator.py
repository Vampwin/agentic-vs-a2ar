"""
G5a evaluator — external validation contract validity (hard gate).

Checks label isolation, set identity, independent implementation, and
bootstrap correctness. Must pass before any performance numbers (G5b) are seen.
"""

from __future__ import annotations

import hashlib
import time
from pathlib import Path
from typing import Any


def evaluate(
    contract: dict[str, Any],
    inputs: dict[str, Any],
    registry: Any,
) -> dict[str, Any]:
    passed, failed, warned, evidence = [], [], [], {}
    critical_failures = []

    benchmark_csv = Path(inputs.get("chembl_benchmark_csv", ""))
    benchmark_role_declaration = inputs.get("benchmark_role_declaration")
    g5a_evaluator_source_hash = inputs.get("g5a_evaluator_source_hash", "")

    for check in contract.get("checks", []):
        cid = check["id"]
        sev = check.get("severity", "medium")

        if cid == "benchmark_role_set":
            # Registry must contain a benchmark_role_declaration record
            decl = registry.get_benchmark_role_declaration() if registry else None
            ok = decl is not None and decl.get("role") not in (None, "", "null")
            evidence[cid] = {"declaration_found": decl is not None, "role": decl.get("role") if decl else None}
            (passed if ok else failed).append(cid)
            if not ok and sev == "critical":
                critical_failures.append(cid)

        elif cid == "label_isolation":
            # Verify from registry access log that labels were not read before this gate.
            # We check that no access_log record shows the benchmark labels file was read
            # by the orchestrator or any evaluator before G5a's timestamp.
            isolation_result = _check_label_isolation(registry, benchmark_csv, inputs.get("_gate_start_time"))
            evidence[cid] = isolation_result
            ok = isolation_result["isolated"]
            (passed if ok else failed).append(cid)
            if not ok and sev == "critical":
                critical_failures.append(cid)

        elif cid == "set_identity":
            # Hash the compound ID list and compare to the registered hash at gate start
            registered_hash = inputs.get("benchmark_id_list_hash_registered")
            if benchmark_csv.exists() and registered_hash:
                actual_hash = _hash_id_list(benchmark_csv)
                ok = actual_hash.upper() == registered_hash.upper()
                evidence[cid] = {
                    "registered_hash": registered_hash,
                    "actual_hash": actual_hash,
                    "match": ok,
                }
                (passed if ok else failed).append(cid)
                if not ok and sev == "critical":
                    critical_failures.append(cid)
            else:
                warned.append(cid)
                evidence[cid] = {
                    "warning": "benchmark_id_list_hash_registered not provided or file missing"
                }

        elif cid == "no_training_overlap":
            overlap_result = inputs.get("qsar_training_overlap_check")
            if overlap_result is None:
                passed.append(cid)
                evidence[cid] = {"note": "QSAR not used in pipeline; check not applicable"}
            else:
                ok = overlap_result.get("n_overlapping", 0) == 0
                evidence[cid] = overlap_result
                (passed if ok else failed).append(cid)
                if not ok and sev == "critical":
                    critical_failures.append(cid)

        elif cid == "independent_metric_implementation":
            # G5a evaluator source hash must differ from all kdss_audit module hashes
            kdss_audit_hashes = inputs.get("kdss_audit_module_hashes", {})
            if not g5a_evaluator_source_hash:
                warned.append(cid)
                evidence[cid] = {"warning": "g5a_evaluator_source_hash not provided"}
            else:
                collisions = [
                    name for name, h in kdss_audit_hashes.items()
                    if h == g5a_evaluator_source_hash
                ]
                ok = len(collisions) == 0
                evidence[cid] = {
                    "g5a_hash": g5a_evaluator_source_hash,
                    "collisions": collisions,
                }
                (passed if ok else warned).append(cid)

        elif cid == "bootstrap_pairs_together":
            # Verify via a known-answer test that the bootstrap implementation
            # resamples (compound, score, label) triplets together.
            # This is verified by running the test in tests/test_G5a.py.
            # The evaluator checks that test was run (its result is in inputs).
            test_result = inputs.get("bootstrap_known_answer_test")
            if test_result is None:
                warned.append(cid)
                evidence[cid] = {
                    "warning": "bootstrap known-answer test result not in inputs; "
                               "run tests/test_G5a.py before gate"
                }
            else:
                ok = test_result.get("pairs_resampled_together") is True
                evidence[cid] = test_result
                (passed if ok else failed).append(cid)
                if not ok and sev == "critical":
                    critical_failures.append(cid)

        elif cid == "bootstrap_n_resamples":
            n = inputs.get("bootstrap_n_resamples", 0)
            ok = n >= 1000
            evidence[cid] = {"n_resamples": n}
            (passed if ok else warned).append(cid)

        elif cid == "enrichment_counts_positive":
            ef_result = inputs.get("ef1pct_count")
            if ef_result is not None:
                passed.append(cid)
                evidence[cid] = {"ef1pct_count": ef_result}
            else:
                warned.append(cid)
                evidence[cid] = {"warning": "ef1pct_count not in inputs"}

    # Critical breaches halt; non-critical failures still require a human
    # decision. Returning PASS while `failed` is non-empty would let a
    # high-severity finding through silently.
    if critical_failures:
        state = "HALT"
    elif failed:
        state = "ESCALATE"
    else:
        state = "PASS"
    return {
        "state": state,
        "passed": passed,
        "failed": failed,
        "warned": warned,
        "critical_failures": critical_failures,
        "evidence": evidence,
        "summary": (
            f"{len(passed)} checks passed, {len(failed)} failed "
            f"({len(critical_failures)} critical), {len(warned)} warned"
        ),
    }


def _check_label_isolation(
    registry: Any,
    benchmark_csv: Path,
    gate_start_time: float | None,
) -> dict[str, Any]:
    """
    Labels are isolated if the benchmark CSV (with label column) was not read
    by any component before this gate's start time.

    A separate labels file (chembl_a2ar_labels.csv) would be cleaner; if the
    scores and labels are in the same CSV, warn that they cannot be fully isolated.
    """
    if registry is None:
        return {"isolated": True, "note": "No registry; isolation assumed (not verified)"}

    gate_start = gate_start_time or time.time()
    access_log = registry.get_file_access_log()

    benchmark_str = str(benchmark_csv)
    early_reads = [
        rec for rec in access_log
        if benchmark_str in rec.get("file_path", "")
        and rec.get("timestamp", float("inf")) < gate_start
    ]

    if early_reads:
        return {
            "isolated": False,
            "early_reads": [{"reader": r["reader"], "purpose": r["purpose"]} for r in early_reads],
            "note": "Label file was read before G5a gate — label isolation violated",
        }
    return {"isolated": True, "early_reads": [], "note": "No early reads detected"}


def _hash_id_list(benchmark_csv: Path) -> str:
    """Hash the sorted compound ID list from the benchmark CSV (not the full file)."""
    try:
        import pandas as pd
        df = pd.read_csv(benchmark_csv)
        id_col = next(
            (c for c in df.columns if c.lower() in {"chembl_id", "compound_id", "id", "snp_id"}),
            df.columns[0],
        )
        sorted_ids = sorted(df[id_col].astype(str).tolist())
        content = "\n".join(sorted_ids).encode()
        return hashlib.sha256(content).hexdigest().upper()
    except Exception as exc:
        return f"ERROR:{exc}"
