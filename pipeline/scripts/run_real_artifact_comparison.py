"""run_real_artifact_comparison.py -- snapshot testing vs contracts, on real runs.

Executes ADR-AIIM-207, frozen before any arm ran.

Two executions of the frozen pipeline on identical input already exist:
`snp_screen_repro_run1_2026-08-15` and `snp_screen_repro_run2_2026-08-15`.
Run 1 is treated as the stored reference and run 2 as a clean rerun. **A clean
rerun that raises an alarm is a false escalation**, which is the whole point:
no fault is injected anywhere in this script.

Four arms, all on the same pair of runs:

  1. byte-exact golden file over every artifact
  2. normalised golden file over every artifact -- the comparator that won on
     the synthetic benchmark (ADR-AIIM-202)
  3. byte-exact golden file scoped to the numeric outputs
  4. the contracts, read from the `contract_outcome` records both runs carry

Nothing is re-executed. Arm 4 reads what those runs recorded rather than
re-running today's evaluators against stored inputs, which would test the wrong
thing.

Output -> 20_new_runs/real_artifact_comparison_2026-09-22/
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kdss_paths import under  # noqa: E402

RUN1 = Path(under("20_new_runs", "snp_screen_repro_run1_2026-08-15"))
RUN2 = Path(under("20_new_runs", "snp_screen_repro_run2_2026-08-15"))

# The artifacts that carry the pipeline's numbers. Arm 3 is scoped to these.
NUMERIC_OUTPUTS = {
    "G3_redocking/ZMA_redocked_4EIY_energies.json",
    "G3_redocking/ZMA_redocked_4EIY.sdf",
    "G1_standardization/snp_source_ledger.csv",
}

# Fields that carry a wall-clock reading, stripped before comparing records.
TIME_FIELDS = ("timestamp", "record_id")


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def rel_files(root: Path) -> set[str]:
    return {str(p.relative_to(root)).replace("\\", "/")
            for p in root.rglob("*") if p.is_file()}


# --------------------------------------------------------------------------
# Arm 2: the normalising comparator, generalised past CSV
# --------------------------------------------------------------------------

def normalise(p: Path) -> str:
    """Strip the things a normalising comparator knows how to strip.

    Row order, float precision and key order are removed where the format
    allows. Nothing format-specific about PDB or PDBQT is known to it, because
    knowing that is the domain knowledge under test.
    """
    raw = p.read_bytes()
    name = p.name.lower()
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError:
        return hashlib.sha256(raw).hexdigest()

    if name.endswith(".json"):
        try:
            obj = json.loads(text)
            return json.dumps(_round(obj), sort_keys=True)
        except json.JSONDecodeError:
            pass
    if name.endswith(".jsonl"):
        rows = []
        for line in text.splitlines():
            if not line.strip():
                continue
            try:
                o = json.loads(line)
                if isinstance(o, dict):
                    o = {k: v for k, v in o.items() if k not in TIME_FIELDS}
                rows.append(json.dumps(_round(o), sort_keys=True))
            except json.JSONDecodeError:
                rows.append(line.strip())
        return "\n".join(sorted(rows))

    # Everything else: sort the lines and collapse float precision.
    lines = [_round_text(l.rstrip()) for l in text.splitlines() if l.strip()]
    return "\n".join(sorted(lines))


_FLOAT = re.compile(r"-?\d+\.\d+")


def _round_text(s: str) -> str:
    return _FLOAT.sub(lambda m: f"{float(m.group()):.3f}", s)


def _round(o):
    if isinstance(o, float):
        return round(o, 3)
    if isinstance(o, dict):
        return {k: _round(v) for k, v in o.items() if k not in TIME_FIELDS}
    if isinstance(o, list):
        return [_round(x) for x in o]
    return o


# --------------------------------------------------------------------------

def contract_outcomes(run: Path) -> list[dict]:
    reg = run / "artifact_registry.jsonl"
    out = []
    if not reg.exists():
        return out
    for line in reg.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        try:
            rec = json.loads(line)
        except json.JSONDecodeError:
            continue
        if rec.get("record_type") == "contract_outcome":
            # Strip wall-clock fields RECURSIVELY: state_history carries a
            # timestamp per transition, and a top-level strip leaves those in,
            # which makes two identical decisions look different.
            out.append(_round(rec))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(
        under("20_new_runs", "real_artifact_comparison_2026-09-22")))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    f1, f2 = rel_files(RUN1), rel_files(RUN2)
    common = sorted(f1 & f2)

    per_file = []
    for rel in common:
        a, b = RUN1 / rel, RUN2 / rel
        byte_same = sha256(a) == sha256(b)
        norm_same = byte_same or (normalise(a) == normalise(b))
        per_file.append({
            "file": rel,
            "byte_identical": byte_same,
            "identical_after_normalisation": norm_same,
            "in_numeric_scope": rel in NUMERIC_OUTPUTS,
        })

    arm1 = [r["file"] for r in per_file if not r["byte_identical"]]
    arm2 = [r["file"] for r in per_file
            if not r["identical_after_normalisation"]]
    arm3 = [r["file"] for r in per_file
            if r["in_numeric_scope"] and not r["byte_identical"]]

    c1, c2 = contract_outcomes(RUN1), contract_outcomes(RUN2)
    key = lambda recs: sorted(json.dumps(r, sort_keys=True) for r in recs)
    arm4_match = key(c1) == key(c2)

    arms = {
        "1_byte_exact_all_artifacts": {
            "files_compared": len(common), "files_flagged": len(arm1),
            "flagged": arm1,
            "false_escalation_on_clean_rerun": bool(arm1)},
        "2_normalised_all_artifacts": {
            "files_compared": len(common), "files_flagged": len(arm2),
            "flagged": arm2,
            "false_escalation_on_clean_rerun": bool(arm2)},
        "3_byte_exact_numeric_outputs_only": {
            "files_compared": sum(1 for r in per_file if r["in_numeric_scope"]),
            "files_flagged": len(arm3), "flagged": arm3,
            "false_escalation_on_clean_rerun": bool(arm3)},
        "4_contracts": {
            "contract_outcomes_run1": len(c1), "contract_outcomes_run2": len(c2),
            "outcomes_match_after_stripping_timestamps": arm4_match,
            "false_escalation_on_clean_rerun": not arm4_match},
    }

    summary = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "pre_registration": "ADR-AIIM-207",
        "nothing_injected": True,
        "nothing_re_executed": True,
        "runs": [str(RUN1.name), str(RUN2.name)],
        "n_files_common": len(common),
        "numeric_scope": sorted(NUMERIC_OUTPUTS),
        "arms": arms,
        "per_file": per_file,
        "scope_limit": "Two executions, one pipeline, three gates. A "
                       "demonstration, not a benchmark; no rate or interval is "
                       "computed from it.",
    }
    with open(out / "real_artifact_comparison.json", "w", encoding="utf-8",
              newline="\n") as f:
        json.dump(summary, f, indent=1, ensure_ascii=False)
        f.write("\n")

    print("=" * 74)
    print("SNAPSHOT TESTING vs CONTRACTS ON TWO REAL RUNS (ADR-AIIM-207)")
    print("=" * 74)
    print(f"  {len(common)} artifacts common to both runs; no fault injected\n")
    print(f"  {'arm':44} {'flagged':>9} {'false escalation':>18}")
    for name, a in arms.items():
        flagged = a.get("files_flagged")
        cell = "n/a" if flagged is None else f"{flagged}/{a.get('files_compared')}"
        print(f"  {name:44} {cell:>9} "
              f"{str(a['false_escalation_on_clean_rerun']):>18}")
    print()
    for r in per_file:
        if not r["byte_identical"]:
            n = "normalised away" if r["identical_after_normalisation"] else "SURVIVES normalisation"
            print(f"    differs: {r['file']:48} {n}")
    print()
    print(f"written -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
