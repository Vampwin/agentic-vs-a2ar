"""The four added checks. Each returns (ok, detail) and fails closed: a missing record or an unreadable file is a failure, not a pass."""
from __future__ import annotations

import re
from pathlib import Path

COVERAGE_MIN = 0.95          # fixed in the pre-registration, before any faulted workspace was evaluated
SITE_RADIUS = 10.0


def _xyz(ln):
    return (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))


# ---------------------------------------------------------------- G2: protonation_matches_record
def his_states(pdbqt: Path) -> dict[str, dict]:
    res: dict[str, dict] = {}
    for ln in Path(pdbqt).read_text(errors="replace").splitlines():
        if ln.startswith("ATOM  ") and ln[17:20] in ("HIS", "HID", "HIE", "HIP"):
            r = res.setdefault(ln[21:27], {"atoms": {}, "types": {}})
            r["atoms"][ln[12:16].strip()] = 1
            r["types"][ln[12:16].strip()] = ln[77:79].strip()
    out = {}
    for k, r in res.items():
        a, t = r["atoms"], r["types"]
        hd1, he2 = "HD1" in a, "HE2" in a
        state = "HIP" if hd1 and he2 else "HID" if hd1 else "HIE" if he2 else "HIS_NO_RING_H"
        typed_ok = (state == "HID" and t.get("ND1") == "N" and t.get("NE2") == "NA") or \
                   (state == "HIE" and t.get("NE2") == "N" and t.get("ND1") == "NA") or \
                   (state == "HIP" and t.get("ND1") == "N" and t.get("NE2") == "N")
        out[k] = {"state": state, "types_consistent": bool(typed_ok)}
    return out


def his_type_inconsistencies(receptor) -> dict[str, str]:
    """Histidines whose ring-nitrogen AutoDock types disagree with the hydrogen they carry. An observation reported over the project's receptors; it
    is not one of the added gate checks (the F06 injector changes the hydrogen and the types together, so this does not detect it)."""
    return {k.strip(): v["state"] for k, v in his_states(Path(receptor)).items() if not v["types_consistent"]}


def protonation_matches_record(receptor, decision):
    if not receptor or not Path(str(receptor)).is_file():
        return False, {"error": "prepared receptor not found"}
    if not decision:
        return False, {"error": "no protonation decision recorded"}
    now = his_states(Path(receptor))
    bad = {}
    for k, v in now.items():
        rec = decision.get(k)
        if rec is None:
            bad[k.strip()] = {"file": v["state"], "recorded": None}
        elif rec != v["state"]:
            bad[k.strip()] = {"file": v["state"], "recorded": rec}
    for k in decision:
        if k not in now:
            bad[k.strip()] = {"file": None, "recorded": "present in the record only"}
    return (not bad), {"n_histidines": len(now), "mismatches": bad}


# ---------------------------------------------------------------- G7: trajectory_covers_plan, binding_site_matches_receptor
def trajectory_covers_plan(eq_csv, planned_ns: float):
    import pandas as pd
    if not eq_csv or not Path(str(eq_csv)).is_file():
        return False, {"error": "no per-frame series"}
    df = pd.read_csv(eq_csv)
    last = df.groupby(["compound_id", "replicate"])["time_ns"].max()
    short = {f"{k[0]}/{k[1]}": round(float(v), 3) for k, v in last.items() if v < COVERAGE_MIN * planned_ns}
    return (not short), {"planned_ns": planned_ns, "coverage_min": COVERAGE_MIN, "n_series": int(len(last)), "too_short": short}


def site_from_receptor(receptor, ligand_sdf, radius: float = SITE_RADIUS) -> list[int]:
    import numpy as np
    from rdkit import Chem
    mol = Chem.MolFromMolFile(str(ligand_sdf), removeHs=False)
    lig = mol.GetConformer().GetPositions()
    out = []
    for ln in Path(receptor).read_text(errors="replace").splitlines():
        if ln.startswith("ATOM  ") and ln[12:16].strip() == "CA":
            if np.min(np.linalg.norm(lig - np.array(_xyz(ln)), axis=1)) <= radius:
                out.append(int(ln[22:26]))
    return sorted(set(out))


def binding_site_matches_receptor(results_csv, receptor, ligand_sdf):
    import pandas as pd
    if not (results_csv and receptor and ligand_sdf) or not all(Path(str(p)).is_file() for p in (results_csv, receptor, ligand_sdf)):
        return False, {"error": "results table, receptor or reference ligand not found"}
    want = set(site_from_receptor(receptor, ligand_sdf))
    df = pd.read_csv(results_csv)
    if "binding_site_residues" not in df.columns:
        return False, {"error": "no binding_site_residues column"}
    bad = {}
    for _, r in df.iterrows():
        got = {int(x) for x in str(r["binding_site_residues"]).split()}
        if got != want:
            bad[f"{r['compound_id']}/{r['replicate']}"] = {"n_missing": len(want - got), "n_extra": len(got - want)}
    return (not bad), {"n_site_residues_expected": len(want), "mismatching_rows": bad}


# ---------------------------------------------------------------- G5a: benchmark_ids_not_exposed
_SPLIT = re.compile(r"[\s,;\"']+")


def benchmark_ids_not_exposed(registry, benchmark_csv, registered_inputs, gate_start):
    import time

    import pandas as pd
    if registry is None or not benchmark_csv or not Path(str(benchmark_csv)).is_file():
        return False, {"error": "no registry or benchmark file"}
    ids = set(pd.read_csv(benchmark_csv).iloc[:, 0].astype(str))
    allowed = {str(Path(p)) for p in registered_inputs}
    start = gate_start or time.time()
    flagged = {}
    for rec in registry.get_file_access_log():
        if rec.get("reader") in ("G5a", "G5b") or rec.get("timestamp", float("inf")) >= start:
            continue
        p = rec.get("file_path", "")
        if str(Path(p)) in allowed:
            continue
        if not Path(p).is_file():
            flagged[p] = {"error": "file read before the gate cannot be opened"}
            continue
        hit = ids & set(_SPLIT.split(Path(p).read_text(errors="replace")))
        if hit:
            flagged[p] = {"reader": rec.get("reader"), "n_benchmark_ids": len(hit)}
    return (not flagged), {"n_benchmark_ids": len(ids), "flagged_files": flagged}


def merge(frozen: dict, new: list) -> dict:
    """Append the added checks to the frozen evaluator's output; the worse state wins."""
    passed, failed = list(frozen.get("passed", [])), list(frozen.get("failed", []))
    crit, warned = list(frozen.get("critical_failures", [])), list(frozen.get("warned", []))
    ev = dict(frozen.get("evidence", {}))
    for cid, ok, sev, det in new:
        ev[cid] = det
        (passed if ok else failed).append(cid)
        if not ok and sev == "critical":
            crit.append(cid)
    state = frozen.get("state")
    new_fail = [c for c, ok, _, _ in new if not ok]
    if new_fail and state in ("PASS", None):
        state = "HALT" if any(c in crit for c in new_fail) else "ESCALATE"
    elif new_fail and state == "ESCALATE" and any(c in crit for c in new_fail):
        state = "HALT"
    return {**frozen, "state": state, "passed": passed, "failed": failed, "warned": warned, "critical_failures": crit, "evidence": ev,
            "summary": f"{len(passed)} checks passed, {len(failed)} failed ({len(crit)} critical), {len(warned)} warned"}
