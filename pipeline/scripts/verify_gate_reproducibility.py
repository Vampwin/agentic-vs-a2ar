"""verify_gate_reproducibility.py -- the check behind the byte-identical claim.

Section 3.1 and Table 1 state that running the gate pipeline twice produced
contract outcomes that were byte-identical after stripping wall-clock
timestamps. The two runs are on disk as `snp_screen_repro_run1_2026-08-15` and
`snp_screen_repro_run2_2026-08-15`. Until 2026-08-29 no script in this project
compared them: the comparison had been performed once, by hand, and its result
was not preserved. A paper whose argument is that checks must leave artifacts
should not rest a reproducibility claim on a check that left none. This is that
script (ADR-AIIM-142).

It answers three questions, in order, because the honest answer to the first is
narrower than the claim as a reader would hear it.

**1. Are the contract outcomes identical?** The registry carries four
`contract_outcome` records per run -- the gate decisions themselves. They are
compared after removing the fields that carry a wall-clock reading: `timestamp`,
and the `record_id`, which embeds one. Nothing else is normalised away.

**2. Is every artifact identical?** No, and this is the part the prose did not
say. Two runs of the same code on the same input produce receptors that differ.
The script classifies every differing line rather than reporting a hash
mismatch, because the classification is the finding: polar-hydrogen placement is
non-deterministic, and in three atoms of one histidine the AutoDock atom *type*
differs between runs at identical coordinates.

**3. Does the difference reach the ligand?** Asked empirically, not argued. The
script measures the distance from every redocked pose to the atoms that differ
in type. If they lie beyond the interaction cutoff for every pose, identical
docking energies follow without any appeal to how the scoring function treats
those types.

Usage:
    python verify_gate_reproducibility.py [--out DIR]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path
from kdss_paths import project_root, resolve, under, dev_dir

# --- pre-registered criteria, transcribed from ADR-AIIM-142 before the run ----
# The manuscript's claim is about contract outcomes. That is the pass/fail
# criterion, and it is the only one. The artifact comparison below is reported
# whatever it says; it is a description, not a test with a threshold.
CONTRACT_RECORD_TYPE = "contract_outcome"
TIMESTAMP_FIELDS = ("timestamp", "record_id")   # record_id embeds a clock reading
VINA_INTERACTION_CUTOFF_A = 8.0                 # Vina's receptor-atom cutoff

# Stripping has to recurse. The first version of this script popped these two
# fields at the top level only and reported the contract outcomes as DIFFERING;
# the whole difference was a `timestamp` nested one level down, inside
# `state_history[0]`. The claim was right and the normaliser was too shallow.
#
# Recursion is the narrow reading of "stripping wall-clock timestamps", not a
# broad one: only keys with these exact names are removed, at any depth.
# Outcomes, state transitions, check counts, evidence lists and artifact hashes
# are all left in place and all take part in the comparison. If the gates ever
# disagree about a decision, this comparison sees it.

RUN1 = "snp_screen_repro_run1_2026-08-15"
RUN2 = "snp_screen_repro_run2_2026-08-15"


def _root() -> Path:
    # Delegates to kdss_paths, which finds whichever layout this file is inside
    # -- the development tree or a deposit -- by walking up from its own
    # location. It previously tried two absolute paths and exited if neither
    # existed, which on the author's machine meant a copy of the deposit read
    # and wrote the development tree instead of itself.
    return project_root()


ROOT = _root()
RUNS = dev_dir("20_new_runs")


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def strip_run_name(s: str) -> str:
    """Absolute paths differ only because the runs live in different folders."""
    return re.sub(r"snp_screen_repro_run[12]_2026-08-15", "<RUN>", s)


# --- 1. contract outcomes ----------------------------------------------------

def strip_timestamps(obj):
    """Remove TIMESTAMP_FIELDS at every depth. Nothing else is touched."""
    if isinstance(obj, dict):
        return {k: strip_timestamps(v) for k, v in obj.items()
                if k not in TIMESTAMP_FIELDS}
    if isinstance(obj, list):
        return [strip_timestamps(v) for v in obj]
    return obj


def contract_outcomes(run: Path) -> list[dict]:
    out = []
    for line in (run / "artifact_registry.jsonl").read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        rec = json.loads(line)
        if rec.get("record_type") != CONTRACT_RECORD_TYPE:
            continue
        rec = strip_timestamps(rec)
        out.append(json.loads(strip_run_name(json.dumps(rec, sort_keys=True))))
    return out


# --- 2. artifact comparison --------------------------------------------------

def classify_pdbqt(a: Path, b: Path) -> dict:
    """Where two prepared receptors differ, and in what."""
    la = a.read_text(encoding="utf-8").splitlines()
    lb = b.read_text(encoding="utf-8").splitlines()
    if len(la) != len(lb):
        return {"comparable": False, "reason": "different line counts"}

    remark_only, coord, typing, other = 0, [], [], []
    for x, y in zip(la, lb):
        if x == y:
            continue
        if x.startswith("REMARK") and strip_run_name(x) == strip_run_name(y):
            remark_only += 1
            continue
        if not (x.startswith(("ATOM", "HETATM")) and y.startswith(("ATOM", "HETATM"))):
            other.append((x, y))
            continue
        name = x[12:16].strip()
        res = f"{x[17:20].strip()} {x[21]}{x[22:26].strip()}"
        xyz_x, xyz_y = x[30:54], y[30:54]
        type_x, type_y = x[77:79].strip(), y[77:79].strip()
        if xyz_x != xyz_y and type_x == type_y:
            coord.append({"atom": name, "residue": res, "autodock_type": type_x})
        elif xyz_x == xyz_y and type_x != type_y:
            typing.append({"atom": name, "residue": res,
                           "type_run1": type_x, "type_run2": type_y,
                           "xyz": " ".join(xyz_x.split())})
        else:
            other.append((x, y))

    by_type: dict[str, int] = {}
    for c in coord:
        by_type[c["autodock_type"]] = by_type.get(c["autodock_type"], 0) + 1

    return {
        "comparable": True,
        "path_only_remark_lines": remark_only,
        "coordinate_differences": len(coord),
        "coordinate_differences_by_autodock_type": by_type,
        "atom_type_differences": typing,
        "unclassified_differences": len(other),
    }


def compare_files(r1: Path, r2: Path) -> tuple[list[dict], list[str]]:
    rows, differing = [], []
    files = sorted(p.relative_to(r1).as_posix() for p in r1.rglob("*") if p.is_file())
    assert files == sorted(p.relative_to(r2).as_posix()
                           for p in r2.rglob("*") if p.is_file()), \
        "the two runs do not contain the same files"
    for rel in files:
        a, b = r1 / rel, r2 / rel
        ha, hb = sha256(a), sha256(b)
        row = {"file": rel, "sha256_run1": ha, "sha256_run2": hb,
               "identical_raw": ha == hb}
        if ha != hb:
            differing.append(rel)
            if rel.endswith(".pdbqt"):
                row["difference"] = classify_pdbqt(a, b)
            elif rel.endswith(".json"):
                da = json.loads(strip_run_name(a.read_text(encoding="utf-8")))
                db = json.loads(strip_run_name(b.read_text(encoding="utf-8")))
                for d in (da, db):
                    for k in list(d):
                        if "timestamp" in k.lower():
                            d.pop(k)
                row["difference"] = {
                    "identical_after_stripping_timestamps_and_paths": da == db,
                    "differing_keys": sorted(k for k in set(da) | set(db)
                                             if da.get(k) != db.get(k)),
                }
        rows.append(row)
    return rows, differing


# --- 3. does the typing difference reach the ligand? -------------------------

def distance_to_poses(run: Path, atoms: list[dict]) -> dict | None:
    """Closest approach of any redocked pose atom to each mis-typed atom."""
    if not atoms:
        return None
    try:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
    except ImportError:
        return {"measured": False, "reason": "rdkit unavailable"}

    sdf = run / "G3_redocking/ZMA_redocked_4EIY.sdf"
    pose_xyz = []
    for mol in Chem.SDMolSupplier(str(sdf), removeHs=True):
        if mol is None:
            continue
        pose_xyz.extend(list(mol.GetConformer().GetPositions()))
    if not pose_xyz:
        return {"measured": False, "reason": f"no poses read from {sdf.name}"}

    import numpy as np
    P = np.asarray(pose_xyz)
    out = []
    for at in atoms:
        q = np.asarray([float(v) for v in at["xyz"].split()])
        out.append({"atom": at["atom"], "residue": at["residue"],
                    "min_distance_to_any_pose_atom": round(
                        float(np.linalg.norm(P - q, axis=1).min()), 2)})
    return {
        "measured": True,
        "n_pose_atoms": int(P.shape[0]),
        "cutoff_used": VINA_INTERACTION_CUTOFF_A,
        "atoms": out,
        "all_beyond_cutoff": all(
            a["min_distance_to_any_pose_atom"] > VINA_INTERACTION_CUTOFF_A
            for a in out),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(under("20_new_runs", "gate_reproducibility_2026-08-29")))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    r1, r2 = RUNS / RUN1, RUNS / RUN2
    for r in (r1, r2):
        if not r.is_dir():
            raise SystemExit(f"missing run directory: {r}")

    co1, co2 = contract_outcomes(r1), contract_outcomes(r2)
    contracts_match = co1 == co2

    rows, differing = compare_files(r1, r2)

    mistyped: list[dict] = []
    for row in rows:
        d = row.get("difference", {})
        if isinstance(d, dict):
            mistyped.extend(d.get("atom_type_differences", []))
    reach = distance_to_poses(r1, mistyped)

    report = {
        "question": "Two executions of the same gate pipeline on the same "
                    "frozen inputs. What is reproducible, and what is not?",
        "runs": {"run1": RUN1, "run2": RUN2},
        "criterion": "The manuscript's claim is about contract outcomes. Those "
                     "are compared after removing `timestamp` and `record_id`, "
                     "the two fields carrying a wall-clock reading, and nothing "
                     "else. Everything below the contract check is reported as "
                     "found, without a threshold.",
        "contract_outcomes": {
            "n_records_run1": len(co1),
            "n_records_run2": len(co2),
            "identical_after_stripping_timestamps": contracts_match,
            "records": co1 if contracts_match else {"run1": co1, "run2": co2},
        },
        "artifacts": {
            "n_files": len(rows),
            "n_identical_raw": sum(1 for r in rows if r["identical_raw"]),
            "n_differing": len(differing),
            "differing": differing,
            "detail": [r for r in rows if not r["identical_raw"]],
        },
        "reach_of_the_typing_difference": reach,
        "reading": "",
    }

    n_coord = sum(r.get("difference", {}).get("coordinate_differences", 0)
                  for r in rows if isinstance(r.get("difference"), dict))
    report["reading"] = (
        f"The four contract outcomes are "
        f"{'identical' if contracts_match else 'NOT identical'} across the two "
        f"runs once timestamps are removed, which is what Section 3.1 claims. "
        f"The artifacts are not: {len(differing)} of {len(rows)} files differ, "
        f"including the prepared receptor, where {n_coord} atom coordinates and "
        f"{len(mistyped)} AutoDock atom types differ between two runs of the "
        f"same code on the same input. The gate decisions are therefore "
        f"invariant under an upstream non-determinism this project had not "
        f"documented. That is a stronger result than bit-identical files would "
        f"have been, and a different one from what the prose implied.")

    dest = out / "gate_reproducibility.json"
    dest.write_text(json.dumps(report, indent=2), encoding="utf-8")

    print(f"contract outcomes ({len(co1)} records) : "
          f"{'IDENTICAL' if contracts_match else 'DIFFER'}")
    print(f"artifacts                          : "
          f"{len(rows) - len(differing)}/{len(rows)} byte-identical")
    for rel in differing:
        print(f"    differs: {rel}")
    for r in rows:
        d = r.get("difference")
        if isinstance(d, dict) and d.get("comparable"):
            print(f"\n  {r['file']}")
            print(f"    path-only REMARK lines : {d['path_only_remark_lines']}")
            print(f"    coordinate differences : {d['coordinate_differences']} "
                  f"{d['coordinate_differences_by_autodock_type']}")
            print(f"    atom-type differences  : {len(d['atom_type_differences'])}")
            for t in d["atom_type_differences"]:
                print(f"        {t['residue']:>12} {t['atom']:<5} "
                      f"{t['type_run1']} -> {t['type_run2']}  (same coordinates)")
            if d["unclassified_differences"]:
                print(f"    UNCLASSIFIED           : {d['unclassified_differences']}")
    if reach and reach.get("measured"):
        print(f"\n  closest approach of any of {reach['n_pose_atoms']} redocked "
              f"pose atoms to the mis-typed atoms:")
        for a in reach["atoms"]:
            print(f"        {a['residue']:>12} {a['atom']:<5} "
                  f"{a['min_distance_to_any_pose_atom']:>6.2f} A")
        print(f"    all beyond the {reach['cutoff_used']} A cutoff: "
              f"{reach['all_beyond_cutoff']}")
    print(f"\nwritten -> {dest}")
    return 0 if contracts_match else 1


if __name__ == "__main__":
    raise SystemExit(main())
