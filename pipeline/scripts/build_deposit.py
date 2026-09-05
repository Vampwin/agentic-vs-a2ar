"""
build_deposit.py -- assemble the public code-and-data deposit for the v6
manuscript as a git-ready tree.

Three rules govern what goes in.

**Nothing that is already archived elsewhere is re-deposited.** LOTUS v11 is
107 MB of gzipped CSV that already carries its own Zenodo DOI; ChEMBL and DUD-E
are public releases. Copying them here would add nothing but weight and a
second, divergeable copy. They are cited with their identifiers and the
SHA-256 this project verified them against, so a reader can confirm they have
the same bytes we did.

**Nothing sealed goes in.** The held-out fault families are excluded. They are
not referenced anywhere in v6, which reports development families only, and
CLAUDE.md requires a separate human approval to expose them. The v4 supplement
carried some of that material; this deposit does not inherit that decision.

**Nothing is regenerated.** Every artifact is copied as the analysis read it,
with its SHA-256 recorded, because a deposit that recomputes its own numbers
demonstrates nothing about the run that produced them.

Every file is scanned before copying: assay-derived values, credentials, and
absolute paths from the author's machine all abort the build rather than
producing a deposit that has to be retracted.

Usage:
    python build_deposit.py [--out DIR]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import sys
from pathlib import Path


def _root() -> Path:
    marker = "00_project_control/08_DECISION_LOG.md"
    for cand in (Path(r"D:\Agentic VS Workflow"), Path("/mnt/d/Agentic VS Workflow")):
        if (cand / marker).is_file():
            return cand
    raise SystemExit("project root not found")


ROOT = _root()
DEFAULT_OUT = ROOT / "40_deposit"

# (source, destination, description)
TREE: list[tuple[str, str, str]] = [
    ("10_new_system/contracts", "pipeline/contracts",
     "the ten frozen G0-G8 audit contracts"),
    ("10_new_system/evaluators", "pipeline/evaluators",
     "the deterministic evaluator modules the agent cannot override"),
    ("10_new_system/fault_injection", "pipeline/fault_injection",
     "the eighteen fault-family definitions (eleven implemented, seven declared stubs)"),
    ("10_new_system/tests", "pipeline/tests",
     "the test suite: 675 tests across 30 modules"),
    ("10_new_system/scripts", "pipeline/scripts",
     "every analysis, figure and audit script referenced in the manuscript"),
    ("10_new_system/kdss_audit", "pipeline/kdss_audit",
     "the artifact registry implementation"),
    # These three were missed in the first build: deposited scripts and tests
    # import them, so without them a clean checkout cannot run its own test
    # suite. `benchmark/heldout_interface` stays excluded by SKIP_DIRS.
    ("10_new_system/risk_monitor", "pipeline/risk_monitor",
     "the learned risk monitor, the component the A/C ablation conditions test"),
    ("10_new_system/telemetry", "pipeline/telemetry",
     "the out-of-domain detector and its telemetry features"),
    ("10_new_system/benchmark", "pipeline/benchmark",
     "SEER scoring and the development-fault benchmark harness"),
    ("30_manuscript/supplement/S10_verification_2026-08-21", "verification",
     "the 2026-08-21 verification evidence behind Sections 2 and 3"),
    ("30_manuscript/figures", "figures",
     "figure sources and their provenance records"),
    ("20_new_runs/verification_summary_2026-08-21", "verification/summary",
     "per-gate mutation scores, structure and protocol checks, the G8 audit"),
    # NOTE: verification/summary and verification/structure_and_protocol share
    # two files; the duplicate is resolved after copying, below.
    ("20_new_runs/metamorphic_2026-08-21", "verification/metamorphic",
     "the four metamorphic invariants, their inputs and their outputs"),
    # DELIBERATELY NOT DEPOSITED: 20_new_runs/heldout_evaluation_2026-08-08.
    # The manuscript reports the held-out RESULT (one of five, with each miss's
    # root cause, §3.1). The faulted workdirs are a different thing: the
    # corrupted artifacts themselves let H01-H05 be reconstructed, which is
    # unsealing, and CLAUDE.md requires separate explicit human approval for
    # that. The Data Availability statement already offers them on request.
    # Depositing them was attempted on 2026-08-23 and reverted (ADR-AIIM-113).
    ("20_new_runs/diffdock_orthogonal_2026-08-23", "verification/orthogonal_engine",
     "the pre-registered DiffDock cross-check: inputs, all predicted poses, "
     "and the machine-readable result the manuscript's numbers derive from"),
]

# Directories never copied, wherever they appear.
SKIP_DIRS = {"__pycache__", ".pytest_cache", ".ipynb_checkpoints", "held_out",
             "heldout_interface", "poses", "fastpass_poses", "_superseded",
             "_tools"}
SKIP_SUFFIX = {".pyc", ".pyo", ".log", ".tpr", ".trr", ".xtc", ".cpt", ".edr"}

# Files never copied, by name. The two held-out scripts matter most: excluding
# the sealed DIRECTORY was not enough, because these live in scripts/ alongside
# everything else and would have shipped. Sealing a folder does not seal a
# topic.
SKIP_NAMES = {
    "run_heldout_evaluation_2026-08-08.py",
    "verify_h01_h02_h03_fixes_2026-08-08.py",
}

# Working files created while debugging, prefixed with an underscore by
# convention in this project. They are not deliverables.
#
# Evaluator helper modules are the exception, and getting this wrong was not
# cosmetic. They follow the same underscore convention but ARE deliverables:
# excluding them removed _domain_gate_bias_probe.py -- which the manuscript
# offers as a contract "directly transferable to any screening campaign" -- and
# _assay_blindness.py, which implements a stated safety boundary, from every
# deposit ever built. A reader following the Data availability statement could
# not obtain either one. Found by the internal review of 2026-08-25
# (ADR-AIIM-127).
_KEEP_UNDERSCORED_IN = {"evaluators"}


def _is_scratch(f: Path) -> bool:
    if f.name == "__init__.py":
        return False
    if f.parent.name in _KEEP_UNDERSCORED_IN and f.suffix == ".py":
        return False
    return f.name.startswith("_")

# Screening lives in its own module: getting these patterns right matters
# more than anything else here, and it carries a self-test.
sys.path.insert(0, str(Path(__file__).parent))
from deposit_screen import (  # noqa: E402
    REVIEWED_AND_CLEARED, screen, self_test as screen_self_test,
)

EXTERNAL = [
    {"name": "LOTUS natural-product knowledge base",
     "version": "v11 frozen snapshot, 13 April 2026",
     "identifier": "doi:10.5281/zenodo.19360665", "licence": "CC BY 4.0",
     "why_not_deposited": "107 MB, already archived with its own DOI",
     "files": ["lotus_260413_frozen.csv.gz", "lotus_260413_frozen_metadata.csv.gz"]},
    {"name": "ChEMBL A2AR bioactivity records",
     "version": "see pipeline/contracts/G0_input.yaml for the frozen extract",
     "identifier": "doi:10.1093/nar/gky1075", "licence": "CC BY-SA 3.0",
     "why_not_deposited": "public database release; the frozen extract's SHA-256 "
                          "is recorded in the run metadata",
     "files": ["chembl_a2a_training_v3.csv", "chembl_a2ar_validation_set.csv"]},
    {"name": "DUD-E A2AR decoy set",
     "version": "as distributed", "identifier": "doi:10.1021/jm300687e",
     "licence": "see source", "why_not_deposited": "public benchmark release",
     "files": ["dude_a2ar_docking_set.csv"]},
    {"name": "A2AR crystal structures",
     "version": "as deposited", "identifier": "PDB 4EIY, 5OLH",
     "licence": "CC0", "why_not_deposited": "public structural archive",
     "files": ["4EIY.pdb", "5OLH.pdb"]},
]


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    args = ap.parse_args()
    if not screen_self_test():
        raise SystemExit('deposit screen failed its own self-test; not building a deposit on a screen that does not work')
    out = Path(args.out)
    # Clear the tree but PRESERVE version-control metadata. The first version
    # of this called shutil.rmtree(out), which deletes .git along with
    # everything else: rebuilding the deposit would have destroyed the commit
    # history that makes it citable. It never fired only because Windows marks
    # git object files read-only and the delete raised instead of succeeding.
    # On a POSIX host it would have succeeded silently.
    # .git carries the citable history; the other four are hand-authored and
    # are NOT regenerated by this script, so a rebuild must not delete them.
    KEEP = {".git", ".gitignore", ".gitattributes", "README.md",
            "CITATION_AND_LICENCE.md"}
    if out.exists():
        for child in out.iterdir():
            if child.name in KEEP:
                continue
            shutil.rmtree(child) if child.is_dir() else child.unlink()
    out.mkdir(parents=True, exist_ok=True)

    manifest: list[dict] = []
    refused: list[tuple[str, str]] = []
    cleared: list[tuple[str, str, str]] = []
    skipped = 0
    total = 0

    for src_rel, dst_rel, desc in TREE:
        src = ROOT / src_rel
        if not src.exists():
            raise SystemExit(f"missing source: {src_rel}")
        for f in sorted(src.rglob("*")):
            if not f.is_file():
                continue
            if any(part in SKIP_DIRS for part in f.parts):
                skipped += 1
                continue
            if f.suffix.lower() in SKIP_SUFFIX:
                skipped += 1
                continue
            if f.name in SKIP_NAMES or _is_scratch(f):
                skipped += 1
                continue
            rel_src = str(f.relative_to(ROOT).as_posix())
            reasons = screen(f)
            if reasons and rel_src in REVIEWED_AND_CLEARED:
                cleared.append((rel_src, "; ".join(reasons),
                                REVIEWED_AND_CLEARED[rel_src]))
                reasons = []
            if reasons:
                refused.append((rel_src, "; ".join(reasons)))
                continue
            rel = f.relative_to(src)
            dst = out / dst_rel / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(f, dst)
            size = dst.stat().st_size
            total += size
            manifest.append({
                "path": str((Path(dst_rel) / rel).as_posix()),
                "source": str(f.relative_to(ROOT).as_posix()),
                "bytes": size, "sha256": sha256(dst)})

    if cleared:
        print(f"flagged and cleared on review ({len(cleared)}):")
        for path, why, because in cleared:
            print(f"   {path}")
            print(f"      flagged: {why}")
            print(f"      cleared: {because[:88]}...")
        print()
    # One artifact, one path. Several files arrive by two routes -- the S10
    # supplement bundle and the run directory both carry them -- and two paths
    # for one artifact is exactly the divergence this deposit exists to rule
    # out. Deduplicate by content hash rather than by name, so a pair that
    # differs is kept and reported instead of silently collapsed.
    by_hash: dict[str, list[dict]] = {}
    for m in manifest:
        by_hash.setdefault(m["sha256"], []).append(m)
    # Prefer the topically-named directory over the catch-all summary one.
    def rank(entry: dict) -> tuple[int, str]:
        return (1 if "/summary/" in entry["path"] else 0, entry["path"])
    dupes = 0
    for entries in by_hash.values():
        if len(entries) < 2:
            continue
        entries.sort(key=rank)
        keep = entries[0]
        for drop in entries[1:]:
            # Two sources can map to the SAME destination path -- the S10
            # bundle and the run directory now both supply the DiffDock
            # artifacts under orthogonal_engine/. The second copy simply
            # overwrote the first, so there is one physical file and two
            # manifest entries. Unlinking "the duplicate" then deleted the
            # very file the surviving entry points at, and the deposit lost
            # orthogonal_check_result.json -- the artifact the manuscript
            # says its Section 3.5 numbers derive from -- while still listing
            # it. Only unlink when the paths actually differ.
            if drop["path"] != keep["path"]:
                (out / drop["path"]).unlink(missing_ok=True)
            manifest.remove(drop)
            dupes += 1

    # A manifest entry pointing at a file that is not there is the failure this
    # deposit exists to prevent, and the dedup step above silently produced
    # exactly that once. Fail loudly rather than write a manifest that lies.
    ghosts = [m["path"] for m in manifest if not (out / m["path"]).is_file()]
    if ghosts:
        raise SystemExit(
            "manifest lists files that are not on disk -- refusing to write it:\n  "
            + "\n  ".join(ghosts))

    if refused:
        print("REFUSED -- these files were not deposited:")
        for path, why in refused:
            print(f"   {path}\n      {why}")

    (out / "MANIFEST.json").write_text(json.dumps({
        "assembled": "2026-08-22",
        "n_files": len(manifest),
        "total_bytes": total,
        "n_files_excluded_by_rule": skipped,
        "n_duplicates_removed": dupes,
        "n_files_refused_by_screening": len(refused),
        "refused": [{"path": p, "reason": w} for p, w in refused],
        "flagged_and_cleared_on_review": [
            {"path": p, "flagged_because": w, "cleared_because": b}
            for p, w, b in cleared],
        "external_inputs_not_deposited": EXTERNAL,
        "files": manifest,
    }, indent=2), encoding="utf-8")

    print(f"\ndeposited {len(manifest)} files, {total / 1024 / 1024:.1f} MB -> {out}")
    print(f"excluded by rule: {skipped}   refused by screening: {len(refused)}")
    for grp in sorted({m['path'].split('/')[0] for m in manifest}):
        n = sum(1 for m in manifest if m["path"].startswith(grp + "/"))
        b = sum(m["bytes"] for m in manifest if m["path"].startswith(grp + "/"))
        print(f"   {grp:14s} {n:4d} files  {b/1024:8.0f} KB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
