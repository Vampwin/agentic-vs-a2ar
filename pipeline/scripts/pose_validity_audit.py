"""pose_validity_audit.py -- are the docked poses physically possible at all?

Every scoring layer in this pipeline answers "which compound is better".
**Nothing in it answers "is this pose chemically valid".** PoseBusters does, with
22 checks that carry no ranking information whatsoever: bond lengths and angles
against reference geometry, aromatic-ring flatness, internal energy, atoms
clashing with themselves, and volume overlap with the receptor.

That distinction is the manuscript's own diagnosis applied to itself. A pose that
fails these checks cannot support any conclusion drawn from its geometry, and no
score in this project would have revealed it.

**This is descriptive and carries no pass/fail rule for compound selection.** The
tool's own defaults are used unmodified; nothing here is tuned. It is run over
every full-protocol pose in the screen, both engines, and the ZM241385 control,
so the answer is a population statement rather than the n = 3 the dynamics give.

Reference: Buttenschoen M, Morris GM, Deane CM. PoseBusters. Chem Sci. 2024.
Installed version is recorded in the output.

Usage:
    python pose_validity_audit.py [--out DIR] [--limit N]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import pandas as pd  # noqa: E402
from posebusters import PoseBusters  # noqa: E402


def _root() -> Path:
    marker = "00_project_control/08_DECISION_LOG.md"
    for cand in (Path(r"D:\Agentic VS Workflow"), Path("/mnt/d/Agentic VS Workflow")):
        if (cand / marker).is_file():
            return cand
    raise SystemExit("project root not found")


ROOT = _root()
RUNS = ROOT / "20_new_runs"
FUNNEL = RUNS / "np_library_funnel_2026-08-15"
DIFFDOCK = RUNS / "diffdock_orthogonal_2026-08-23"
RECEPTOR = DIFFDOCK / "4EIY_protein_only.pdb"
CONTROL_POSE = RUNS / "snp_screen_2026-08-04/G3_redocking/ZMA_redocked_4EIY.sdf"
CANDIDATES = ["NPLIB01185", "NPLIB01643", "NPLIB03056"]

# Checks that describe the file rather than the geometry. Kept in the output but
# excluded from the "physically valid" verdict so a loading quirk cannot be
# reported as a chemistry failure.
LOADING = {"mol_pred_loaded", "mol_cond_loaded", "sanitization", "inchi_convertible"}


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def compound_name(p: Path) -> str:
    """DiffDock writes every pose as `rank1.sdf`, so the stem is the rank and the
    compound is the directory. Taking the stem regardless labelled all four
    DiffDock poses "rank1" and collapsed them into one compound in the summary."""
    return p.parent.name if p.stem.lower().startswith("rank") else p.stem


def bust(pb, paths, tag):
    """Run PoseBusters over a list of SDFs, returning a tidy frame."""
    frames = []
    for i, p in enumerate(paths, 1):
        name = compound_name(p)
        try:
            df = pb.bust([p], None, RECEPTOR)
        except Exception as exc:            # a broken file is data, not a crash
            frames.append(pd.DataFrame([{"compound": name, "engine": tag,
                                         "pose_rank": None, "error": str(exc)[:200]}]))
            continue
        df = df.reset_index(drop=True)
        df.insert(0, "pose_rank", range(1, len(df) + 1))
        df.insert(0, "engine", tag)
        df.insert(0, "compound", name)
        frames.append(df)
        if i % 50 == 0:
            print(f"  {tag}: {i}/{len(paths)} files", flush=True)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def summarise(df: pd.DataFrame, checks: list[str]) -> dict:
    if df.empty:
        return {"n_poses": 0}
    ok = df[checks].apply(lambda r: bool(r.all()), axis=1)
    per_check = {c: int((~df[c].astype("boolean").fillna(False)).sum()) for c in checks}
    return {
        "n_poses": int(len(df)),
        "n_compounds": int(df["compound"].nunique()),
        "n_poses_all_checks_pass": int(ok.sum()),
        "fraction_valid": round(float(ok.mean()), 4),
        "failures_by_check": {k: v for k, v in sorted(
            per_check.items(), key=lambda kv: -kv[1]) if v},
        "compounds_with_any_failing_pose": sorted(
            df.loc[~ok, "compound"].unique().tolist())[:50],
        "n_compounds_with_any_failing_pose": int(df.loc[~ok, "compound"].nunique()),
        "rank1_fraction_valid": round(float(
            ok[df["pose_rank"] == 1].mean()), 4) if (df["pose_rank"] == 1).any() else None,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(RUNS / "pose_validity_2026-09-03"))
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    if not RECEPTOR.is_file():
        raise SystemExit(f"receptor missing: {RECEPTOR}")

    vina_poses = sorted(FUNNEL.glob("full_poses_chunk*/*.sdf"))
    if args.limit:
        vina_poses = vina_poses[:args.limit]
    dd_poses = sorted(DIFFDOCK.glob("results/*/rank1.sdf"))
    print(f"receptor {RECEPTOR.name}")
    print(f"vina full-protocol poses: {len(vina_poses)} files")
    print(f"diffdock poses:           {len(dd_poses)} files")

    import posebusters
    pb = PoseBusters("dock")

    frames = {"vina_screen": bust(pb, vina_poses, "vina"),
              "diffdock": bust(pb, dd_poses, "diffdock")}
    if CONTROL_POSE.is_file():
        frames["control_vina_redock"] = bust(pb, [CONTROL_POSE], "vina_control")

    allf = pd.concat([f for f in frames.values() if not f.empty], ignore_index=True)
    checks = [c for c in allf.columns
              if c not in ("compound", "engine", "pose_rank", "error")
              and c not in LOADING]

    allf.to_csv(out / "pose_validity_per_pose.csv", index=False)

    report = {
        "question": "Are the docked poses physically valid? A separate question "
                    "from whether they rank well, and one nothing else in this "
                    "project asks.",
        "tool": {"name": "posebusters", "version": posebusters.__version__,
                 "config": "dock", "defaults": "unmodified",
                 "reference": "Buttenschoen, Morris & Deane, Chem Sci 2024"},
        "descriptive_only": "No pass/fail rule is applied to compound selection. "
                            "No threshold was tuned. The tool's own defaults "
                            "decide every check.",
        "receptor": {"file": RECEPTOR.relative_to(ROOT).as_posix(),
                     "sha256": sha256(RECEPTOR)},
        "checks_applied": checks,
        "checks_excluded_from_verdict": sorted(LOADING),
        "by_set": {k: summarise(v, checks) for k, v in frames.items() if not v.empty},
        "candidates": {},
    }
    for cid in CANDIDATES:
        sub = allf[allf.compound == cid]
        if sub.empty:
            continue
        report["candidates"][cid] = {
            eng: summarise(sub[sub.engine == eng], checks)
            for eng in sub.engine.unique()
        }

    dest = out / "pose_validity.json"
    dest.write_text(json.dumps(report, indent=2), encoding="utf-8")

    print()
    for k, s in report["by_set"].items():
        print(f"{k:<22} {s['n_poses']:>6} poses  {s['n_compounds']:>4} compounds  "
              f"valid {s['fraction_valid']:.4f}"
              + (f"  rank1 {s['rank1_fraction_valid']:.4f}"
                 if s.get("rank1_fraction_valid") is not None else ""))
        if s["failures_by_check"]:
            for c, n in list(s["failures_by_check"].items())[:6]:
                print(f"     {c:<42} {n} pose(s) fail")
    print(f"\nwritten -> {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
