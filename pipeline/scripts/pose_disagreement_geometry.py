"""pose_disagreement_geometry.py -- m10. What kind of disagreement is 8.80 A?

Section 3.6 reports that Vina and DiffDock place the three Tier-A candidates
6.43, 8.80 and 9.50 A apart, and Figure 5B is captioned "8.80 A apart". A
reviewer reads that as a displacement, looks at the panel, and sees two
molecules on top of each other. The panel was built to be checked that way, so
the check was run.

An in-place RMSD folds two different things together: how far the molecule
moved, and how differently it is oriented and folded once it gets there. The
centroid separation separates them. For NPLIB01643 the centroids are 3.78 A
apart -- less than the molecule's own radius of gyration -- against an RMSD of
8.80 A. The two poses occupy nearly the same space in nearly opposite
orientations. There is no gap for any viewpoint to reveal, which is why an
orthogonal render of that panel was built and then discarded: it showed the
same overlap from a second angle, correctly.

The three candidates are not alike in this respect, and the article had been
reporting them as though they were.

Usage:
    python pose_disagreement_geometry.py [--out DIR]
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from kdss_paths import project_root, resolve, under, dev_dir

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolAlign

RDLogger.DisableLog("rdApp.*")

CANDIDATES = ["NPLIB01185", "NPLIB01643", "NPLIB03056"]


def _root() -> Path:
    # Delegates to kdss_paths, which finds whichever layout this file is inside
    # -- the development tree or a deposit -- by walking up from its own
    # location. It previously tried two absolute paths and exited if neither
    # existed, which on the author's machine meant a copy of the deposit read
    # and wrote the development tree instead of itself.
    return project_root()


ROOT = _root()
RUNS = dev_dir("20_new_runs")
FUNNEL = under("20_new_runs", "np_library_funnel_2026-08-15")
DD = under("20_new_runs", "diffdock_orthogonal_2026-08-23/results")


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def vina_pose(cid: str) -> Path:
    for c in (1, 2, 3):
        p = FUNNEL / f"full_poses_chunk{c}/{cid}.sdf"
        if p.exists():
            return p
    raise SystemExit(f"no Vina pose archived for {cid}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(under("20_new_runs", "pose_disagreement_2026-08-29")))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    rows = {}
    for cid in CANDIDATES:
        vp, dp = vina_pose(cid), DD / cid / "rank1.sdf"
        v = Chem.MolFromMolFile(str(vp), removeHs=True)
        d = Chem.MolFromMolFile(str(dp), removeHs=True)
        if v is None or d is None:
            raise SystemExit(f"could not read a pose for {cid}")
        cv = np.array(v.GetConformer().GetPositions())
        cd = np.array(d.GetConformer().GetPositions())
        rmsd = float(rdMolAlign.CalcRMS(d, v))
        cdist = float(np.linalg.norm(cv.mean(0) - cd.mean(0)))
        rg = float(np.sqrt(((cv - cv.mean(0)) ** 2).sum(1).mean()))
        rows[cid] = {
            "vina_pose": str(vp.relative_to(ROOT).as_posix()),
            "vina_sha256": sha256(vp),
            "diffdock_pose": str(dp.relative_to(ROOT).as_posix()),
            "diffdock_sha256": sha256(dp),
            "n_heavy_atoms": int(v.GetNumHeavyAtoms()),
            "rmsd_in_place": round(rmsd, 2),
            "centroid_separation": round(cdist, 2),
            "radius_of_gyration_vina": round(rg, 2),
            "translated_fraction": round(cdist / rmsd, 2),
        }

    report = {
        "method": "In-place symmetry-corrected RMSD (rdMolAlign.CalcRMS), the "
                  "same convention used throughout, beside the separation of the "
                  "two heavy-atom centroids. Neither pose is superposed.",
        "reading": "The ratio is the part of the disagreement that is a "
                   "displacement. NPLIB01643, the compound Figure 5B shows and "
                   "the one DiffDock scores highest, is the least displaced of "
                   "the three: its centroids are closer together than the "
                   "molecule's own radius of gyration, so the two engines put it "
                   "in the same place facing differently. NPLIB03056 is the most "
                   "displaced. Reporting one RMSD each had made them look alike.",
        "candidates": rows,
    }
    (out / "pose_disagreement_geometry.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")

    print(f"{'compound':>12} {'RMSD':>7} {'centroid':>9} {'Rg':>6} {'displaced':>10}")
    for cid, r in rows.items():
        print(f"{cid:>12} {r['rmsd_in_place']:>6.2f}A {r['centroid_separation']:>8.2f}A "
              f"{r['radius_of_gyration_vina']:>5.2f}A {r['translated_fraction']:>10.2f}")
    print(f"\nwritten -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
