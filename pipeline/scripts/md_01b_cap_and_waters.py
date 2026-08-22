"""
md_01b_cap_and_waters.py — cap the artificial chain break and reinstate the
crystallographic waters, producing the production-ready receptor.

Two corrections to the benchmark build:

1. **Capping.** Removing BRIL leaves an artificial break in ICL3. Left as-is,
   tleap gives each side a charged terminus (NH3+/COO-), placing two formal
   charges in the middle of an intracellular loop that does not have them in
   the real receptor. ACE (acetyl) and NME (N-methylamide) caps neutralise the
   break. The true N- and C-termini of the construct are capped for the same
   reason: they are truncation points, not real chain ends.

   Only the linking atom of each cap is placed here, extrapolated from local
   backbone geometry; tleap rebuilds the remaining cap atoms from its own
   internal coordinates. No side-chain or backbone coordinate of the receptor
   is invented.

2. **Crystallographic waters.** Liu et al. (Science 2012) report 57 ordered
   waters in three clusters inside 4EIY, including the cluster harbouring the
   allosteric Na+. Discarding them and re-solvating from scratch would throw
   away experimentally resolved structure. Waters within a cutoff of the
   receptor are reinstated.

Usage:
    python md_01b_cap_and_waters.py --out <dir>
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def _root() -> Path:
    marker = "00_project_control/08_DECISION_LOG.md"
    for cand in (Path("/mnt/d/Agentic VS Workflow"), Path(r"D:\Agentic VS Workflow")):
        if (cand / marker).is_file():
            return cand
    raise SystemExit("project root not found")


ROOT = _root()
BENCH = ROOT / "20_new_runs/md_benchmark_2026-08-21/01_receptor"
WATER_CUTOFF = 5.0      # A from any receptor heavy atom
BOND_MAX = 2.0          # A; a real peptide C-N bond is ~1.33 A


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def parse(line: str) -> dict:
    return {
        "record": line[:6].strip(), "serial": line[6:11], "name": line[12:16].strip(),
        "resname": line[17:20].strip(), "chain": line[21], "resseq": int(line[22:26]),
        "xyz": np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
        "raw": line,
    }


def fmt(record: str, serial: int, name: str, resname: str, chain: str,
        resseq: int, xyz: np.ndarray, element: str) -> str:
    nm = f" {name:<3}" if len(name) < 4 else name
    return (f"{record:<6}{serial:>5} {nm}{'':1}{resname:>3} {chain}{resseq:>4}    "
            f"{xyz[0]:>8.3f}{xyz[1]:>8.3f}{xyz[2]:>8.3f}{1.00:>6.2f}{0.00:>6.2f}"
            f"{'':>10}{element:>2}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--receptor", default=str(BENCH / "4eiy_prep_keepion.pdb"))
    ap.add_argument("--waters", default=str(BENCH / "crystal_waters.pdb"))
    ap.add_argument("--out", default=str(ROOT / "20_new_runs/md_production_2026-08-21/01_receptor"))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    lines = Path(args.receptor).read_text(encoding="utf-8").splitlines()
    atoms = [parse(l) for l in lines if l.startswith("ATOM")]
    hets = [l for l in lines if l.startswith("HETATM")]

    def get(resseq: int, name: str) -> np.ndarray | None:
        for a in atoms:
            if a["resseq"] == resseq and a["name"] == name:
                return a["xyz"]
        return None

    # Segment on PEPTIDE-BOND DISTANCE, not on residue numbering. pdb4amber
    # renumbers consecutively, so the BRIL excision leaves no gap in the
    # numbers -- only a ~10.8 A gap in space between the C of one residue and
    # the N of the next. A numbering-based split silently finds nothing.
    resids = sorted({a["resseq"] for a in atoms})
    segments: list[tuple[int, int]] = []
    breaks: list[tuple[int, int, float]] = []
    s = resids[0]
    for i in range(len(resids) - 1):
        c = get(resids[i], "C")
        n = get(resids[i + 1], "N")
        d = float(np.linalg.norm(n - c)) if (c is not None and n is not None) else 999.0
        if d > BOND_MAX:
            segments.append((s, resids[i]))
            breaks.append((resids[i], resids[i + 1], round(d, 3)))
            s = resids[i + 1]
    segments.append((s, resids[-1]))

    def perp(v: np.ndarray) -> np.ndarray:
        """Any unit vector perpendicular to v."""
        a = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(a, v)) > 0.9:
            a = np.array([0.0, 1.0, 0.0])
        w = np.cross(v, a)
        return w / np.linalg.norm(w)

    # Build COMPLETE caps (all heavy atoms). A cap carrying only its linking
    # atom leaves a fractional formal charge, which pdb2pqr rejects outright
    # ("non-integer charge ... exceeding error tolerance") and which would
    # otherwise be discovered only after the whole build had run.
    caps: list[tuple[int, str, list[tuple[str, np.ndarray, str]]]] = []
    for (first, last) in segments:
        n, ca = get(first, "N"), get(first, "CA")
        if n is not None and ca is not None:
            v = (n - ca) / np.linalg.norm(n - ca)
            w = perp(v)
            c_at = n + 1.33 * v                 # carbonyl C, bonded to backbone N
            dn = -v                             # from C back toward N
            o_at = c_at + 1.23 * (-0.5 * dn + 0.866 * w)     # trigonal planar
            me_at = c_at + 1.50 * (-0.5 * dn - 0.866 * w)
            caps.append((first, "ACE", [("C", c_at, "C"), ("O", o_at, "O"),
                                        ("CH3", me_at, "C")]))
        c, ca_l = get(last, "C"), get(last, "CA")
        if c is not None and ca_l is not None:
            v = (c - ca_l) / np.linalg.norm(c - ca_l)
            w = perp(v)
            n_at = c + 1.33 * v                 # amide N, bonded to backbone C
            me_at = n_at + 1.47 * (0.5 * v + 0.866 * w)
            caps.append((last, "NME", [("N", n_at, "N"), ("CH3", me_at, "C")]))

    # Waters within the cutoff of any receptor heavy atom.
    prot_xyz = np.array([a["xyz"] for a in atoms])
    kept_wat: list[str] = []
    wat_lines = [l for l in Path(args.waters).read_text(encoding="utf-8").splitlines()
                 if l.startswith(("ATOM", "HETATM"))]
    for wl in wat_lines:
        w = parse(wl)
        if w["name"] not in ("O", "OW"):
            continue
        if np.min(np.linalg.norm(prot_xyz - w["xyz"], axis=1)) <= WATER_CUTOFF:
            kept_wat.append(wl)

    # Emit: segments with caps interleaved, then hetero, then waters.
    serial = 1
    outl: list[str] = []
    ace = {c[0]: c for c in caps if c[1] == "ACE"}
    nme = {c[0]: c for c in caps if c[1] == "NME"}
    for (first, last) in segments:
        if first in ace:
            for nm, xyz, el in ace[first][2]:
                outl.append(fmt("ATOM", serial, nm, "ACE", "A", first - 1, xyz, el))
                serial += 1
        for a in atoms:
            if first <= a["resseq"] <= last:
                outl.append(a["raw"][:6] + f"{serial:>5}" + a["raw"][11:])
                serial += 1
        if last in nme:
            for nm, xyz, el in nme[last][2]:
                outl.append(fmt("ATOM", serial, nm, "NME", "A", last + 1, xyz, el))
                serial += 1
        outl.append("TER")
    for h in hets:
        outl.append(h[:6] + f"{serial:>5}" + h[11:]); serial += 1
    outl.append("TER")
    for wl in kept_wat:
        outl.append(wl[:6] + f"{serial:>5}" + wl[11:]); serial += 1
    outl.append("END")

    p_out = out / "4eiy_capped_waters.pdb"
    p_out.write_text("\n".join(outl) + "\n", encoding="utf-8")

    meta = {
        "input_receptor": args.receptor,
        "protein_segments": segments,
        "chain_breaks_detected": [{"after": a, "before": b, "C_N_distance_A": d} for a,b,d in breaks],
        "caps_added": [{"kind": k, "anchor_resseq": r, "atoms": [a[0] for a in ats]} for r, k, ats in caps],
        "n_caps": len(caps),
        "cap_rationale": "BRIL removal leaves an artificial ICL3 break; charged "
                         "termini there would place formal charges the real "
                         "receptor does not have. Construct N/C termini capped "
                         "for the same reason. Only the linking atom is placed; "
                         "tleap rebuilds the rest.",
        "water_cutoff_A": WATER_CUTOFF,
        "crystal_waters_available": len(wat_lines),
        "crystal_waters_retained": len(kept_wat),
        "water_rationale": "Liu et al. Science 2012 (10.1126/science.1219218) report "
                           "57 ordered waters in three clusters in 4EIY, including the "
                           "cluster around the allosteric Na+.",
        "membrane_decision": "Pure POPC, no cholesterol. Ng, Laughton & Doughty "
                             "(J Chem Inf Model 2014;54:573-81, 10.1021/ci400463z) "
                             "simulated A2AR in POPC (300 ns, triplicate), giving a "
                             "direct precedent for this receptor. Cholesterol was "
                             "considered and rejected for this experiment: Yang & Lyman "
                             "(Biochemistry 2019;58:4096-4105, 10.1021/acs.biochem.9b00607) "
                             "show the local lipid environment around A2AR, cholesterol "
                             "especially, requires microseconds to equilibrate, which a "
                             "50-100 ns pose-persistence run cannot deliver.",
        "output": {"path": str(p_out), "sha256": sha256(p_out)},
    }
    (out / "receptor_prep_v2_metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

    print(f"segments        : {segments}")
    print(f"caps added      : {len(caps)}  ({[k for _,k,_ in caps]})")
    print(f"crystal waters  : {len(kept_wat)} retained of {len(wat_lines)} (<= {WATER_CUTOFF} A)")
    print(f"hetero kept     : {len(hets)}")
    print(f"\nwritten -> {p_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
