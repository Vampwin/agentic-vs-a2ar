"""md_18_recognition_motif.py -- does the predicted pose hold the interaction
that defines antagonist binding at A2A?

Criteria are module constants transcribed from ADR-AIIM-161 **before the
trajectories were read**, which is the convention R2, R3 and R4 followed and the
reason their outcomes could be written up without argument.

**Why not RMSD.** §3.6 rests on ligand RMSD, which is a whole-molecule average.
It cannot separate a ligand that rotated a flexible tail from one that released
the interaction that makes it an antagonist, and for NPLIB01185 it failed
outright on a periodic image (ADR-AIIM-156). This measures the interaction
directly.

**Why the standard is trustworthy.** The control is ZM241385 -- the co-crystal
antagonist -- simulated in the same receptor, same force field, same protocol,
only the ligand differing. The criterion comes from the native ligand measured
in this same system, not from the literature and not from our own docking.

**Two traps this script is built around.**

1. `gmx mindist` applies the minimum-image convention itself. Nothing here reads
   a trajectory representation, which is the failure that produced a fictitious
   dissociating ligand in ADR-AIIM-156.

2. **Residue numbers in this topology are not unique.** `resid 36` is LEU at
   atom 516 and ASN at atom 3806; 91 residue numbers are reused. Selecting an
   anchor by residue number measures the wrong residue and returns a
   plausible-looking number. Every group here is built from explicit atom
   indices, and the anchor is located by *sequence position* with its identity
   asserted afterwards.

Usage:
    python md_18_recognition_motif.py [--out DIR]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
from pathlib import Path
from kdss_paths import project_root, resolve, under, wsl_root, dev_dir

# --- pre-registered, ADR-AIIM-161 -------------------------------------------
HBOND_CUTOFF_A = 3.5          # heavy-atom donor-acceptor, conventional
CONTACT_CUTOFF_A = 4.5        # ring contact
SECONDARY_WINDOW_NS = 20.0    # discarded for the secondary window only
RHO_REPRODUCES = 0.5
RHO_DOES_NOT = 0.1
# Sidechain atoms that carry each interaction.
ANCHOR_ATOMS = {
    "N6.55":  ("ASN", ["OD1", "ND2"]),
    "E-ECL2": ("GLU", ["OE1", "OE2"]),
    "F-ECL2": ("PHE", ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]),
}
CUTOFF_FOR = {"N6.55": HBOND_CUTOFF_A, "E-ECL2": HBOND_CUTOFF_A,
              "F-ECL2": CONTACT_CUTOFF_A}
# The ligand side: polar heavy atoms for the H-bonds, all heavy for the contact.
POLAR_ONLY = {"N6.55", "E-ECL2"}
# TM6 fingerprint used to prove the asparagine really is 6.55: the spacing
# 6.48 -> 6.51 -> 6.52 -> 6.55 is 3, 1, 3 and no other asparagine in the pocket
# sits at the end of that run.
TM6_SPACING = {-7: "TRP", -4: "LEU", -3: ("HIE", "HIS", "HID", "HIP")}

GMX = "/opt/gromacs/2026.3-cuda13.1-sm120/bin/gmx"
WSL_DISTRO = "A2A-MD"
NON_PROTEIN = {"WAT", "SOL", "HOH", "POPC", "PA", "PC", "OL", "CHL",
               "Na+", "Cl-", "NA", "CL", "K+", "ZMA", "LIG", "MOL", "UNL"}
LIGAND_NAMES = ("LIG", "ZMA", "MOL", "UNL")
CANDIDATES = ["NPLIB01185", "NPLIB01643", "NPLIB03056"]
REPS = [1, 2, 3]


def _root() -> Path:
    # Delegates to kdss_paths, which finds whichever layout this file is inside
    # -- the development tree or a deposit -- by walking up from its own
    # location. It previously tried two absolute paths and exited if neither
    # existed, which on the author's machine meant a copy of the deposit read
    # and wrote the development tree instead of itself.
    return project_root()


ROOT = _root()
RUNS = dev_dir("20_new_runs")
CAND_BASE = under("20_new_runs", "md_candidates_2026-08-22")
CONTROL = under("20_new_runs", "md_production_2026-08-21/04_run")


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def read_gro(path: Path):
    """(resid, resname, atomname, atom_index_1based, x, y, z) in angstrom."""
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    n = int(lines[1])
    out = []
    for i in range(n):
        s = lines[2 + i]
        out.append((int(s[0:5]), s[5:10].strip(), s[10:15].strip(), i + 1,
                    float(s[20:28]) * 10, float(s[28:36]) * 10,
                    float(s[36:44]) * 10))
    return out


def residue_instances(atoms):
    """Protein residues in sequence order. A residue *instance* is a run of
    atoms sharing (resid, resname); resid alone is not unique in this topology."""
    inst, prev = [], None
    for ri, rn, an, idx, x, y, z in atoms:
        if rn in NON_PROTEIN:
            continue
        key = (ri, rn)
        if key != prev:
            inst.append({"resid": ri, "resname": rn, "atoms": []})
            prev = key
        inst[-1]["atoms"].append((an, idx, (x, y, z)))
    return inst


def ligand_atoms(atoms, polar_only=False):
    name = next((n for n in LIGAND_NAMES
                 if any(a[1] == n for a in atoms)), None)
    if name is None:
        raise SystemExit("no ligand residue found")
    out = []
    for ri, rn, an, idx, x, y, z in atoms:
        if rn != name or an.startswith("H"):
            continue
        if polar_only and an[0] not in "NO":
            continue
        out.append((an, idx, (x, y, z)))
    return name, out


def locate_anchors(inst, lig_all):
    """Find the three anchors in the control, by contact, then prove the
    asparagine is 6.55 from the TM6 spacing rather than from its number."""
    found = {}
    for tag, (resname, sc_atoms) in ANCHOR_ATOMS.items():
        best = None
        for pos, r in enumerate(inst):
            if r["resname"] != resname:
                continue
            sel = [a for a in r["atoms"] if a[0] in sc_atoms]
            if not sel:
                continue
            d = min(math.dist(c, lc) for _, _, c in sel for _, _, lc in lig_all)
            if best is None or d < best[0]:
                best = (d, pos)
        if best is None:
            raise SystemExit(f"no {resname} found for anchor {tag}")
        found[tag] = {"position": best[1], "start_distance_A": round(best[0], 2)}

    pos = found["N6.55"]["position"]
    for off, want in TM6_SPACING.items():
        got = inst[pos + off]["resname"]
        want_t = (want,) if isinstance(want, str) else want
        if got not in want_t:
            raise SystemExit(
                f"TM6 spacing check FAILED: position {pos}{off:+d} is {got}, "
                f"expected one of {want_t}. The asparagine at sequence position "
                f"{pos} is not 6.55 and every downstream number would be wrong.")
    found["N6.55"]["tm6_spacing_verified"] = (
        f"{inst[pos-7]['resname']}(6.48) {inst[pos-4]['resname']}(6.51) "
        f"{inst[pos-3]['resname']}(6.52) {inst[pos]['resname']}(6.55), "
        "spacing 3-1-3")
    return found


def write_index(path: Path, groups: dict[str, list[int]]):
    txt = []
    for name, idxs in groups.items():
        txt.append(f"[ {name} ]")
        for k in range(0, len(idxs), 15):
            txt.append(" ".join(str(i) for i in idxs[k:k + 15]))
    path.write_text("\n".join(txt) + "\n", encoding="utf-8", newline="\n")


def run_mindist(tpr: Path, xtc: Path, ndx: Path, group: str, out_xvg: Path,
                log: Path) -> bool:
    if out_xvg.is_file() and out_xvg.stat().st_size > 0:
        return True
    rel = lambda p: p.relative_to(ROOT).as_posix()
    cmd = (f"cd {wsl_root()!r} && printf '{group}\\nLigand\\n' | "
           f"{GMX} mindist -s '{rel(tpr)}' -f '{rel(xtc)}' -n '{rel(ndx)}' "
           f"-od '{rel(out_xvg)}' -tu ns > '{rel(log)}' 2>&1")
    subprocess.run(["wsl", "-d", WSL_DISTRO, "--", "bash", "-c", cmd],
                   capture_output=True, text=True, timeout=7200)
    return out_xvg.is_file() and out_xvg.stat().st_size > 0


def read_xvg(p: Path):
    t, v = [], []
    for line in p.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = line.strip()
        if not line or line[0] in "#@":
            continue
        parts = line.split()
        if len(parts) >= 2:
            t.append(float(parts[0]))            # ns; -tu ns was passed
            v.append(float(parts[1]) * 10.0)     # nm -> angstrom
    return t, v


def persistence(t, v, cutoff, t_min=0.0):
    sel = [(tt, vv) for tt, vv in zip(t, v) if tt >= t_min]
    if not sel:
        return None
    hits = sum(1 for _, vv in sel if vv <= cutoff)
    return {"n_frames": len(sel), "n_satisfying": hits,
            "persistence": round(hits / len(sel), 4),
            "median_A": round(sorted(vv for _, vv in sel)[len(sel) // 2], 2),
            "min_A": round(min(vv for _, vv in sel), 2)}


def analyse_system(tag, gro: Path, tpr: Path, xtc: Path, workdir: Path,
                   anchor_pos: dict):
    atoms = read_gro(gro)
    inst = residue_instances(atoms)
    lig_name, lig_all = ligand_atoms(atoms)
    _, lig_polar = ligand_atoms(atoms, polar_only=True)
    workdir.mkdir(parents=True, exist_ok=True)

    result = {"ligand_resname": lig_name, "n_ligand_heavy": len(lig_all),
              "n_ligand_polar": len(lig_polar), "anchors": {}}
    for tag_a, pos in anchor_pos.items():
        resname, sc_atoms = ANCHOR_ATOMS[tag_a]
        r = inst[pos]
        if r["resname"] != resname:
            raise SystemExit(
                f"{tag}: sequence position {pos} is {r['resname']}, expected "
                f"{resname}. The receptor differs between systems and no "
                f"anchor mapping is safe.")
        a_idx = [a[1] for a in r["atoms"] if a[0] in sc_atoms]
        if not a_idx:
            raise SystemExit(f"{tag}: {resname} at {pos} has none of {sc_atoms}")
        lig = lig_polar if tag_a in POLAR_ONLY else lig_all
        ndx = workdir / f"{tag_a.replace('.', '_')}.ndx"
        write_index(ndx, {"Anchor": a_idx, "Ligand": [a[1] for a in lig]})
        xvg = workdir / f"{tag_a.replace('.', '_')}.xvg"
        log = workdir / f"{tag_a.replace('.', '_')}.log"
        if not run_mindist(tpr, xtc, ndx, "Anchor", xvg, log):
            result["anchors"][tag_a] = {"available": False,
                                        "reason": "gmx mindist produced nothing"}
            continue
        t, v = read_xvg(xvg)
        cut = CUTOFF_FOR[tag_a]
        result["anchors"][tag_a] = {
            "available": True, "resname": resname, "resid_in_topology": r["resid"],
            "cutoff_A": cut,
            "full_trajectory": persistence(t, v, cut),
            "last_80ns": persistence(t, v, cut, t_min=SECONDARY_WINDOW_NS),
            "series_ns": [round(t[0], 3), round(t[-1], 3)],
        }
    # secondary: both clamp contacts in the same frame
    a, b = result["anchors"].get("N6.55"), result["anchors"].get("E-ECL2")
    if a and b and a.get("available") and b.get("available"):
        ta, va = read_xvg(workdir / "N6_55.xvg")
        tb, vb = read_xvg(workdir / "E-ECL2.xvg")
        n = min(len(va), len(vb))
        both = sum(1 for i in range(n)
                   if va[i] <= HBOND_CUTOFF_A and vb[i] <= HBOND_CUTOFF_A)
        result["bidentate_clamp"] = {
            "definition": "both N6.55 and E-ECL2 within 3.5 A in the same frame",
            "n_frames": n, "n_satisfying": both,
            "persistence": round(both / n, 4) if n else None}
    return result


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(under("20_new_runs", "md_recognition_motif_2026-09-02")))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    # 1. Locate the anchors in the control, and prove the asparagine is 6.55.
    ctrl_gro = CONTROL / "em.gro"
    atoms = read_gro(ctrl_gro)
    inst = residue_instances(atoms)
    _, lig_all = ligand_atoms(atoms)
    anchors = locate_anchors(inst, lig_all)
    print("anchors located in the ZM241385 control:")
    for tag, a in anchors.items():
        rn = inst[a["position"]]["resname"]
        print(f"  {tag:<8} {rn} at sequence position {a['position']:>4} "
              f"(topology resid {inst[a['position']]['resid']}), "
              f"{a['start_distance_A']:.2f} A in the starting structure")
    print(f"  spacing check: {anchors['N6.55']['tm6_spacing_verified']}\n")
    anchor_pos = {t: a["position"] for t, a in anchors.items()}

    systems = {"control_ZM241385": (ctrl_gro, CONTROL / "prod.tpr",
                                    CONTROL / "prod.xtc",
                                    out / "control_ZM241385")}
    for cid in CANDIDATES:
        run = CAND_BASE / cid / "04_run"
        for r in REPS:
            systems[f"{cid}_r{r}"] = (run / "em.gro", run / f"prod_r{r}.tpr",
                                      run / f"prod_r{r}.xtc",
                                      out / f"{cid}_r{r}")

    results = {}
    for tag, (gro, tpr, xtc, wd) in systems.items():
        if not (tpr.is_file() and xtc.is_file()):
            print(f"{tag}: trajectory missing, skipped")
            continue
        print(f"{tag} ... ", end="", flush=True)
        results[tag] = analyse_system(tag, gro, tpr, xtc, wd, anchor_pos)
        p = results[tag]["anchors"].get("N6.55", {})
        ft = p.get("full_trajectory") or {}
        print(f"N6.55 persistence {ft.get('persistence')}")

    # 2. Outcome, against criteria fixed in ADR-AIIM-161 before this ran.
    ctrl = results.get("control_ZM241385", {}).get("anchors", {}) \
                  .get("N6.55", {}).get("full_trajectory", {})
    p_ctrl = ctrl.get("persistence")
    outcomes = {}
    if p_ctrl:
        for cid in CANDIDATES:
            ps = [results[f"{cid}_r{r}"]["anchors"]["N6.55"]["full_trajectory"]
                  ["persistence"] for r in REPS if f"{cid}_r{r}" in results]
            if not ps:
                continue
            mean_p = sum(ps) / len(ps)
            rho = mean_p / p_ctrl
            outcomes[cid] = {
                "replicate_persistence": ps,
                "mean_persistence": round(mean_p, 4),
                "rho_vs_control": round(rho, 3),
                "outcome": ("REPRODUCES" if rho >= RHO_REPRODUCES else
                            "DOES_NOT_REPRODUCE" if rho < RHO_DOES_NOT
                            else "PARTIAL"),
            }

    report = {
        "question": "Does the predicted pose hold the interaction that defines "
                    "antagonist binding at A2A, over 100 ns?",
        "authority": "ADR-AIIM-161, pre-registered before the trajectories were read",
        "skill": {"name": "molecular-dynamics", "version": "1.0",
                  "influence": "contact-analysis framing and the "
                               "equilibration-discarded secondary window; its "
                               "OpenMM/MDAnalysis code paths were not used"},
        "criteria": {
            "hbond_cutoff_A": HBOND_CUTOFF_A,
            "contact_cutoff_A": CONTACT_CUTOFF_A,
            "primary_endpoint": "N6.55 H-bond persistence, full trajectory",
            "rho_reproduces_at_or_above": RHO_REPRODUCES,
            "rho_does_not_reproduce_below": RHO_DOES_NOT,
            "no_statistical_test": "the control is one trajectory; rho is "
                                   "descriptive",
        },
        "anchor_identification": {
            "method": "located by contact in the ZM241385 control, then proved "
                      "by TM6 helical spacing 6.48-6.51-6.52-6.55 = 3-1-3",
            "why_not_resid": "residue numbers are not unique in this topology "
                             "(resid 36 is both LEU and ASN); every group is "
                             "built from explicit atom indices",
            "anchors": {t: {**a, "resname": inst[a["position"]]["resname"]}
                        for t, a in anchors.items()},
        },
        "control": {
            "ligand": "ZM241385, the co-crystal antagonist",
            "n_replicates": 1,
            "tpr_sha256": sha256(CONTROL / "prod.tpr"),
            "N6.55_persistence": p_ctrl,
        },
        "by_system": results,
        "outcomes": outcomes,
    }
    dest = out / "recognition_motif.json"
    dest.write_text(json.dumps(report, indent=2), encoding="utf-8")

    print(f"\ncontrol ZM241385 N6.55 persistence: {p_ctrl}")
    print(f"{'compound':<13}{'r1':>7}{'r2':>7}{'r3':>7}{'mean':>8}{'rho':>7}  outcome")
    for cid, o in outcomes.items():
        ps = o["replicate_persistence"] + [None] * (3 - len(o["replicate_persistence"]))
        cells = "".join(f"{p:>7.3f}" if p is not None else f"{'-':>7}" for p in ps)
        print(f"{cid:<13}{cells}{o['mean_persistence']:>8.3f}"
              f"{o['rho_vs_control']:>7.2f}  {o['outcome']}")
    print(f"\nwritten -> {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
