"""(ADR-AIIM-396) Replicate 3 only: the same script as pocket_start_snplib.py, unchanged except for the replicate list (3,).

Pocket-from-start distance for SNP00157 and SNP00163, as md_17_pocket_from_start.py makes it for the
NP-library compounds, so that panel A of the molecular-dynamics figure can show all five compounds
by one measure.

Same definition: the binding site is every receptor C-alpha within 10 A of the ligand in the
energy-minimised starting structure (em.gro); the series is the gmx mindist between that group and the
ligand, minimum-image applied by gmx itself. The frozen md_17 is imported for its two helpers
(read_gro, pocket_indices) and is not modified; its main() is NP-specific and is not used, because it
loops over three fixed compounds and its first shell string lacks an f-prefix. Output goes where md_17's
goes, `analysis_prod_r<k>/pocket_start_dist.xvg` in each compound's 04_run.
"""
import subprocess, sys, hashlib, json
from pathlib import Path
sys.path.insert(0, r"D:\Agentic VS Workflow\10_new_system\scripts")
import md_17_pocket_from_start as m17

ROOT = Path(r"D:\Agentic VS Workflow")
BASE = ROOT / "20_new_runs" / "md_candidates_2026-08-22"
GMX = m17.GMX
OUT = ROOT / "20_new_runs" / "md_snplib_replicate3_2026-10-10"
WSL_ROOT = "/mnt/d/Agentic VS Workflow"
res = {}
for cid in ("SNP00157", "SNP00163"):
    run = BASE / cid / "04_run"
    atoms = m17.read_gro(run / "em.gro")
    pocket = m17.pocket_indices(atoms, 1.0)                  # 10 A in nm
    lig = [a[0] for a in atoms if a[1] == "LIG"]
    print(f"{cid}: pocket = {len(pocket)} C-alpha within 10 A of the docked ligand ({len(lig)} atoms)")
    for r in (3,):
        d = run / f"analysis_prod_r{r}"
        ndx = d / "pocket_start.ndx"
        ndx.write_text("[ Pocket ]\n" + "\n".join(" ".join(map(str, pocket[k:k + 15])) for k in range(0, len(pocket), 15))
                       + "\n[ Ligand ]\n" + "\n".join(" ".join(map(str, lig[k:k + 15])) for k in range(0, len(lig), 15)) + "\n",
                       encoding="utf-8", newline="\n")
        xvg = d / "pocket_start_dist.xvg"
        rel = lambda p: p.relative_to(ROOT).as_posix()
        cmd = (f"cd '{WSL_ROOT}' && export OMP_NUM_THREADS=4 GMX_MAXBACKUP=-1 && printf 'Pocket\nLigand\n' | "
               f"{GMX} mindist -s '{rel(run / f'prod_r{r}.tpr')}' -f '{rel(run / f'prod_r{r}.xtc')}' "
               f"-n '{rel(ndx)}' -od '{rel(xvg)}' -tu ns > '{rel(d / 'pocket_start.log')}' 2>&1")
        subprocess.run(["wsl", "-d", "A2A-MD", "--", "bash", "-c", cmd], capture_output=True, text=True, timeout=3600)
        v = [float(l.split()[1]) * 10 for l in xvg.read_text().splitlines() if l and l[0] not in "@#"]
        t = [float(l.split()[0]) for l in xvg.read_text().splitlines() if l and l[0] not in "@#"]
        res[f"{cid}_r{r}"] = {"frames": len(v), "t_end_ns": round(t[-1], 2), "min_A": round(min(v), 2),
                              "median_A": round(sorted(v)[len(v) // 2], 2), "max_A": round(max(v), 2),
                              "frames_beyond_10A": sum(x > 10 for x in v),
                              "sha256": hashlib.sha256(xvg.read_bytes()).hexdigest()}
        print(f"  r{r}: {res[f'{cid}_r{r}']}")
(OUT / "pocket_start_r3.json").write_text(json.dumps(res, indent=1), "utf-8")
