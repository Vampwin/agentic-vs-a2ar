"""Figure of the post hoc RMSF profiles (ADR-AIIM-394): receptor C-alpha RMSF per residue, three runs per system, site residues marked."""
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

H = Path(__file__).resolve().parent
d = json.load(open(H / "posthoc_rmsf_hbond.json"))
res = d["_residues"]
x = np.arange(len(res))
fig, axs = plt.subplots(3, 1, figsize=(10, 7.2), sharex=True, sharey=True)
cols = {"r1": "#1f77b4", "r2": "#d62728", "r3": "#2ca02c"}
for ax, (sysn, lab) in zip(axs, (("CONTROL", "ZM241385 control"), ("SNP00157", "SNP00157"), ("SNP00163", "SNP00163"))):
    site = set(d[sysn]["site_residues"])
    for i, r in enumerate(res):
        if r in site:
            ax.axvspan(i - 0.5, i + 0.5, color="#bbbbbb", alpha=0.6, lw=0)
    for r in ("r1", "r2", "r3"):
        ax.plot(x, d[sysn]["rmsf"][r]["profile_A"], color=cols[r], lw=0.9, label=r)
    ax.set_ylabel("C-alpha RMSF (A)")
    ax.set_title(lab + "  (grey: residues within 4.5 A of the ligand in the starting structure)", fontsize=9, loc="left")
    ax.legend(ncol=3, fontsize=8, loc="upper right", frameon=False)
    ax.set_ylim(0, 4)
chb = next(i for i, r in enumerate(res) if r.startswith("B:"))
for ax in axs:
    ax.axvline(chb - 0.5, color="k", lw=0.6, ls=":")
axs[-1].set_xlabel("receptor residue, chain A then chain B (dotted line); numbers as in the simulation structure; y axis cut at 4 A")
fig.tight_layout()
fig.savefig(H / "posthoc_rmsf_profiles.png", dpi=160)
print("ok", max(max(d[s]["rmsf"][r]["profile_A"]) for s in ("CONTROL", "SNP00157", "SNP00163") for r in ("r1", "r2", "r3")))
