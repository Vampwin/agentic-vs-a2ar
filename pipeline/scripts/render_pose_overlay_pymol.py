"""
render_pose_overlay_pymol.py -- the pose-agreement figure, run inside PyMOL.

Two panels, not one pose each. The point is not "here is a docked ligand" --
every other figure in this manuscript avoids depicting a result it cannot
back, and a single confident-looking pose would assert a binding mode
Section 3.5 explicitly says cannot be distinguished between engines. So:

  A  Positive control, ZM241385. Crystallographic pose (thin grey reference)
     against the Vina redock (bold blue, 0.379 A -- near-total overlap) and
     the DiffDock pose (bold red, 1.92 A -- visibly offset but still within
     the 2.0 A gate). This is what agreement looks like when both engines
     are competent.

  B  NPLIB01643. Vina rank-1 (blue) against DiffDock rank-1 (red), 8.80 A
     apart -- the largest disagreement of the three candidates, and the one
     where DiffDock's confidence (+0.47, "high") was nonetheless highest.
     Same pocket, different molecule entirely. Drawn at the same scale as
     panel A so the two are directly comparable.

Every coordinate is loaded verbatim from the frozen artifacts already
verified against the manuscript's own RMSD figures (0.379, 1.92, 8.80 A) --
nothing here is re-fit or adjusted for a better-looking overlay.

Three rendering choices, and why.

Colour is one flat colour per pose, not CPK. Standard element colouring
paints every nitrogen blue and every oxygen red on all three molecules at
once, which makes overlapping poses blend into one indistinguishable mass.
One flat colour per pose lets the eye follow a single trace through the
overlap -- the standard convention in docking-comparison figures.

The reference (crystal) pose is drawn thin, the compared poses bold. Vina's
near-total overlap with the crystal pose (0.379 A) simply vanished under a
crystal stick of comparable thickness in an earlier attempt. A thin
reference with bold prediction on top shows overlap as a bold tube with a
thin grey pencil-outline at its edge -- visible agreement rather than an
invisible one.

The cartoon is restricted to the pocket neighbourhood rather than the whole
receptor. Three small black marks appeared in early renders of the full
receptor; chased down rather than papered over, they turned out to be real:
two unresolved loop stretches (missing electron density -- routine in GPCR
crystal structures) and the boundary left by excising the BRIL fusion
domain, confirmed against the residue numbering itself
(gaps at 208->219, 307->1001, 1042->1058) and independently ruled out as
disulfide-bond rendering (real Cys-Cys pairs exist here, but unbonding them
changed nothing). None of the three sits near the orthosteric site, so a
local-pocket view is not a concession to a rendering artifact -- it is the
right framing for a pose-comparison figure on its own terms, and it happens
to exclude all three discontinuities along with the empty helix bulk a
whole-receptor view carried.

Run inside PyMOL (not as a plain Python script):
    pymol -cq render_pose_overlay_pymol.py
"""
from pymol import cmd

ROOT = "/mnt/d/Agentic VS Workflow"
RECEPTOR = f"{ROOT}/20_new_runs/snp_screen_2026-08-04/G2_preparation/4EIY_protein_plus_waters.pdb"
ZMA_CRYSTAL = f"{ROOT}/20_new_runs/snp_screen_2026-08-04/G2_preparation/ZMA_reference_4EIY.sdf"
ZMA_VINA = f"{ROOT}/20_new_runs/snp_screen_2026-08-04/G3_redocking/ZMA_redocked_4EIY.sdf"
ZMA_DIFFDOCK = f"{ROOT}/20_new_runs/diffdock_orthogonal_2026-08-23/results/ZM241385_control/rank1.sdf"
NPLIB_VINA = f"{ROOT}/20_new_runs/np_library_funnel_2026-08-15/full_poses_chunk1/NPLIB01643.sdf"
NPLIB_DIFFDOCK = f"{ROOT}/20_new_runs/diffdock_orthogonal_2026-08-23/results/NPLIB01643/rank1.sdf"
OUT = f"{ROOT}/30_manuscript/figures"

POCKET_RADIUS = 13  # Angstrom, around all loaded ligand poses in a panel

CARTOON_GREY = "0xe4e7ec"
REF_GREY = "0x7a7a7a"
VINA_BLUE = "0x1f5fa8"
DIFFDOCK_RED = "0xd6403f"


def setup_common():
    cmd.reinitialize()
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 0)
    cmd.set("antialias", 2)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("ray_trace_mode", 0)
    cmd.set("depth_cue", 0)
    cmd.set("specular", 0.1)
    cmd.set("valence", 1)


def load_receptor():
    cmd.load(RECEPTOR, "receptor")
    cmd.remove("receptor and resn HOH")
    cmd.hide("everything", "receptor")
    cmd.show("cartoon", "receptor")
    cmd.color(CARTOON_GREY, "receptor")
    cmd.set("cartoon_transparency", 0.55, "receptor")


def restrict_cartoon_to_pocket(ligand_names: list[str]):
    """Local-pocket view; see module docstring for why this is the right
    framing, not just a way to dodge the three chain discontinuities."""
    around = " or ".join(f"({n})" for n in ligand_names)
    cmd.hide("cartoon",
             f"receptor and not byres (receptor within {POCKET_RADIUS} of ({around}))")


def load_ligand_flat(path, name, color, radius):
    """One flat colour for the whole molecule -- not CPK. See module docstring."""
    cmd.load(path, name)
    cmd.hide("everything", name)
    cmd.show("sticks", name)
    cmd.color(color, name)
    cmd.set("stick_radius", radius, name)


def frame_and_ray(width, height, zoom_sel, zoom_buffer, out_png):
    cmd.orient(zoom_sel)
    cmd.zoom(zoom_sel, zoom_buffer)
    cmd.turn("x", -15)
    cmd.ray(width, height)
    cmd.png(out_png, width=width, height=height, dpi=600)


def panel_a_control():
    setup_common()
    load_receptor()
    # thin reference first, bold predictions on top -- see module docstring
    load_ligand_flat(ZMA_CRYSTAL, "zma_crystal", REF_GREY, radius=0.09)
    load_ligand_flat(ZMA_VINA, "zma_vina", VINA_BLUE, radius=0.20)
    load_ligand_flat(ZMA_DIFFDOCK, "zma_diffdock", DIFFDOCK_RED, radius=0.20)
    cmd.order("zma_vina zma_diffdock zma_crystal", location="top")
    restrict_cartoon_to_pocket(["zma_crystal", "zma_vina", "zma_diffdock"])
    frame_and_ray(2400, 2100, "zma_crystal or zma_vina or zma_diffdock", 5,
                 f"{OUT}/_panelA_control.png")


def panel_b_candidate():
    setup_common()
    load_receptor()
    load_ligand_flat(NPLIB_VINA, "n1643_vina", VINA_BLUE, radius=0.20)
    load_ligand_flat(NPLIB_DIFFDOCK, "n1643_diffdock", DIFFDOCK_RED, radius=0.20)
    restrict_cartoon_to_pocket(["n1643_vina", "n1643_diffdock"])
    frame_and_ray(2400, 2100, "n1643_vina or n1643_diffdock", 5,
                 f"{OUT}/_panelB_candidate.png")


panel_a_control()
panel_b_candidate()
print("PYMOL_PANELS_DONE")
