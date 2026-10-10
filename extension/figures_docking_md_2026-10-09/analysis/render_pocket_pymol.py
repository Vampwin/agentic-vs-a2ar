"""Panels A, B and C of the redesigned Figure 6 (ADR-AIIM-366): ZM241385, SNP00157 and SNP00163 in the site, run inside PyMOL.

Each panel shows the crystallographic ZM241385 pose (thin grey, reference), the Vina rank-1 pose from the 22 A site box
(blue) and the DiffDock-Pocket rank-1 pose of the main run (red). One camera for the three panels: orientation set on the
crystal pose, then turned -15 degrees about x as in the earlier renders, orthographic, the same zoom on the crystal pose, so
the three panels are at the same scale and the same view. Receptor 4EIY_protein_plus_waters (waters removed), grey cartoon
(transparency 0.35) drawn as WHOLE secondary-structure segments (helix, strand or loop) of which any residue is within 9 A of the
crystal pose, so that every ribbon is continuous (the same set in every panel); the Asn253 side chain is drawn and its label is
added in the composite (matplotlib, the figure font), at the position calibrated here (_pocket_calib_<id>.json).
Coordinates are loaded verbatim; nothing is re-fit. Vina poses are model 1 of the ten written per compound.

    pymol -cq render_pocket_pymol.py
"""
import glob
from pymol import cmd

ROOT = "/mnt/d/Agentic VS Workflow"
RUNS = f"{ROOT}/20_new_runs"
RECEPTOR = f"{RUNS}/snp_screen_2026-08-04/G2_preparation/4EIY_protein_plus_waters.pdb"
ZMA_CRYSTAL = f"{RUNS}/snp_screen_2026-08-04/G2_preparation/ZMA_reference_4EIY.sdf"
POCKET = f"{RUNS}/docking_matched_filter_pocket_2026-10-09/pocket_C"
OUT = f"{RUNS}/figures_docking_md_2026-10-09/source"
POCKET_RADIUS = 9   # a residue segment is drawn whole if any of its residues is within this distance of the crystal pose
CARTOON_GREY, VINA_BLUE, DD_RED, REF_GREY = "0xc4cad4", "0x2c5f8a", "0x9b2226", "0x7a7a7a"   # the palette of figstyle.py
ZOOM_BUFFER = 11


def setup():
    cmd.reinitialize()
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 0)
    cmd.set("antialias", 2)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("ray_trace_mode", 0)
    cmd.set("depth_cue", 0)
    cmd.set("specular", 0.1)
    cmd.set("valence", 1)
    cmd.set("orthoscopic", 1)
    cmd.load(RECEPTOR, "receptor")
    cmd.remove("receptor and resn HOH")
    cmd.hide("everything", "receptor")
    cmd.show("cartoon", "receptor")
    cmd.color(CARTOON_GREY, "receptor")
    cmd.set("cartoon_transparency", 0.35, "receptor")


def lig(path, name, color, radius, first_model_only=False):
    if first_model_only:
        cmd.load(path, "_tmp_all")
        cmd.create(name, "_tmp_all", 1, 1)
        cmd.delete("_tmp_all")
    else:
        cmd.load(path, name)
    cmd.remove(f"{name} and hydro")
    assert cmd.count_states(name) == 1 and cmd.count_atoms(name) > 0, name
    cmd.hide("everything", name)
    cmd.show("sticks", name)
    cmd.color(color, name)
    cmd.set("stick_radius", radius, name)


def asn253_stick():
    sel = "receptor and resi 253 and name CA"
    if cmd.count_atoms(sel) == 0:
        raise SystemExit("Asn253 CA not found")
    cmd.show("sticks", "receptor and resi 253 and sidechain")
    cmd.color("0x6b7280", "receptor and resi 253 and sidechain")
    cmd.set("stick_radius", 0.07, "receptor and resi 253")


def whole_segments(radius):
    """Residues of every secondary-structure segment that has a residue within `radius` of the crystal pose."""
    near = set()
    cmd.iterate(f"byres (receptor within {radius} of zma) and name CA", "near.add(int(resi))", space={"near": near})
    cas = []
    cmd.iterate("receptor and name CA", "cas.append((int(resi), ss))", space={"cas": cas})
    cas.sort()
    keep, seg, prev = set(), [], None
    for resi, ss in cas + [(10**9, "#")]:
        if prev is not None and (resi != prev[0] + 1 or ss != prev[1]):
            if any(r in near for r in seg):
                keep.update(seg)
            seg = []
        seg.append(resi); prev = (resi, ss)
    return sorted(keep)


def calibrate(tag):
    """Where the Asn253 C-alpha falls in the image: render a lone marker with the same camera and find it."""
    sel = "receptor and resi 253 and name CA"
    cmd.pseudoatom("calib", pos=cmd.get_atom_coords(sel))
    cmd.hide("everything", "all")
    cmd.show("spheres", "calib"); cmd.set("sphere_scale", 1.2, "calib"); cmd.color("magenta", "calib")
    cmd.ray(2400, 2100)
    cmd.png(f"{OUT}/_pocket_calib_{tag}.png", width=2400, height=2100, dpi=72)


def dd_pocket_pose(name):
    d = glob.glob(f"{POCKET}/index*___4EIY_protein_only.pdb___{name}.sdf")
    assert len(d) == 1, (name, d)
    return f"{d[0]}/rank1.sdf"


def render(tag, vina_path, vina_first_model):
    setup()
    lig(ZMA_CRYSTAL, "zma", REF_GREY, 0.09)
    lig(vina_path, "vina", VINA_BLUE, 0.20, first_model_only=vina_first_model)
    lig(dd_pocket_pose(tag), "dd", DD_RED, 0.20)
    cmd.order("vina dd zma", location="top")
    cmd.dss("receptor")   # secondary structure from the coordinates, so that segments are real helices/strands/loops
    keep = whole_segments(POCKET_RADIUS)
    cmd.hide("cartoon", "receptor")
    cmd.show("cartoon", "receptor and resi " + "+".join(str(r) for r in keep))
    asn253_stick()
    cmd.orient("zma")
    cmd.turn("x", -15)
    cmd.zoom("zma or (receptor and resi 253 and name CA)", ZOOM_BUFFER)
    cmd.ray(2400, 2100)
    cmd.png(f"{OUT}/_pocket_panel_{tag}.png", width=2400, height=2100, dpi=600)
    calibrate(tag)
    print("RENDERED", tag, "segments residues", len(keep))


render("ZM241385", f"{RUNS}/snplib_dock_2026-10-06/poses_zma/ZM241385.sdf", True)
for snp in ("SNP00157", "SNP00163"):
    render(snp, f"{RUNS}/snplib_dock_2026-10-06/poses/{snp}.sdf", True)
print("PYMOL_POCKET_PANELS_DONE")
