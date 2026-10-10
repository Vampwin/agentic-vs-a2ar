"""Panels A, B and C of Figure 6: SNP00157 and SNP00163, Vina rank-1 against DiffDock rank-1, run inside PyMOL.

Same conventions as 10_new_system/scripts/render_pose_overlay_pymol.py, whose conventions are reused, with one camera for all three panels: receptor 4EIY_protein_plus_waters.pdb, grey cartoon (transparency 0.35) restricted to residues within 13 A of
the displayed poses, one flat colour per pose (Vina blue, DiffDock red, stick radius 0.20), Asn253 labelled,
orient on the crystal pose, turn x -15, zoom(buffer 9), ray 2400 x 2100. Coordinates are loaded verbatim; nothing is re-fit.
The crystallographic ZM241385 pose is added as a thin grey reference so that the reader can see where the pocket is;
the zoom covers both poses plus that reference. Vina poses are model 1 of the ten written per compound (the rank-1 pose that was scored).

    pymol -cq render_snplib_pymol.py
"""
from pymol import cmd

ROOT = "/mnt/d/Agentic VS Workflow"
RUNS = f"{ROOT}/20_new_runs"
RECEPTOR = f"{RUNS}/snp_screen_2026-08-04/G2_preparation/4EIY_protein_plus_waters.pdb"
ZMA_CRYSTAL = f"{RUNS}/snp_screen_2026-08-04/G2_preparation/ZMA_reference_4EIY.sdf"
OUT = f"{RUNS}/figures_docking_md_2026-10-09/source"
POCKET_RADIUS = 9   # a ribbon segment is drawn whole if any of its residues is within this distance of the displayed poses
CARTOON_GREY, VINA_BLUE, DIFFDOCK_RED = "0xc4cad4", "0x2c5f8a", "0x9b2226"   # the palette of figstyle.py


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
    cmd.set("orthoscopic", 1)   # no perspective, so the same camera gives the same view at any zoom
    cmd.load(RECEPTOR, "receptor")
    cmd.remove("receptor and resn HOH")
    cmd.hide("everything", "receptor")
    cmd.show("cartoon", "receptor")
    cmd.color(CARTOON_GREY, "receptor")
    cmd.set("cartoon_transparency", 0.35, "receptor")


def lig(path, name, color, first_model_only=False):
    if first_model_only:
        cmd.load(path, "_tmp_all")
        cmd.create(name, "_tmp_all", 1, 1)
        cmd.delete("_tmp_all")
    else:
        cmd.load(path, name)
    assert cmd.count_states(name) == 1 and cmd.count_atoms(name) > 0, name
    cmd.hide("everything", name)
    cmd.show("sticks", name)
    cmd.color(color, name)
    cmd.set("stick_radius", 0.20, name)


def asn253_stick():
    sel = "receptor and resi 253 and name CA"
    if cmd.count_atoms(sel) == 0:
        raise SystemExit("Asn253 CA not found")
    cmd.show("sticks", "receptor and resi 253 and sidechain")
    cmd.color("0x6b7280", "receptor and resi 253 and sidechain")
    cmd.set("stick_radius", 0.07, "receptor and resi 253")


def whole_segments(radius, around):
    """Residues of every secondary-structure segment that has a residue within `radius` of the selection `around`."""
    near = set()
    cmd.iterate(f"byres (receptor within {radius} of ({around})) and name CA", "near.add(int(resi))", space={"near": near})
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
    """Where the Asn253 C-alpha falls in the image: a lone marker rendered with the same camera (the label is drawn in the composite)."""
    cmd.pseudoatom("calib", pos=cmd.get_atom_coords("receptor and resi 253 and name CA"))
    cmd.hide("everything", "all")
    cmd.show("spheres", "calib"); cmd.set("sphere_scale", 1.2, "calib"); cmd.color("magenta", "calib")
    cmd.ray(2400, 2100)
    cmd.png(f"{OUT}/_calib_{tag}.png", width=2400, height=2100, dpi=72)


ZMA_VINA = f"{RUNS}/snp_screen_2026-08-04/G3_redocking/ZMA_redocked_4EIY.sdf"
ZMA_DD = f"{RUNS}/diffdock_orthogonal_2026-08-23/results/ZM241385_control/rank1.sdf"


def render(tag, vina_path, dd_path, first_model_only):
    """All three panels share one camera orientation: orient on the crystallographic ZM241385 pose alone, then
    turn x -15; only the zoom differs, so that each panel fits its own poses."""
    setup()
    lig(vina_path, "vina", VINA_BLUE, first_model_only=first_model_only)
    lig(dd_path, "dd", DIFFDOCK_RED)
    lig(ZMA_CRYSTAL, "zma", "0x7a7a7a")
    cmd.set("stick_radius", 0.09, "zma")
    cmd.order("vina dd zma", location="top")
    ALL = "vina or dd or zma"
    cmd.dss("receptor")
    keep = whole_segments(POCKET_RADIUS, ALL)
    cmd.hide("cartoon", "receptor")
    cmd.show("cartoon", "receptor and resi " + "+".join(str(r) for r in keep))
    asn253_stick()
    cmd.orient("zma")
    cmd.turn("x", -15)
    cmd.zoom(ALL + " or (receptor and resi 253 and name CA)", 11)
    cmd.ray(2400, 2100)
    cmd.png(f"{OUT}/_panel_{tag}.png", width=2400, height=2100, dpi=600)
    calibrate(tag)
    print("RENDERED", tag)


render("control", ZMA_VINA, ZMA_DD, True)
for snp in ("SNP00157", "SNP00163"):
    render(snp, f"{RUNS}/snplib_dock_2026-10-06/poses/{snp}.sdf", f"{RUNS}/diffdock_snplib_2026-10-07/results/{snp}/rank1.sdf", True)
print("PYMOL_SNP_PANELS_DONE")
