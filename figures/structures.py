#!/usr/bin/env python
"""Structure renderings for Figures 1 and 2 (PyMOL, ray-traced, white background).

Fig 1 (same view, all superposed on Vasa RecA1)
  struct_vasa_phe   Vasa-AMPPNP (2DB3 A): Phe anchor F247, Q-motif stacker Y265,
                    Q-motif Pro P268, Gln Q272, adenine; dashes = anchor ring edge->N6, Gln->N6
  struct_prp5_trp   Prp5-ADP (4LJY A): Trp anchor W257, stacker F276, L279 (Q-4), Q283
  struct_mtr4_none  Mtr4-ADP (2XGJ A), Ski2-like: Q-motif Gln Q154 reads N6, but the anchor
                    slot is empty (sphere = where Vasa F247's ring sits); Mtr4's Q-25 Tyr129
                    lies ~27 A from N6
Fig 2
  struct_ddx19_pivot  DDX19B closed (3FHT A, AMPPNP+RNA) vs open (6B4K A, AMPPNP, N-terminal
                      helix), superposed on RecA1: residues before the anchor re-route, F94 stays
  struct_ejc_y14      eIF4AIII in the EJC (2J0Q A with its Y14, chain G): Y14 R68/W73/E122-Y124
                      contact D41 (UA+1) next to the buried anchor F40
Colours follow sc_style: Trp Crimson, Phe Deep Forest, stacker Sand, Gln Moss, Q-4 Sage,
nucleotide Terracotta, partner Sand, protein cartoon light grey.
"""
from pathlib import Path

from pymol import cmd

ROOT = Path(__file__).resolve().parents[1]
S0 = ROOT / "stage0" / "data"
OUT = ROOT / "figures" / "out"
C = {"crimson": "0x9B2335", "terracotta": "0xC96A3A", "sand": "0xDFB166", "sage": "0xB8C4A2",
     "moss": "0x5E7A52", "forest": "0x2F4A3E", "grey": "0xDCDCDC", "annot": "0x4A4A4A"}
W, H = 2000, 1600
TIP = {"PHE": "CZ", "TYR": "OH", "TRP": "CZ2", "PRO": "CG", "GLN": "CD", "LEU": "CG", "ASP": "CG",
       "ARG": "CZ", "GLU": "CD", "THR": "OG1"}


def setup():
    cmd.reinitialize()
    for k, v in C.items():
        cmd.set_color(f"sc_{k}", [int(v[2:4], 16) / 255, int(v[4:6], 16) / 255, int(v[6:8], 16) / 255])
    cmd.bg_color("white")
    for k, v in {"ray_opaque_background": 1, "antialias": 2, "ray_shadows": 0, "orthoscopic": 1,
                 "cartoon_transparency": 0.45, "cartoon_fancy_helices": 1, "stick_radius": 0.2,
                 "dash_gap": 0.25, "dash_radius": 0.07, "dash_color": "sc_annot", "label_size": 24,
                 "label_font_id": 5, "label_color": "black", "float_labels": 1, "depth_cue": 0,
                 "specular": 0.1, "ambient": 0.45, "sphere_transparency": 0.55}.items():
        cmd.set(k, v)


def load(name, pdb):
    path = next(p for p in (S0 / "cif" / f"{pdb}.cif.gz", S0 / "loop_cif" / f"{pdb}.cif.gz",
                            S0 / "partner_cif" / f"{pdb}.cif.gz") if p.exists())
    cmd.load(str(path), name)
    cmd.remove(f"{name} and (solvent or hydro)")


def sticks(sel, color):
    cmd.show("sticks", f"({sel}) and not name N+C+O")
    cmd.color(color, f"({sel}) and elem C")
    cmd.color("atomic", f"({sel}) and not elem C")


def label(sel, text, offset=(1.5, 1.0, 3.0)):
    """Label the side-chain tip, pushed toward the viewer and off the atoms."""
    m = cmd.get_model(f"({sel}) and name CA")
    atom = TIP.get(m.atom[0].resn, "CA") if m.atom else "CA"
    target = f"({sel}) and name {atom}"
    if not cmd.count_atoms(target):
        target = f"({sel}) and name CA"
    cmd.label(target, f'"{text}"')
    cmd.set("label_position", offset, target)


def nearest_dash(name, a_sel, b_sel):
    """One dash from b_sel to the nearest atom of a_sel."""
    best = min(((cmd.get_distance(f"({a_sel}) and name {n}", b_sel), n)
                for n in [at.name for at in cmd.get_model(a_sel).atom]), default=None)
    if best:
        cmd.distance(name, f"({a_sel}) and name {best[1]}", b_sel)
        cmd.hide("labels", name)
        return best[0]


def pocket(obj, anchor, stacker, pro, gln, lig, anchor_color, names, offsets):
    base = f"{obj} and chain A"
    cmd.hide("everything", obj)
    cmd.show("cartoon", base)
    cmd.color("sc_grey", base)
    if anchor:
        sticks(f"{base} and resi {anchor}", anchor_color)
    sticks(f"{base} and resi {stacker}", "sc_sand")
    sticks(f"{base} and resi {pro}", "sc_sage")
    sticks(f"{base} and resi {gln}", "sc_moss")
    # the nucleotide next to this chain's Q-motif Gln (ligand chain labels vary between entries)
    nuc = f"({obj} and resn {lig} and byres ({obj} within 6 of ({base} and resi {gln})))"
    cmd.show("sticks", nuc)
    cmd.color("sc_terracotta", f"{nuc} and elem C")
    n6 = f"{nuc} and name N6"
    out = {"gln": nearest_dash(f"{obj}_gln", f"{base} and resi {gln} and name OE1+NE2", n6)}
    if anchor:
        out["edge"] = nearest_dash(f"{obj}_edge", f"{base} and resi {anchor} and name CZ+CZ2+CE1+CE2+CH2+CZ3", n6)
    for resi, text in names.items():
        label(f"{base} and resi {resi}", text, offsets.get(resi, (1.5, 1.0, 3.0)))
    return out


def render(stem):
    OUT.mkdir(parents=True, exist_ok=True)
    cmd.png(str(OUT / f"{stem}.png"), width=W, height=H, dpi=600, ray=1)
    print("wrote", stem)


def fig1():
    setup()
    load("vasa", "2db3")
    load("prp5", "4ljy")
    load("mtr4", "2xgj")
    cmd.super("prp5 and chain A and resi 270-480 and name CA", "vasa and chain A and resi 262-470 and name CA")
    cmd.super("mtr4 and chain A and resi 140-340 and name CA", "vasa and chain A and resi 262-470 and name CA")
    focus = "vasa and chain A and (resi 247+265+268+272 or (resn ANP and name N9+C8+N7+C5+C6+N6+N1+C2+N3+C4))"
    cmd.orient(focus)
    cmd.turn("x", -20)
    cmd.zoom(focus, 5)
    view = cmd.get_view()
    near = lambda o: f"{o} and not ({o} within 14 of ({focus}))"

    d = pocket("vasa", 247, 265, 268, 272, "ANP", "sc_forest",
               {247: "F247 (anchor)", 265: "Y265", 268: "P268", 272: "Q272"},
               {247: (2.5, 2.0, 4.0), 265: (-2.5, 2.0, 4.0), 268: (2.5, -1.0, 4.0), 272: (2.5, -2.0, 4.0)})
    print("Vasa", {k: round(v, 1) for k, v in d.items()})
    cmd.hide("everything", "prp5 or mtr4")
    cmd.hide("cartoon", near("vasa"))
    cmd.set_view(view)
    render("struct_vasa_phe")
    # where Vasa's anchor ring sits, for the Mtr4 panel
    ring = cmd.get_coords("vasa and chain A and resi 247 and name CG+CD1+CD2+CE1+CE2+CZ").mean(0).tolist()

    cmd.hide("everything", "vasa")
    cmd.delete("vasa_*")
    d = pocket("prp5", 257, 276, 279, 283, "ADP", "sc_crimson",
               {257: "W257 (anchor)", 276: "F276", 279: "L279", 283: "Q283"},
               {257: (2.5, 2.0, 4.0), 276: (-2.5, 2.0, 4.0), 279: (2.5, -1.0, 4.0), 283: (2.5, -2.0, 4.0)})
    print("Prp5", {k: round(v, 1) for k, v in d.items()})
    cmd.hide("cartoon", near("prp5"))
    cmd.set_view(view)
    render("struct_prp5_trp")

    cmd.hide("everything", "prp5")
    cmd.delete("prp5_*")
    d = pocket("mtr4", None, 148, 150, 154, "ADP", "sc_sage",
               {148: "F148", 154: "Q154"}, {148: (-2.5, 2.0, 4.0), 154: (2.5, -2.0, 4.0)})
    cmd.pseudoatom("slot", pos=ring)
    cmd.show("spheres", "slot")
    cmd.set("sphere_scale", 1.4, "slot")
    cmd.color("sc_annot", "slot")
    cmd.label("slot", '"no anchor"')
    cmd.set("label_position", (2.5, 2.0, 4.0), "slot")
    cmd.hide("cartoon", near("mtr4"))
    cmd.set_view(view)
    render("struct_mtr4_none")
    a = cmd.get_coords("mtr4 and chain A and resi 129 and name OH")
    b = cmd.get_coords("mtr4 and resn ADP and byres (mtr4 within 6 of (mtr4 and chain A and resi 154)) and name N6")
    y129 = min(float(((x - y) ** 2).sum() ** 0.5) for x in a for y in b)  # ligand may carry alt confs
    print(f"Mtr4 Gln-N6 {d['gln']:.1f} A; Y129 OH -> N6 {y129:.1f} A")


def fig2():
    setup()
    load("closed", "3fht")
    load("open", "6b4k")
    cmd.super("open and chain A and resi 102-285 and name CA", "closed and chain A and resi 102-285 and name CA")
    for obj in ("closed", "open"):
        cmd.hide("everything", obj)
        cmd.show("cartoon", f"{obj} and chain A and resi 54-290")
        cmd.color("sc_grey", f"{obj} and chain A")
    cmd.hide("cartoon", "open and chain A and resi 95-290")  # one copy of the shared core
    cmd.set("cartoon_transparency", 0.0, "closed and chain A and resi 76-93")
    cmd.set("cartoon_transparency", 0.0, "open and chain A and resi 54-93")
    cmd.color("sc_sand", "closed and chain A and resi 76-93")
    cmd.color("sc_moss", "open and chain A and resi 54-93")
    sticks("closed and chain A and resi 94", "sc_forest")
    sticks("open and chain A and resi 94", "sc_forest")
    nuc = "closed and resn ANP and byres (closed within 6 of (closed and chain A and resi 119))"
    cmd.show("sticks", nuc)
    cmd.color("sc_terracotta", f"{nuc} and elem C")
    label("closed and chain A and resi 94", "F94 (anchor)", (3.0, -2.5, 5.0))
    label("closed and chain A and resi 80", "closed, RNA-bound", (-7.0, 4.0, 5.0))
    label("open and chain A and resi 62", "open", (3.0, 0.0, 5.0))
    focus = "(closed or open) and chain A and resi 54-100"
    cmd.orient(focus)
    cmd.zoom(focus, 4)
    cmd.hide("cartoon", "closed and chain A and not (closed within 18 of (closed and chain A and resi 94))")
    render("struct_ddx19_pivot")

    setup()
    load("ejc", "2j0q")
    cmd.hide("everything", "ejc")
    cmd.show("cartoon", "ejc and chain A+G")
    cmd.color("sc_grey", "ejc and chain A")
    cmd.color("sc_sand", "ejc and chain G")
    sticks("ejc and chain A and resi 40", "sc_forest")
    sticks("ejc and chain A and resi 41", "sc_terracotta")
    sticks("ejc and chain G and resi 68+73+122+123+124", "sc_sand")
    cmd.color("0xB08A3E", "ejc and chain G and resi 68+73+122+123+124 and elem C")  # darker sand for contrast
    nuc = "ejc and resn ANP and byres (ejc within 6 of (ejc and chain A and resi 65))"
    cmd.show("sticks", nuc)
    cmd.color("sc_terracotta", f"{nuc} and elem C")
    label("ejc and chain A and resi 40", "F40 (anchor)", (-3.0, 1.5, 5.0))
    label("ejc and chain A and resi 41", "D41 (UA+1)", (-1.0, -3.5, 5.0))
    label("ejc and chain G and resi 68", "Y14 R68", (2.5, 1.5, 5.0))
    label("ejc and chain G and resi 73", "Y14 W73", (4.5, -1.0, 5.0))
    focus = "ejc and ((chain A and resi 38-43) or (chain G and resi 68+73+122-124))"
    cmd.orient(focus)
    cmd.zoom(focus, 8)
    cmd.hide("cartoon", "ejc and not (ejc within 18 of (" + focus + "))")
    render("struct_ejc_y14")


if __name__ == "__main__":
    fig1()
    fig2()
