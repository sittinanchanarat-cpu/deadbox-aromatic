#!/usr/bin/env python
"""Characterise the upstream anchor (UA) -- adenine N6 contact.

Tanner et al. (2003) saw the upstream Phe touch the Q-motif Pro; Stage 0 found
its ring edge also 3.3-4.1 A from adenine N6. For every DEAD-box site with an
inserted UA (results/ua.tsv, register WK-44..52), measure:
  n6_to_ring       N6 -> ring centroid distance
  approach_angle   angle between the ring normal and the centroid->N6 vector
                   (~0 deg = N6 over the ring face, NH-pi; ~90 deg = in-plane edge)
  n6_to_ring_atom  closest ring atom and distance
  q_hbond_n6       Q-motif Gln OE1/NE2 -> N6 distance (same N6 read twice?)
plus nucleotide identity (ATP analogue / ADP-state / AMP) and whether RNA is
present in the entry (proxy for the closed, RNA-bound state).

Output
  results/n6_geometry.tsv
"""
import sys
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA, RES = ROOT / "data", ROOT / "results"
RINGS = {
    "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TRP": ["CG", "CD1", "NE1", "CE2", "CD2", "CE3", "CZ2", "CZ3", "CH2"],
    "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"],
}
STATE = {"ATP": "ATP-like", "ANP": "ATP-like", "ACP": "ATP-like", "AGS": "ATP-like",
         "APC": "ATP-like", "ADP": "ADP", "AMP": "AMP", "ADX": "other", "AN2": "other"}
TS_MIMICS = {"ALF", "AF3", "BEF", "VO4", "MGF"}


def plane(xyz):
    c = xyz.mean(0)
    return c, np.linalg.svd(xyz - c)[2][2]


def main():
    ua = pd.read_csv(RES / "ua.tsv", sep="\t")
    sites = pd.read_csv(RES / "sites_labelled.tsv", sep="\t")
    ua = ua.merge(sites[["site", "family", "protein", "q_gln"]], on="site")
    ua = ua[(ua.ua_state == "inserted") & ua.ua_to_walkerK.between(44, 52)]
    ent = pd.read_csv(DATA / "entries.tsv", sep="\t", dtype=str).fillna("")
    ligs = ent.groupby("pdb").ligands.first().to_dict()
    rows = []
    for pdb, g in ua.groupby("pdb"):
        st = gemmi.read_structure(str(DATA / "cif" / f"{pdb.lower()}.cif.gz"))
        st.setup_entities()
        model = st[0]
        has_rna = any(e.polymer_type == gemmi.PolymerType.Rna for e in st.entities)
        for r in g.itertuples():
            _, lch, lid = r.site.split(":")
            lig = next((x for x in model[lch] if x.name == lid[:3] and str(x.seqid.num) == lid[3:]), None)
            if lig is None or lig.find_atom("N6", "*") is None:
                continue
            n6 = np.array(lig.find_atom("N6", "*").pos.tolist())
            # the UA residue: the ring nearest N6 with this name/number
            num, name = int(r.ua[3:]), r.ua[:3]
            cands = [x for ch in model for x in ch if x.seqid.num == num and x.name == name]
            best = None
            for x in cands:
                try:
                    xyz = np.array([x[a][0].pos.tolist() for a in RINGS[name]])
                except (RuntimeError, IndexError):
                    continue
                c, nrm = plane(xyz)
                d = float(np.linalg.norm(n6 - c))
                if best is None or d < best[0]:
                    best = (d, c, nrm, x, xyz)
            if best is None:
                continue
            d, c, nrm, x, xyz = best
            v = (n6 - c) / d
            ang = float(np.degrees(np.arccos(abs(np.dot(v, nrm)))))
            names = RINGS[name]
            j = int(np.argmin(np.linalg.norm(xyz - n6, axis=1)))
            qd = None
            if isinstance(r.q_gln, str) and r.q_gln:
                qc, qn = r.q_gln.split(":")
                qres = next((y for y in model[qc] if y.seqid.num == int(qn[3:])), None)
                if qres is not None:
                    qd = min((a.pos.dist(lig.find_atom("N6", "*").pos) for a in qres if a.name in ("OE1", "NE2")), default=None)
            entry_ligs = set(ligs.get(pdb, "").split(","))
            rows.append({
                "site": r.site, "pdb": pdb, "protein": r.protein, "ua": r.ua,
                "nucleotide": lig.name, "state": STATE.get(lig.name, "other") + ("+TS-mimic" if entry_ligs & TS_MIMICS and lig.name == "ADP" else ""),
                "rna_present": has_rna,
                "n6_to_ring": round(d, 2), "approach_angle": round(ang, 1),
                "n6_to_ring_atom": f"{names[j]} {np.linalg.norm(xyz[j] - n6):.2f}",
                "min_ring_atom_dist": round(float(np.linalg.norm(xyz[j] - n6)), 2),
                "q_hbond_n6": round(qd, 2) if qd is not None else None,
            })
    df = pd.DataFrame(rows)
    df.to_csv(RES / "n6_geometry.tsv", sep="\t", index=False)
    print(f"{len(df)} inserted-UA sites", file=sys.stderr)


if __name__ == "__main__":
    main()
