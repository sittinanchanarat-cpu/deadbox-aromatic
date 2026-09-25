#!/usr/bin/env python
"""Stage 0, step 2b: the upstream aromatic anchor (UA; Prp5 W257).

The UA is not an adenine stacker: it sits on a flexible loop ~25 residues
upstream of the Q (~48 upstream of the Walker A Lys) and inserts behind the
Q-motif aromatic, touching the adenine edge (Prp5 W257 CZ2-N6 3.6 A).

For every nucleotide site from sites.tsv with a Walker K, scan aromatics in the
window WK-60..WK-38 (register-free: works with or without a Q motif) and call
  inserted  = ring centroid <= 9.5 A from the purine centroid AND
              closest atom <= 5.0 A from the nucleotide
  displaced = an aromatic is in the window but none meets that
  absent    = no aromatic in a mostly modelled window
  not_modelled = half or more of the window unmodelled (no evidence either way)
Also records unmodelled residues in the loop and the anchor's B-factor z-score,
because the loop's mobility is the point.

Output
  results/ua.tsv   one row per nucleotide site
"""
import sys
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA, RES = ROOT / "data", ROOT / "results"

PURINE = ["N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9"]
RINGS = {
    "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TRP": ["CG", "CD1", "NE1", "CE2", "CD2", "CE3", "CZ2", "CZ3", "CH2"],
    "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"],
}
WINDOW = (-60, -38)          # relative to Walker K label_seq
RING_MAX, CONTACT_MAX = 9.5, 5.0


def xyz(res, names):
    try:
        return np.array([res[n][0].pos.tolist() for n in names])
    except (RuntimeError, IndexError):
        return None


def one(name):
    return gemmi.find_tabulated_residue(name).one_letter_code.upper()


def analyse(pdb, rows):
    st = gemmi.read_structure(str(DATA / "cif" / f"{pdb.lower()}.cif.gz"))
    st.setup_entities()
    model = st[0]
    out = []
    for r in rows.itertuples():
        # Prefer the Q as register anchor: it is ~25 residues from the UA, whereas
        # the Q->Walker K spacer can carry insertions (yeast Prp28: 75 instead of
        # 23). Positions are still reported on the Walker K scale (Q = WK-23).
        # only a genuine Q-motif Gln can serve: right spacer to Walker K, or the
        # Q-motif stacker 5-8 residues upstream of it (covers spacer insertions)
        q_motif = isinstance(r.q_gln, str) and bool(r.q_gln) and (
            (r.q_to_walkerK == r.q_to_walkerK and 21 <= r.q_to_walkerK <= 25)
            or (r.stacker_to_q == r.stacker_to_q and 5 <= r.stacker_to_q <= 8)
            or not (isinstance(r.walker_lys, str) and r.walker_lys))
        if q_motif:
            chain_name, anc = r.q_gln.split(":")
            offset = 23
        else:
            chain_name, anc = r.walker_lys.split(":")
            offset = 0
        chain = model[chain_name]
        wk = next((x for x in chain if x.seqid.num == int(anc[3:])), None)
        # ligand residue: site id is PDB:chain:NAMEnum
        _, lch, lid = r.site.split(":")
        lig = next((x for x in model[lch] if x.name == lid[:3] and str(x.seqid.num) == lid[3:]), None)
        if wk is None or lig is None or wk.label_seq is None:
            continue
        pur = xyz(lig, PURINE)
        if pur is None:
            continue
        pc = pur.mean(0)
        by_ls = {x.label_seq: x for x in chain if x.label_seq}
        protein_b = np.array([np.mean([a.b_iso for a in x]) for x in chain if x.het_flag != "H"])
        k_ls = wk.label_seq + offset       # (virtual) Walker K position
        lo, hi = k_ls + WINDOW[0], k_ls + WINDOW[1]
        missing = sum(1 for i in range(lo, hi + 1) if i not in by_ls)
        cands = []
        for i in range(lo, hi + 1):
            x = by_ls.get(i)
            if x is None or x.name not in RINGS:
                continue
            ring = xyz(x, RINGS[x.name])
            if ring is None:
                continue
            d_ring = float(np.linalg.norm(ring.mean(0) - pc))
            d_min = min(a.pos.dist(b.pos) for a in x for b in lig)
            bz = (np.mean([a.b_iso for a in x]) - protein_b.mean()) / protein_b.std()
            cands.append((x, d_ring, d_min, bz))
        ins = [c for c in cands if c[1] <= RING_MAX and c[2] <= CONTACT_MAX]
        best = min(ins, key=lambda c: c[1]) if ins else (min(cands, key=lambda c: c[1]) if cands else None)
        width = hi - lo + 1
        if ins:
            state = "inserted"
        elif cands:
            state = "displaced"
        elif missing >= 0.5 * width:
            state = "not_modelled"   # construct starts later or loop disordered: no evidence either way
        else:
            state = "absent"
        # the segment itself, for eyeballing register and for later logos
        seg = "".join(one(by_ls[i].name) if i in by_ls else "-" for i in range(lo, hi + 1))
        out.append({
            "site": r.site, "pdb": pdb, "ua_state": state,
            "ua": f"{best[0].name}{best[0].seqid.num}" if best else "",
            "ua_aa": one(best[0].name) if best else "",
            "ua_to_walkerK": (k_ls - best[0].label_seq) if best else None,
            "anchor": "walkerK" if offset == 0 else "Q",
            "ua_ring_to_adenine": round(best[1], 2) if best else None,
            "ua_min_to_nucleotide": round(best[2], 2) if best else None,
            "ua_bfactor_z": round(best[3], 2) if best else None,
            "loop_missing": missing,
            "n_window_aromatics": len(cands),
            "window_seq": seg,
        })
    return out


def main():
    sites = pd.read_csv(RES / "sites.tsv", sep="\t")
    # need a register anchor: Walker K, or failing that the Q-motif Gln
    has_k = sites.walker_lys.notna() & (sites.walker_lys != "")
    has_q = sites.q_gln.notna() & (sites.q_gln != "")
    sites = sites[has_k | has_q]
    rows = []
    for pdb, g in sites.groupby("pdb"):
        try:
            rows += analyse(pdb, g)
        except Exception as ex:
            print(f"{pdb}: {ex}", file=sys.stderr)
    df = pd.DataFrame(rows)
    df.to_csv(RES / "ua.tsv", sep="\t", index=False)
    print(df.ua_state.value_counts().to_string(), file=sys.stderr)


if __name__ == "__main__":
    main()
