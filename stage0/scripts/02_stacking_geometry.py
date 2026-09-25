#!/usr/bin/env python
"""Stage 0, step 2: define the adenine-stacking residue geometrically.

For every adenine nucleotide bound to an SF2 helicase chain, record:
  * aromatic side chains (F/W/Y/H) near the purine, with ring-centroid distance
    and interplanar angle; stacked = dist <= 5.0 A and angle <= 30 deg
  * Gln side chains H-bonding adenine N6/N7 (the Q-motif glutamine)
  * the motif I (Walker A) Lys whose NZ contacts the phosphates -- a
    register anchor that exists whether or not a Q motif does

Sequence offsets use label_seq_id, i.e. positions in the deposited entity
sequence, so unmodelled gaps do not distort spacer lengths.

Outputs
  results/contacts.tsv  one row per (nucleotide site, nearby aromatic)
  results/sites.tsv     one row per nucleotide site, with the call
"""
import sys
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA, RES = ROOT / "data", ROOT / "results"

ADENINE_LIGS = {"ATP", "ADP", "ANP", "ACP", "AGS", "APC", "AMP", "ADX", "AN2"}
PURINE = ["N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9"]
PHOSPHATE_ATOMS = {"PA", "PB", "PG", "O1A", "O2A", "O3A", "O1B", "O2B", "O3B",
                   "O1G", "O2G", "O3G", "N3B", "C3B", "S1G"}
RINGS = {
    "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TRP": ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"],  # benzo ring; indole checked too
    "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"],
}
TRP_INDOLE = ["CG", "CD1", "NE1", "CE2", "CD2", "CE3", "CZ2", "CZ3", "CH2"]

STACK_DIST, STACK_ANGLE = 5.0, 30.0
SCAN_DIST = 7.0      # report aromatics out to here for manual inspection
QHB_DIST = 3.5       # Gln OE1/NE2 to adenine N6/N7
WALKER_DIST = 4.5    # Lys NZ to phosphate atoms


def plane(coords):
    c = coords.mean(axis=0)
    _, _, vt = np.linalg.svd(coords - c)
    return c, vt[2]


def angle(n1, n2):
    a = np.degrees(np.arccos(np.clip(abs(np.dot(n1, n2)), 0, 1)))
    return min(a, 180 - a)


def atoms_xyz(res, names):
    out = []
    for n in names:
        a = res.find_atom(n, "*")
        if a is None:
            return None
        out.append(a.pos.tolist())
    return np.array(out)


def residue_label_seq(res):
    return res.label_seq if res.label_seq is not None else None


def analyse(path, sf2_chains):
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    model = st[0]
    ns = gemmi.NeighborSearch(model, st.cell, 8).populate()
    sites, contacts = [], []
    for chain in model:
        for lig in chain:
            if lig.name not in ADENINE_LIGS:
                continue
            pur = atoms_xyz(lig, PURINE)
            if pur is None:
                continue
            pc, pn = plane(pur)
            center = gemmi.Position(*pc)
            site_id = f"{path.name[:4].upper()}:{chain.name}:{lig.name}{lig.seqid.num}"

            # aromatic candidates
            seen = set()
            cands = []
            for m in ns.find_atoms(center, "\0", radius=SCAN_DIST + 3):
                cra = m.to_cra(model)
                r = cra.residue
                if r.name not in RINGS or cra.chain.name not in sf2_chains:
                    continue
                key = (cra.chain.name, r.seqid.num, r.seqid.icode)
                if key in seen:
                    continue
                seen.add(key)
                ring = atoms_xyz(r, RINGS[r.name])
                if ring is None:
                    continue
                rc, rn = plane(ring)
                if r.name == "TRP":
                    ind = atoms_xyz(r, TRP_INDOLE)
                    if ind is not None:
                        ic, inn = plane(ind)
                        if np.linalg.norm(ic - pc) < np.linalg.norm(rc - pc):
                            rc, rn = ic, inn
                d = float(np.linalg.norm(rc - pc))
                if d > SCAN_DIST:
                    continue
                ang = angle(pn, rn)
                offset = float(np.linalg.norm((rc - pc) - np.dot(rc - pc, pn) * pn))
                cands.append({
                    "site": site_id, "chain": cra.chain.name, "resname": r.name,
                    "auth_num": r.seqid.num, "label_seq": residue_label_seq(r),
                    "centroid_dist": round(d, 2), "plane_angle": round(ang, 1),
                    "lateral_offset": round(offset, 2),
                    "stacked": d <= STACK_DIST and ang <= STACK_ANGLE,
                })

            # Q-motif glutamine: side chain H-bond to N6 or N7
            qhits = []
            for aname in ("N6", "N7"):
                a = lig.find_atom(aname, "*")
                if a is None:
                    continue
                for m in ns.find_atoms(a.pos, "\0", radius=QHB_DIST):
                    cra = m.to_cra(model)
                    if (cra.residue.name == "GLN" and cra.atom.name in ("OE1", "NE2")
                            and cra.chain.name in sf2_chains):
                        qhits.append((cra.chain.name, cra.residue.seqid.num,
                                      residue_label_seq(cra.residue), aname,
                                      round(cra.atom.pos.dist(a.pos), 2)))

            # Walker A lysine
            khits = []
            for a in lig:
                if a.name not in PHOSPHATE_ATOMS:
                    continue
                for m in ns.find_atoms(a.pos, "\0", radius=WALKER_DIST):
                    cra = m.to_cra(model)
                    if (cra.residue.name == "LYS" and cra.atom.name == "NZ"
                            and cra.chain.name in sf2_chains):
                        khits.append((cra.chain.name, cra.residue.seqid.num,
                                      residue_label_seq(cra.residue),
                                      round(cra.atom.pos.dist(a.pos), 2)))
            # keep closest contact per Lys
            kbest = {}
            for c, num, ls, d in khits:
                if (c, num) not in kbest or d < kbest[(c, num)][2]:
                    kbest[(c, num)] = (ls, num, d)

            if not cands and not qhits and not kbest:
                continue  # nucleotide not bound to an SF2 chain

            stacked = sorted((c for c in cands if c["stacked"]), key=lambda c: c["centroid_dist"])
            # prefer a Lys inside a GK[ST] motif I; a Lys elsewhere can also touch
            # the phosphates (yeast Prp28 4W7S picked one 75 residues from the Q)
            def in_motif_I(item):
                (c, num), _ = item
                ch = model[c]
                idx = next((i for i, x in enumerate(ch) if x.seqid.num == num), None)
                if idx is None or idx == 0 or idx + 1 >= len(ch):
                    return False
                return ch[idx - 1].name == "GLY" and ch[idx + 1].name in ("SER", "THR")
            motif_I = [kv for kv in kbest.items() if in_motif_I(kv)]
            pool = motif_I or list(kbest.items())
            walker = min(pool, key=lambda kv: kv[1][2]) if pool else None
            q = min(qhits, key=lambda h: h[4]) if qhits else None
            # RecA1 stacker = same chain as Walker K and upstream of it; the Q-motif
            # aromatic lives here. Anything else (RecA2 Phe in DEAH, other chain) is
            # reported separately so it can't masquerade as the RecA1 residue.
            wk = walker[1][0] if walker else None
            wchain = walker[0][0] if walker else None

            def domain(c):
                if c["label_seq"] is None:
                    return "unassigned"
                if wk is not None and c["chain"] == wchain:
                    return "RecA1" if c["label_seq"] < wk else "RecA2"
                # no Walker K contact (e.g. Mg-free ADP): fall back to the Q anchor
                if q and c["chain"] == q[0] and q[2] is not None:
                    return "RecA1" if 0 < q[2] - c["label_seq"] <= 40 else "unassigned"
                return "unassigned"

            for c in cands:
                c["domain"] = domain(c)
            recA1 = [c for c in stacked if c["domain"] == "RecA1"]
            other = [c for c in stacked if c["domain"] != "RecA1"]
            best = recA1[0] if recA1 else (stacked[0] if stacked and wk is None else None)

            def diff(a, b):
                return (a - b) if a is not None and b is not None else None

            sites.append({
                "site": site_id, "pdb": path.name[:4].upper(), "ligand": lig.name,
                "n_stacked": len(stacked), "has_recA1_stacker": bool(recA1),
                "stacker": f"{best['resname']}{best['auth_num']}" if best else "",
                "stacker_chain": best["chain"] if best else "",
                "stacker_dist": best["centroid_dist"] if best else None,
                "stacker_angle": best["plane_angle"] if best else None,
                "all_stackers": ";".join(f"{c['chain']}:{c['resname']}{c['auth_num']}({c['domain']})" for c in stacked),
                "other_stacker": ";".join(f"{c['chain']}:{c['resname']}{c['auth_num']}({c['domain']})" for c in other),
                "q_gln": f"{q[0]}:GLN{q[1]}" if q else "",
                "q_contacts": ";".join(f"{h[3]}={h[4]}" for h in qhits if q and h[1] == q[1]),
                "walker_lys": f"{walker[0][0]}:LYS{walker[1][1]}" if walker else "",
                # spacer: stacker to Q (paper's Q-n register) and stacker to Walker K
                "stacker_to_q": diff(q[2] if q else None, best["label_seq"] if best else None)
                                if best and q and best["chain"] == q[0] else None,
                "stacker_to_walkerK": diff(walker[1][0] if walker else None, best["label_seq"] if best else None)
                                      if best and walker and best["chain"] == walker[0][0] else None,
                "q_to_walkerK": diff(walker[1][0] if walker else None, q[2] if q else None)
                                if q and walker and q[0] == walker[0][0] else None,
            })
            contacts.extend(cands)
    return sites, contacts


def main():
    ent = pd.read_csv(DATA / "entries.tsv", sep="\t", dtype=str)
    chains = ent.groupby("pdb")["auth_chains"].apply(
        lambda s: {c for v in s.dropna() for c in v.split(",")}).to_dict()
    all_sites, all_contacts = [], []
    for pdb, ch in sorted(chains.items()):
        path = DATA / "cif" / f"{pdb.lower()}.cif.gz"
        if not path.exists():
            continue
        try:
            s, c = analyse(path, ch)
        except Exception as ex:
            print(f"{pdb}: {ex}", file=sys.stderr)
            continue
        all_sites += s
        all_contacts += c
    RES.mkdir(exist_ok=True)
    pd.DataFrame(all_sites).to_csv(RES / "sites.tsv", sep="\t", index=False)
    pd.DataFrame(all_contacts).to_csv(RES / "contacts.tsv", sep="\t", index=False)
    print(f"{len(all_sites)} nucleotide sites on SF2 chains", file=sys.stderr)


if __name__ == "__main__":
    main()
