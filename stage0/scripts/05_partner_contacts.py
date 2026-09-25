#!/usr/bin/env python
"""Partner-complex census, step 2: who touches the upstream anchor (UA)?

For every DEAD-box chain in data/partner_entities.tsv:
  * UA state, judged without needing a nucleotide:
      inserted     UA ring centroid <= 8.0 A from the Q-motif aromatic ring
                   (Stage 0 inserted anchors sit 6.3-7.1 A behind it)
      displaced    UA modelled but further away
      not_modelled UA residue absent from the model
      no_UA        no aromatic at Q-21..Q-30 in the sequence
  * other chains within 4.5 A of the UA residue, and of its loop (UA-8..UA+3)
    - "partner"  = a different entity (protein or nucleic acid)
    - "self"     = another copy of the same entity (often crystal packing)
    Only the deposited coordinates count: crystal-symmetry images are ignored,
    and so is any partner chain that is not the cognate copy of its entity (the
    copy with the largest interface to this helicase chain), because in
    multi-complex asymmetric units the other copies touch it through packing.
    (v1 counted symmetry images; that made lattice contacts in DDX19B-Gle1 and
    DDX19B-Nup214 crystals look like loop binding. Kept as
    results/partner_contacts.v1_with_symmetry.tsv.)

Residue positions come from the entity sequence (label_seq), so they match
the sequence-based UA/Q calls in 04_partner_fetch.py.

Output
  results/partner_contacts.tsv   one row per DEAD-box chain
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
ADENINE = {"ATP", "ADP", "ANP", "ACP", "AGS", "APC", "AMP", "ADX", "AN2"}
CONTACT, INSERTED = 4.5, 8.0
LOOP = (-8, 3)


def ring_centre(res):
    try:
        return np.array([res[a][0].pos.tolist() for a in RINGS[res.name]]).mean(0)
    except (KeyError, RuntimeError, IndexError):
        return None


def cognate_chains(model, ns, own_chain, entity_of):
    """For each other entity, the chain with most atom contacts to own_chain."""
    size = {}
    for res in model[own_chain]:
        for atom in res:
            for mk in ns.find_atoms(atom.pos, "\0", radius=CONTACT):
                cra = mk.to_cra(model)
                # image_idx is 0 for pure lattice translations too, so test the
                # deposited coordinates directly to drop every symmetry image
                if cra.atom.pos.dist(atom.pos) > CONTACT:
                    continue
                if cra.chain.name == own_chain or cra.chain.name not in entity_of:
                    continue
                size[cra.chain.name] = size.get(cra.chain.name, 0) + 1
    best = {}
    for c, n in size.items():
        ent = entity_of[c][0]
        if ent not in best or n > size[best[ent]]:
            best[ent] = c
    return set(best.values())


def contacts(model, ns, residues, own_chain, entity_of, allowed):
    """{(kind, entity_description): n_atom_contacts} for atoms of `residues`."""
    out = {}
    for res in residues:
        for atom in res:
            for mk in ns.find_atoms(atom.pos, "\0", radius=CONTACT):
                cra = mk.to_cra(model)
                # image_idx is 0 for pure lattice translations too, so test the
                # deposited coordinates directly to drop every symmetry image
                if cra.atom.pos.dist(atom.pos) > CONTACT:
                    continue
                if cra.chain.name == own_chain or cra.residue.het_flag == "H":
                    continue
                if cra.chain.name not in allowed:
                    continue
                ent = entity_of.get(cra.chain.name)
                if ent is None:
                    continue
                kind = "self" if ent[0] == entity_of[own_chain][0] else "partner"
                key = (kind, ent[1])
                out[key] = out.get(key, 0) + 1
    return out


def analyse(pdb, ents):
    st = gemmi.read_structure(str(DATA / "partner_cif" / f"{pdb.lower()}.cif.gz"))
    st.setup_entities()
    model = st[0]
    ns = gemmi.NeighborSearch(model, st.cell, 6).populate()
    entity_of = {}
    for e in ents.itertuples():
        for c in str(e.chains).split(","):
            entity_of[c] = (e.entity, f"{e.type}:{e.description}")
    rows = []
    for e in ents[ents.is_dead].itertuples():
        for cname in str(e.chains).split(","):
            chain = model.find_chain(cname)
            if chain is None:
                continue
            by = {r.label_seq: r for r in chain if r.label_seq}
            row = {"pdb": pdb, "chain": cname, "entity": e.entity, "protein": e.description,
                   "uniprot": e.uniprot, "organism": e.organism, "method": e.method,
                   "resolution": e.resolution, "ua_aa": e.ua_aa if isinstance(e.ua_aa, str) else "",
                   "ua_seq": e.ua_seq, "q_seq": e.q_seq}
            if not isinstance(e.ua_aa, str) or not e.ua_aa:
                row["ua_state"] = "no_UA"
                rows.append(row)
                continue
            ua_ls, q_ls = int(e.ua_seq), int(e.q_seq)
            ua = by.get(ua_ls)
            row["loop_modelled"] = sum(1 for i in range(ua_ls + LOOP[0], ua_ls + LOOP[1] + 1) if i in by)
            if ua is None or ua.name not in RINGS:
                row["ua_state"] = "not_modelled"
                rows.append(row)
                continue
            row["ua_auth"] = f"{ua.name}{ua.seqid.num}"
            # backbone-only measure for models without side chains (low-res cryo-EM)
            qres = by.get(q_ls)
            if qres is not None and ua.find_atom("CA", "*") and qres.find_atom("CA", "*"):
                row["ua_ca_to_q_ca"] = round(ua.find_atom("CA", "*").pos.dist(qres.find_atom("CA", "*").pos), 2)
            row["ua_sidechain"] = ring_centre(ua) is not None
            # Q-motif aromatic: Q-7 or Q-6
            qa = next((by[q_ls - o] for o in (7, 6) if q_ls - o in by and by[q_ls - o].name in RINGS), None)
            uc = ring_centre(ua)
            if qa is not None and uc is not None and ring_centre(qa) is not None:
                d = float(np.linalg.norm(uc - ring_centre(qa)))
                row["ua_to_qaromatic"] = round(d, 2)
                # indole centroid sits further back than a phenyl centroid
                cut = INSERTED + (0.5 if ua.name == "TRP" else 0.0)
                row["ua_state"] = "inserted" if d <= cut else "displaced"
            elif "ua_ca_to_q_ca" in row:
                # side chains not built: backbone rule calibrated on the ring-based
                # calls (inserted 9.2 +/- 0.5 A, 8.6-11.8; displaced 9.5-11.8)
                # Trp sits further out when inserted (9.9-10.4 A, n=6)
                bb = row["ua_ca_to_q_ca"]
                lo, hi = (10.6, 12.0) if ua.name == "TRP" else (10.0, 12.0)
                row["ua_state"] = ("inserted(bb)" if bb <= lo else
                                   "displaced(bb)" if bb >= hi else "ambiguous(bb)")
            else:
                row["ua_state"] = "q_missing"
            # nucleotide in this chain's pocket, if any
            nuc = [r for ch in model for r in ch if r.name in ADENINE and r.find_atom("N6", "*")]
            if nuc and uc is not None:
                n6 = min((r.find_atom("N6", "*").pos for r in nuc), key=lambda p: np.linalg.norm(np.array(p.tolist()) - uc))
                row["ua_ring_to_N6"] = round(float(np.linalg.norm(np.array(n6.tolist()) - uc)), 2)
            loop = [by[i] for i in range(ua_ls + LOOP[0], ua_ls + LOOP[1] + 1) if i in by]
            allowed = cognate_chains(model, ns, cname, entity_of)
            c_ua = contacts(model, ns, [ua], cname, entity_of, allowed)
            c_loop = contacts(model, ns, loop, cname, entity_of, allowed)
            fmt = lambda c, kind: "; ".join(f"{k[1][:50]}({n})" for k, n in sorted(c.items(), key=lambda x: -x[1]) if k[0] == kind)
            row.update({
                "ua_partner": fmt(c_ua, "partner"), "ua_self": fmt(c_ua, "self"),
                "loop_partner": fmt(c_loop, "partner"), "loop_self": fmt(c_loop, "self"),
                "n_ua_partner_atoms": sum(n for k, n in c_ua.items() if k[0] == "partner"),
                "n_loop_partner_atoms": sum(n for k, n in c_loop.items() if k[0] == "partner"),
            })
            rows.append(row)
    return rows


def main():
    ents = pd.read_csv(DATA / "partner_entities.tsv", sep="\t")
    rows = []
    for pdb, g in ents.groupby("pdb"):
        if not (DATA / "partner_cif" / f"{pdb.lower()}.cif.gz").exists():
            continue
        try:
            rows += analyse(pdb, g)
        except Exception as ex:
            print(f"{pdb}: {ex}", file=sys.stderr)
    df = pd.DataFrame(rows)
    df.to_csv(RES / "partner_contacts.tsv", sep="\t", index=False)
    print(df.ua_state.value_counts().to_string(), file=sys.stderr)


if __name__ == "__main__":
    main()
