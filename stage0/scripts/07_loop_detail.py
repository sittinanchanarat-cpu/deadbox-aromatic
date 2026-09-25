#!/usr/bin/env python
"""Anchor-loop detail: which loop residues do partners touch, and does the loop move?

For each helicase in PANELS:
  * contacts   every residue in the window UA-12..UA+6: partner residues within
               4.5 A in the deposited coordinates (no symmetry images; only the
               cognate copy of each partner entity), and the SASA it loses when
               the partner binds (Shrake-Rupley, chain alone vs with partners)
  * conformation  each chain is superposed on a free reference chain using RecA1
               core CA atoms (UA+8..UA+190, outliers >2 A dropped iteratively);
               per window residue CA shift, UA ring-centroid shift, UA chi1/chi2.
               Each chain is compared with its nearest free chain from another
               entry (lowest loop RMSD), so a nucleotide-state difference is not
               read as a partner effect; free-vs-free gives the baseline.

Numbering follows the reference chain (the first "free" entry of each panel);
every other chain is mapped onto it by pairwise sequence alignment, because old
entries use a different register (e.g. yeast eIF4A 1FUU F23 = 2VSX F24).

Output
  results/loop_contacts.tsv      one row per (chain, window residue)
  results/loop_conformation.tsv  one row per (chain, window residue)
  results/loop_summary.tsv       one row per chain
"""
import io
import urllib.request
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser
from Bio.PDB.SASA import ShrakeRupley

ROOT = Path(__file__).resolve().parents[1]
CIF, RES = ROOT / "data" / "loop_cif", ROOT / "results"
WIN = (-12, 6)
CORE = (8, 190)
CONTACT = 4.5

# protein: (UA number in the reference chain, [(pdb, helicase chain, label)])
# The reference is the first "free" chain (no protein partner at the loop;
# nucleotide/RNA allowed). Partners are found automatically: the cognate chain
# of every other entity (see cognate_chains), crystal-symmetry images excluded.
PANELS = {
    "DDX19B (human)": (94, [
        ("3FHT", "A", "free"), ("3FHT", "B", "free"), ("3G0H", "A", "free"),
        ("3EWS", "A", "free"), ("3EWS", "B", "free"), ("6B4K", "A", "free"), ("6B4K", "B", "free"),
        ("6B4I", "E", "Gle1"), ("6B4I", "F", "Gle1"), ("6B4J", "E", "Gle1"), ("6B4J", "F", "Gle1"),
        ("3FHC", "B", "Nup214"), ("3FMO", "B", "Nup214"), ("3FMP", "B", "Nup214"), ("3FMP", "D", "Nup214"),
    ]),
    "Dbp5 (yeast)": (94, [
        ("3PEW", "A", "free"), ("3PEY", "A", "free"), ("5ELX", "A", "free"),
        ("3RRM", "A", "Gle1+Nup159"), ("3RRN", "A", "Gle1"),
    ]),
    "eIF4A (yeast)": (23, [
        ("1FUU", "A", "free"), ("1FUU", "B", "free"), ("1QDE", "A", "free"), ("1QVA", "A", "free"),
        ("2VSO", "A", "eIF4G"), ("2VSO", "B", "eIF4G"), ("2VSX", "A", "eIF4G"), ("2VSX", "B", "eIF4G"),
    ]),
    "eIF4A1 (human)": (34, [
        ("2G9N", "A", "free"), ("2G9N", "B", "free"), ("5ZC9", "A", "free"),
        ("9DTS", "A", "free"), ("9DTS", "B", "free"), ("9DTS", "C", "free"), ("9DTS", "D", "free"),
        ("3EIQ", "A", "PDCD4"), ("2ZU6", "C", "PDCD4"), ("8HUJ", "A", "eIF4G"), ("6ZMW", "j", "eIF4G(48S)"),
    ]),
    "eIF4AIII (human)": (40, [
        ("2HXY", "A", "free"), ("2HXY", "B", "free"), ("2HXY", "C", "free"), ("2HXY", "D", "free"),
        ("2J0U", "A", "free"), ("2J0U", "B", "free"),
        ("2J0S", "A", "EJC"), ("2J0Q", "A", "EJC"), ("2J0Q", "B", "EJC"), ("2XB2", "A", "EJC+UPF3B"),
        ("2XB2", "X", "EJC+UPF3B"), ("3EX7", "C", "EJC"), ("3EX7", "H", "EJC"),
        ("7ZNJ", "A", "EJC(cryoEM)"), ("7ZNJ", "F", "EJC(cryoEM)"), ("4C9B", "A", "CWC22"),
    ]),
    "DDX6 (human)": (98, [
        ("4CT5", "A", "free"), ("1VEC", "A", "free"),
        ("4CT4", "B", "CNOT1"), ("4CT4", "D", "CNOT1"), ("5ANR", "B", "CNOT1+4E-T"),
    ]),
    # no free structure: contacts only, conformation relative to the first chain
    "Spb4 (yeast)": (7, [("7NAC", "x", "pre-60S"), ("7NAD", "x", "pre-60S")]),
    "DDX54 (human)": (98, [("8FKV", "NI", "pre-60S"), ("8FKW", "NI", "pre-60S")]),
}
# reference numbering -> UniProt (1FUU drops Met1, Tanner 2003 numbering: F23 = UniProt F24)
UNIPROT_OFFSET = {"eIF4A (yeast)": 1}
RING = ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]


def load(pdb, cache={}):
    if pdb not in cache:
        path = CIF / f"{pdb.lower()}.cif.gz"
        if not path.exists():
            CIF.mkdir(parents=True, exist_ok=True)
            urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pdb.upper()}.cif.gz", path)
        st = gemmi.read_structure(str(path))
        st.setup_entities()
        st.remove_hydrogens()
        cache[pdb] = st
    return cache[pdb]


def residues(chain):
    return {r.seqid.num: r for r in chain if r.het_flag == "A" and r.seqid.icode == " "}


def mapped(chain_res, ref_res):
    """Re-key chain residues by the reference numbering via sequence alignment.
    Unmodelled stretches are filled with X so numbering gaps are respected."""
    from Bio import Align
    al = Align.PairwiseAligner(mode="global", open_gap_score=-10, extend_gap_score=-0.5)
    fill = lambda d: list(range(min(d), max(d) + 1))
    ks, kr = fill(chain_res), fill(ref_res)
    one = lambda d, ks: "".join((gemmi.find_tabulated_residue(d[k].name).one_letter_code.upper() or "X") if k in d else "X" for k in ks)
    a = al.align(one(chain_res, ks), one(ref_res, kr))[0]
    out = {}
    for (s0, s1), (r0, r1) in zip(*a.aligned):
        for i, j in zip(range(s0, s1), range(r0, r1)):
            if ks[i] in chain_res and kr[j] in ref_res:
                out[kr[j]] = chain_res[ks[i]]
    return out


def ca(res):
    a = res.find_atom("CA", "*")
    return np.array(a.pos.tolist()) if a else None


def ring(res):
    try:
        return np.array([res[a][0].pos.tolist() for a in RING]).mean(0)
    except (RuntimeError, IndexError):
        return None


def dihedral(res, names):
    try:
        p = [res[n][0].pos for n in names]
    except (RuntimeError, IndexError):
        return np.nan
    return round(float(np.degrees(gemmi.calculate_dihedral(*p))), 0)


def kabsch(P, Q):
    """Rotation R, translation t so that P @ R.T + t ~ Q."""
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - pc @ R.T


def superpose(mob, ref, keys):
    keys = [k for k in keys if k in mob and k in ref and ca(mob[k]) is not None and ca(ref[k]) is not None]
    P = np.array([ca(mob[k]) for k in keys])
    Q = np.array([ca(ref[k]) for k in keys])
    use = np.ones(len(keys), bool)
    for _ in range(5):
        R, t = kabsch(P[use], Q[use])
        d = np.linalg.norm(P @ R.T + t - Q, axis=1)
        new = d < 2.0
        if new.sum() < 30 or (new == use).all():
            break
        use = new
    return R, t, int(use.sum()), float(np.sqrt((d[use] ** 2).mean()))


def sasa_by_res(pdb, chains):
    """Per-residue SASA of the first chain, computed with only `chains` present."""
    st = load(pdb).clone()
    model = st[0]
    for name in [c.name for c in model]:
        if name not in chains:
            model.remove_chain(name)
    model.remove_ligands_and_waters()
    doc = st.make_mmcif_document()
    s = MMCIFParser(QUIET=True).get_structure(pdb, io.StringIO(doc.as_string()))
    ShrakeRupley().compute(s[0], level="R")
    out = {}
    for r in s[0][chains[0]]:
        if r.id[0] == " " and r.id[2] == " ":
            out[r.id[1]] = r.sasa
    return out


def cognate_chains(st, own):
    """{chain: entity description} for the chain of each other polymer entity
    with the most atom contacts to `own`; symmetry images are ignored."""
    model = st[0]
    ns = gemmi.NeighborSearch(model, st.cell, 6).populate()
    ent_of = {}
    for e in st.entities:
        if e.entity_type == gemmi.EntityType.Polymer:
            for sub in e.subchains:
                for ch in model:
                    if ch.get_subchain(sub):
                        ent_of[ch.name] = (e.name, e.polymer_type.name, e.full_sequence and len(e.full_sequence))
    size = {}
    for res in model[own]:
        for atom in res:
            for mk in ns.find_atoms(atom.pos, "\0", radius=CONTACT):
                cra = mk.to_cra(model)
                if cra.chain.name == own or cra.chain.name not in ent_of or cra.residue.het_flag != "A":
                    continue
                if cra.atom.pos.dist(atom.pos) > CONTACT:  # lattice image
                    continue
                size[cra.chain.name] = size.get(cra.chain.name, 0) + 1
    own_ent = ent_of.get(own, (None,))[0]
    best = {}
    for c, n in size.items():
        e = ent_of[c][0]
        if e != own_ent and (e not in best or n > size[best[e]]):
            best[e] = c
    return {c: e for e, c in best.items()}, ns


ENTITY_NAMES = {}


def entity_names(pdb):
    if pdb not in ENTITY_NAMES:
        b = gemmi.cif.read(str(CIF / f"{pdb.lower()}.cif.gz")).sole_block()
        ENTITY_NAMES[pdb] = dict(zip(b.find_values("_entity.id"), b.find_values("_entity.pdbx_description")))
    return ENTITY_NAMES[pdb]


def loop_rmsd(res, other, win, core):
    R, t, n, rms = superpose(res, other, core)
    d = [float(np.linalg.norm(ca(res[k]) @ R.T + t - ca(other[k])))
         for k in win if k in res and k in other and ca(res[k]) is not None and ca(other[k]) is not None]
    return (float(np.sqrt(np.mean(np.square(d)))) if len(d) >= 5 else np.nan), len(d), R, t


def main():
    crows, frows, srows = [], [], []
    for prot, (ua_num, entries) in PANELS.items():
        ref_pdb, ref_ch = next(((p, c) for p, c, lab in entries if lab == "free"), entries[0][:2])
        ref = residues(load(ref_pdb)[0][ref_ch])
        assert ref[ua_num].name in ("PHE", "TRP"), (prot, ref_pdb, ref[ua_num].name)
        win = list(range(ua_num + WIN[0], ua_num + WIN[1] + 1))
        core = list(range(ua_num + CORE[0], ua_num + CORE[1] + 1))
        chains = {}
        for pdb, ch, label in entries:
            chain = load(pdb)[0].find_chain(ch)
            if chain is None:
                print("missing", pdb, ch)
                continue
            chains[(pdb, ch)] = (mapped(residues(chain), ref), label)
        free = [k for k, (_, lab) in chains.items() if lab == "free"]
        for (pdb, ch), (res, label) in chains.items():
            st = load(pdb)
            ua = res.get(ua_num)
            R, t, ncore, rms = superpose(res, ref, core)
            srow = {"protein": prot, "pdb": pdb, "chain": ch, "label": label, "resolution": st.resolution,
                    "ref": f"{ref_pdb}:{ref_ch}", "core_n": ncore, "core_rmsd": round(rms, 2),
                    "window_modelled": sum(k in res for k in win),
                    "ua": f"{ua.name}{ua.seqid.num}" if ua else ""}
            if ua is not None and ring(ua) is not None and ring(ref[ua_num]) is not None:
                srow["ua_ring_shift_vs_ref"] = round(float(np.linalg.norm(ring(ua) @ R.T + t - ring(ref[ua_num]))), 2)
                srow["ua_chi1"] = dihedral(ua, ["N", "CA", "CB", "CG"])
                # Phe chi2 is 180-degree symmetric
                c2 = dihedral(ua, ["CA", "CB", "CG", "CD1"])
                srow["ua_chi2"] = (c2 + 90) % 180 - 90 if ua.name == "PHE" else c2
            # nearest free chain (a different entry), so nucleotide state is not read as a partner effect
            best = (np.inf, None, 0, None, None)
            for fk in free:
                if fk[0] == pdb:
                    continue
                r_, n_, R2, t2 = loop_rmsd(res, chains[fk][0], win, core)
                if r_ < best[0]:
                    best = (r_, fk, n_, R2, t2)
            if best[1] is not None:
                fres = chains[best[1]][0]
                srow.update({"nearest_free": ":".join(best[1]), "loop_rmsd_nearest_free": round(best[0], 2),
                             "loop_n_compared": best[2]})
                R2, t2 = best[3], best[4]
                if ua is not None and ring(ua) is not None and ring(fres.get(ua_num)) is not None if fres.get(ua_num) else False:
                    srow["ua_ring_shift_nearest_free"] = round(float(np.linalg.norm(ring(ua) @ R2.T + t2 - ring(fres[ua_num]))), 2)
                for k in win:
                    r = res.get(k)
                    row = {"protein": prot, "pdb": pdb, "chain": ch, "label": label, "pos": k - ua_num,
                           "resnum": k + UNIPROT_OFFSET.get(prot, 0), "auth": r.seqid.num if r else None,
                           "resname": r.name if r else "", "nearest_free": srow["nearest_free"]}
                    if r is not None and k in fres and ca(r) is not None and ca(fres[k]) is not None:
                        row["ca_shift"] = round(float(np.linalg.norm(ca(r) @ R2.T + t2 - ca(fres[k]))), 2)
                    frows.append(row)
            # partner contacts (deposited coordinates only)
            if label != "free":
                partners, ns = cognate_chains(st, ch)
                names = entity_names(pdb)
                touched, buried, pnames = 0, 0.0, set()
                hits_by = {}
                for k in win:
                    r = res.get(k)
                    if r is None:
                        continue
                    hits = set()
                    for atom in r:
                        for mk in ns.find_atoms(atom.pos, "\0", radius=CONTACT):
                            cra = mk.to_cra(st[0])
                            if cra.chain.name in partners and cra.residue.het_flag == "A" \
                                    and cra.atom.pos.dist(atom.pos) <= CONTACT:
                                hits.add((cra.chain.name, f"{cra.residue.name}{cra.residue.seqid.num}"))
                    hits_by[k] = hits
                pchains = sorted({c for h in hits_by.values() for c, _ in h})
                alone = sasa_by_res(pdb, [ch]) if pchains else {}
                bound = sasa_by_res(pdb, [ch] + pchains) if pchains else {}
                for k, hits in hits_by.items():
                    r = res[k]
                    n = r.seqid.num
                    d = round(alone.get(n, 0) - bound.get(n, 0), 1)
                    touched += bool(hits)
                    buried += d
                    for c, _ in hits:
                        pnames.add(names.get(partners[c], partners[c])[:40])
                    crows.append({"protein": prot, "pdb": pdb, "chain": ch, "label": label,
                                  "pos": k - ua_num, "resnum": k + UNIPROT_OFFSET.get(prot, 0), "auth": n,
                                  "resname": r.name, "sasa_free": round(alone.get(n, np.nan), 1) if pchains else np.nan,
                                  "d_sasa": d if pchains else 0.0,
                                  "partner_residues": "; ".join(f"{names.get(partners[c], c)[:20]}:{x}" for c, x in sorted(hits))})
                srow.update({"window_touched": touched, "window_buried_A2": round(buried, 0),
                             "loop_partners": "; ".join(sorted(pnames))})
            srows.append(srow)
    pd.DataFrame(crows).to_csv(RES / "loop_contacts.tsv", sep="\t", index=False)
    pd.DataFrame(frows).to_csv(RES / "loop_conformation.tsv", sep="\t", index=False)
    s = pd.DataFrame(srows)
    s.to_csv(RES / "loop_summary.tsv", sep="\t", index=False)
    with pd.option_context("display.width", 250, "display.max_colwidth", 40):
        print(s.drop(columns=["ref"]).to_string())


if __name__ == "__main__":
    main()
