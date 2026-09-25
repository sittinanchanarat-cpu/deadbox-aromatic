#!/usr/bin/env python
"""Stage 2, step 2: is the upstream aromatic inserted at the adenine in AlphaFold models?

AlphaFold models carry no nucleotide, so one is transplanted: the reference
Vasa-AMPPNP complex (2DB3 chain A) is superposed on each model using RecA1 motif
C-alpha atoms only (Q motif Q-8..Q+1, motif I GxGK[TS], motif II DExD; 19 atoms),
and its adenine is carried into the model frame. The Q-motif stacker (aromatic at
Q-7/Q-6) is then checked against the transplanted adenine as a sanity test.

For every aromatic (F/W/Y) at Q-30..Q-20 of each model:
  n6_centroid      N6 -> ring centroid (Stage 0 PDB, inserted Phe: 5.1 +/- 0.34 A)
  n6_closest       closest ring atom to N6 and its distance (PDB: CZ, 3.3-4.1 A)
  approach_angle   ring normal vs centroid->N6 (PDB: 82 +/- 5 deg, edge-on)
  to_qarom         ring centroid -> Q-motif stacker ring centroid (PDB inserted 6.3-7.1 A)
  rsasa            relative side-chain-inclusive SASA in the model (Tien 2013 maxima)
  plddt            pLDDT of the residue (CA B-factor)
  inserted         n6_closest <= 5.0 A and n6_centroid <= 7.0 A

Per protein, the structural anchor is the inserted aromatic with the shortest
n6_closest; it is compared with the sequence call from Stage 1.

Output
  results/af_aromatics.tsv   one row per aromatic in the window
  results/af_anchor.tsv      one row per protein
"""
import gzip
import io
import re
import sys
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser
from Bio.PDB.SASA import ShrakeRupley

ROOT = Path(__file__).resolve().parents[1]
STAGE0 = ROOT.parent / "stage0"
REF_PDB, REF_CHAIN, REF_Q, REF_LIG = "2db3", "A", 272, "ANP"
RINGS = {
    "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TRP": ["CG", "CD1", "NE1", "CE2", "CD2", "CE3", "CZ2", "CZ3", "CH2"],
}
ADENINE = ["N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4"]
MAX_ASA = {"PHE": 240.0, "TRP": 285.0, "TYR": 263.0}
ONE = {"PHE": "F", "TRP": "W", "TYR": "Y"}
# Q-motif Gln for references missing from the census (no InterPro Q-motif match)
Q_FALLBACK = {"P23394": 201}  # yeast Prp28, curated Q motif 172..202


def seq_of(res):
    return "".join(gemmi.find_tabulated_residue(r.name).one_letter_code.upper() for r in res)


def motif_idx(res, q_idx):
    """Indices (into res) of the 19 superposition residues, or None."""
    s = seq_of(res)
    m1 = re.compile(r"G.GK[TS]").search(s, q_idx + 10, q_idx + 70)
    if not m1:
        return None
    m2 = re.compile(r"DE.D").search(s, m1.end() + 60, m1.end() + 300)
    if not m2:
        return None
    idx = list(range(q_idx - 8, q_idx + 2)) + list(range(m1.start(), m1.end())) + list(range(m2.start(), m2.end()))
    return idx if min(idx) >= 0 else None


def ca(r):
    return np.array(r["CA"][0].pos.tolist())


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - pc @ R.T


def ring(r):
    try:
        xyz = np.array([r[a][0].pos.tolist() for a in RINGS[r.name]])
    except (RuntimeError, IndexError, KeyError):
        return None, None, None
    c = xyz.mean(0)
    return c, np.linalg.svd(xyz - c)[2][2], xyz


def load_reference():
    st = gemmi.read_structure(str(STAGE0 / "data" / "cif" / f"{REF_PDB}.cif.gz"))
    ch = st[0][REF_CHAIN]
    res = [r for r in ch if r.het_flag == "A"]
    q_idx = next(i for i, r in enumerate(res) if r.seqid.num == REF_Q)
    idx = motif_idx(res, q_idx)
    lig = next(r for r in ch if r.name == REF_LIG)
    ade = {a: np.array(lig[a][0].pos.tolist()) for a in ADENINE}
    return np.array([ca(res[i]) for i in idx]), ade


def sasa(path):
    """Per-residue SASA (A^2) of the model, keyed by residue number."""
    with gzip.open(path, "rt") as fh:
        s = MMCIFParser(QUIET=True).get_structure("m", io.StringIO(fh.read()))
    ShrakeRupley().compute(s[0], level="R")
    return {r.id[1]: r.sasa for r in s[0].get_residues()}


def analyse(acc, q_pos, ref_ca, ref_ade):
    path = ROOT / "data" / "af" / f"AF-{acc}.cif.gz"
    st = gemmi.read_structure(str(path))
    res = [r for r in st[0][0] if r.het_flag == "A"]
    num = {r.seqid.num: i for i, r in enumerate(res)}
    if q_pos not in num or res[num[q_pos]].name != "GLN":
        return None, "q_mismatch"
    q_idx = num[q_pos]
    idx = motif_idx(res, q_idx)
    if idx is None:
        return None, "motif_not_found"
    P = np.array([ca(res[i]) for i in idx])
    R, t = kabsch(ref_ca, P)                                   # reference -> model frame
    rms = float(np.sqrt(((ref_ca @ R.T + t - P) ** 2).sum(1).mean()))
    ade = {a: x @ R.T + t for a, x in ref_ade.items()}
    n6, ade_c = ade["N6"], np.mean([ade[a] for a in ADENINE[1:]], 0)
    # Q-motif stacker sanity check
    qa = next((res[q_idx - o] for o in (7, 6) if res[q_idx - o].name in RINGS), None)
    qa_c = ring(qa)[0] if qa is not None else None
    qa_to_ade = float(np.linalg.norm(qa_c - ade_c)) if qa_c is not None else np.nan
    asa = sasa(path)
    rows = []
    for off in range(20, 31):
        i = q_idx - off
        if i < 0 or res[i].name not in RINGS:
            continue
        r = res[i]
        c, nrm, xyz = ring(r)
        v = n6 - c
        k = int(np.argmin(np.linalg.norm(xyz - n6, axis=1)))
        rows.append({
            "acc": acc, "resnum": r.seqid.num, "aa": ONE[r.name], "q_offset": off,
            "n6_centroid": round(float(np.linalg.norm(v)), 2),
            "n6_closest_atom": RINGS[r.name][k], "n6_closest": round(float(np.linalg.norm(xyz[k] - n6)), 2),
            "approach_angle": round(float(np.degrees(np.arccos(abs(v @ nrm) / np.linalg.norm(v)))), 1),
            "to_qarom": round(float(np.linalg.norm(c - qa_c)), 2) if qa_c is not None else np.nan,
            "rsasa": round(asa.get(r.seqid.num, np.nan) / MAX_ASA[r.name], 3),
            "plddt": round(r["CA"][0].b_iso, 1),
        })
    return {"motif_rmsd": round(rms, 2), "qarom": f"{qa.name}{qa.seqid.num}" if qa is not None else "",
            "qarom_to_adenine": round(qa_to_ade, 2), "rows": rows}, "ok"


def main():
    sel = pd.read_csv(ROOT / "data" / "selection.tsv", sep="\t")
    allc = pd.read_csv(ROOT.parent / "stage1" / "results" / "census_all.tsv", sep="\t",
                       usecols=["acc", "q_pos"]).set_index("acc").q_pos
    sel["q_pos"] = sel.q_pos.fillna(sel.acc.map(allc)).fillna(sel.acc.map(Q_FALLBACK))
    ref_ca, ref_ade = load_reference()
    arom, prot = [], []
    for s in sel.itertuples():
        base = {"acc": s.acc, "subfamily": s.subfamily, "label": s.label if isinstance(s.label, str) else "",
                "organism": s.organism, "class": sel.at[s.Index, "class"], "census_anchor": s.anchor,
                "census_pos": s.ua_pos}
        if not str(s.af).startswith(("ok", "cached")):
            prot.append({**base, "status": "no_model"})
            continue
        if pd.isna(s.q_pos):
            prot.append({**base, "status": "no_q"})
            continue
        try:
            out, status = analyse(s.acc, int(s.q_pos), ref_ca, ref_ade)
        except Exception as ex:
            out, status = None, f"error:{ex}"
        if out is None:
            prot.append({**base, "status": status})
            continue
        rows = pd.DataFrame(out["rows"])
        if len(rows):
            rows["inserted"] = (rows.n6_closest <= 5.0) & (rows.n6_centroid <= 7.0)
            arom.append(rows.assign(subfamily=s.subfamily))
            ins = rows[rows.inserted].sort_values("n6_closest")
        else:
            ins = rows
        a = ins.iloc[0] if len(ins) else None
        prot.append({**base, "status": "ok", "motif_rmsd": out["motif_rmsd"], "qarom": out["qarom"],
                     "qarom_to_adenine": out["qarom_to_adenine"],
                     "struct_anchor": f"{a.aa}{a.resnum}" if a is not None else "none",
                     "struct_aa": a.aa if a is not None else "-", "struct_offset": a.q_offset if a is not None else np.nan,
                     "matches_census": (a is not None and a.resnum == s.ua_pos),
                     **({k: a[k] for k in ["n6_centroid", "n6_closest_atom", "n6_closest", "approach_angle",
                                           "to_qarom", "rsasa", "plddt"]} if a is not None else {})})
    pd.concat(arom).to_csv(ROOT / "results" / "af_aromatics.tsv", sep="\t", index=False)
    p = pd.DataFrame(prot)
    p.to_csv(ROOT / "results" / "af_anchor.tsv", sep="\t", index=False)
    print(p.status.value_counts().to_dict(), file=sys.stderr)


if __name__ == "__main__":
    main()
