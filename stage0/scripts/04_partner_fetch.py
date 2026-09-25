#!/usr/bin/env python
"""Partner-complex census, step 1: DEAD-box proteins in complex with anything.

Tanner et al. (2003) proposed that partner proteins bound near the N-terminus
could regulate DEAD-box ATPase through the Q motif + upstream aromatic. This
collects every PDB entry with a DEAD-box chain plus at least one other polymer
(protein or nucleic acid), nucleotide-bound or not.

DEAD-box entities are called from sequence, not annotation: a Q motif
(..P..[IVLM]Q) 12-40 residues before motif I (GKT/GKS) and a DEAD/DEAH-type
motif II after it, with motif II reading DExD (DEAD family, incl. variants such as DEVD). The upstream
anchor (UA) is the aromatic 21-30 residues before the Q, nearest to Q-25.

Outputs
  data/partner_entities.tsv   one row per polymer entity in the kept entries
  data/partner_cif/*.cif.gz
"""
import argparse, re, sys, time
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
CIF = DATA / "partner_cif"
SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL = "https://data.rcsb.org/graphql"
ANNOT = ["PF00270", "IPR011545", "IPR014014", "IPR000629"]

MOTIF_I = re.compile(r"[AGST]..[GST][ST]GK[ST]")
Q_MOTIF = re.compile(r"(?:P..|..P)[IVLMF]Q")   # Pro at Q-4 or Q-2 (Prp5 has Leu at Q-4)
MOTIF_II = re.compile(r"DE.[DH]")    # DExD covers DEAD variants (DDX21: DEVD)


def locate(seq):
    """Return (q_idx, ua_idx, ua_aa, motif_II) 0-based, or None if not DEAD-like."""
    for m1 in MOTIF_I.finditer(seq):
        k = m1.start() + 6                      # the K of GKT
        qs = [m.end() - 1 for m in Q_MOTIF.finditer(seq, max(0, k - 45), k)
              if 12 <= k - (m.end() - 1) <= 40]
        if not qs:
            continue
        # canonical Q sits 23 residues before the Walker K; take the candidate
        # closest to that (human Prp28 has a second "..PIGLQ" 14 before the K)
        q = min(qs, key=lambda i: abs((k - i) - 23))
        m2 = MOTIF_II.search(seq, k + 60, k + 200)
        if not m2:
            continue
        ua = [i for i in range(max(0, q - 30), q - 20) if seq[i] in "FWYH"]
        ua_i = min(ua, key=lambda i: abs((q - i) - 25)) if ua else None
        return q, ua_i, (seq[ua_i] if ua_i is not None else ""), m2.group()
    return None


def search_ids():
    nodes = [{"type": "terminal", "service": "text",
              "parameters": {"attribute": "rcsb_polymer_entity_annotation.annotation_id",
                             "operator": "exact_match", "value": a}} for a in ANNOT]
    q = {"query": {"type": "group", "logical_operator": "or", "nodes": nodes},
         "return_type": "entry", "request_options": {"return_all_hits": True}}
    r = requests.post(SEARCH, json=q, timeout=60)
    r.raise_for_status()
    return sorted(h["identifier"] for h in r.json()["result_set"])


GQL = """
query($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    exptl { method }
    rcsb_entry_info { resolution_combined deposited_atom_count }
    struct { title }
    nonpolymer_entities { pdbx_entity_nonpoly { comp_id } }
    polymer_entities {
      rcsb_id
      entity_poly { pdbx_seq_one_letter_code_can rcsb_entity_polymer_type }
      rcsb_polymer_entity { pdbx_description }
      rcsb_polymer_entity_container_identifiers { auth_asym_ids }
      uniprots { rcsb_id }
      rcsb_entity_source_organism { scientific_name }
    }
  }
}"""


def metadata(ids):
    rows = []
    for i in range(0, len(ids), 40):
        for attempt in range(3):
            try:
                r = requests.post(GRAPHQL, json={"query": GQL, "variables": {"ids": ids[i:i + 40]}}, timeout=180)
                r.raise_for_status()
                break
            except Exception:
                time.sleep(5)
        for e in r.json()["data"]["entries"]:
            if not e or not e.get("polymer_entities"):
                continue
            info = e["rcsb_entry_info"] or {}
            ligs = ",".join(sorted({n["pdbx_entity_nonpoly"]["comp_id"] for n in (e["nonpolymer_entities"] or [])}))
            for pe in e["polymer_entities"]:
                seq = pe["entity_poly"]["pdbx_seq_one_letter_code_can"] or ""
                ptype = pe["entity_poly"]["rcsb_entity_polymer_type"]
                loc = locate(seq) if ptype == "Protein" else None
                rows.append({
                    "pdb": e["rcsb_id"], "entity": pe["rcsb_id"], "type": ptype,
                    "chains": ",".join(pe["rcsb_polymer_entity_container_identifiers"]["auth_asym_ids"] or []),
                    "description": pe["rcsb_polymer_entity"]["pdbx_description"],
                    "uniprot": ((pe["uniprots"] or [{}])[0]).get("rcsb_id", ""),
                    "organism": ((pe["rcsb_entity_source_organism"] or [{}])[0]).get("scientific_name", ""),
                    "is_dead": bool(loc and loc[3].endswith("D")),   # Q motif + DExD
                    "is_dead_like": bool(loc),
                    "motif_II": loc[3] if loc else "",
                    "q_seq": loc[0] + 1 if loc else None,          # 1-based = label_seq
                    "ua_seq": (loc[1] + 1) if loc and loc[1] is not None else None,
                    "ua_aa": loc[2] if loc else "",
                    "length": len(seq),
                    "method": (e["exptl"] or [{}])[0].get("method", ""),
                    "resolution": (info.get("resolution_combined") or [None])[0],
                    "atoms": info.get("deposited_atom_count") or 0,
                    "ligands": ligs, "title": (e["struct"] or {}).get("title", ""),
                })
        time.sleep(0.2)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--download", action="store_true")
    args = ap.parse_args()
    ids = search_ids()
    print(f"{len(ids)} entries with DEAD/DEAH annotations", file=sys.stderr)
    df = metadata(ids)
    dead_entries = set(df[df.is_dead].pdb)
    n_poly = df.groupby("pdb").entity.nunique()
    keep = sorted(p for p in dead_entries if n_poly[p] > 1)
    df = df[df.pdb.isin(keep)]
    DATA.mkdir(exist_ok=True)
    df.to_csv(DATA / "partner_entities.tsv", sep="\t", index=False)
    atoms = df.drop_duplicates("pdb").atoms.sum()
    print(f"{len(dead_entries)} entries with a DEAD-box chain; {len(keep)} have another polymer; "
          f"{atoms/1e6:.1f} M atoms total", file=sys.stderr)
    if args.download:
        CIF.mkdir(exist_ok=True)
        for n, p in enumerate(keep, 1):
            out = CIF / f"{p.lower()}.cif.gz"
            if out.exists() and out.stat().st_size:
                continue
            try:
                r = requests.get(f"https://files.rcsb.org/download/{p}.cif.gz", timeout=600)
                r.raise_for_status()
                out.write_bytes(r.content)
            except Exception as ex:
                print(f"  {p}: {ex}", file=sys.stderr)
            if n % 50 == 0:
                print(f"  {n}/{len(keep)}", file=sys.stderr)


if __name__ == "__main__":
    main()
