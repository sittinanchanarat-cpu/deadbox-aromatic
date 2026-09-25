#!/usr/bin/env python
"""Stage 0, step 1: find nucleotide-bound SF2 helicase structures in the PDB.

Query = (SF2 helicase Pfam on any polymer entity) AND (adenine nucleotide or
analogue bound in the entry). Writes entry metadata and downloads mmCIFs.

Outputs
  data/entries.tsv   one row per polymer entity carrying an SF2 Pfam
  data/cif/*.cif.gz  coordinates
"""
import argparse, sys, time
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
CIF = DATA / "cif"

SF2_PFAMS = {
    "PF00270": "DEAD",        # DEAD/DEAH box helicase (DEAD, DEAH, Ski2, RecQ, ...)
    "PF04851": "ResIII",      # type III restriction enzyme res subunit (RIG-I-like, ...)
    "PF07652": "Flavi_DEAD",  # flavivirus NS3 helicase
    "PF00271": "Helicase_C",  # C-terminal helicase domain (catches the rest)
}
# Adenine nucleotides and common analogues. Transition-state mimics (ADP-AlF4,
# ADP-BeF3, ADP-VO4, ADP-MgF) are caught through ADP.
LIGANDS = ["ATP", "ADP", "ANP", "ACP", "AGS", "APC", "AMP", "ADX", "AN2", "ATPS"]

SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL = "https://data.rcsb.org/graphql"


def search_entries():
    pfam_nodes = [
        {"type": "terminal", "service": "text",
         "parameters": {"attribute": "rcsb_polymer_entity_annotation.annotation_id",
                        "operator": "exact_match", "value": p}}
        for p in SF2_PFAMS]
    query = {
        "query": {"type": "group", "logical_operator": "and", "nodes": [
            {"type": "group", "logical_operator": "or", "nodes": pfam_nodes},
            {"type": "terminal", "service": "text",
             "parameters": {"attribute": "rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id",
                            "operator": "in", "value": LIGANDS}},
        ]},
        "return_type": "entry",
        "request_options": {"return_all_hits": True},
    }
    r = requests.post(SEARCH, json=query, timeout=60)
    r.raise_for_status()
    return sorted(h["identifier"] for h in r.json()["result_set"])


GQL = """
query($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    struct { title }
    exptl { method }
    rcsb_entry_info { resolution_combined }
    rcsb_accession_info { initial_release_date }
    nonpolymer_entities { pdbx_entity_nonpoly { comp_id } }
    polymer_entities {
      rcsb_id
      rcsb_polymer_entity { pdbx_description }
      rcsb_polymer_entity_container_identifiers { auth_asym_ids asym_ids }
      rcsb_entity_source_organism { scientific_name ncbi_taxonomy_id }
      uniprots { rcsb_id rcsb_uniprot_protein { name { value } } }
      rcsb_polymer_entity_annotation { annotation_id type }
    }
  }
}"""


def fetch_metadata(ids):
    rows = []
    for i in range(0, len(ids), 50):
        r = requests.post(GRAPHQL, json={"query": GQL, "variables": {"ids": ids[i:i + 50]}}, timeout=120)
        r.raise_for_status()
        for e in r.json()["data"]["entries"]:
            ligs = sorted({n["pdbx_entity_nonpoly"]["comp_id"] for n in (e["nonpolymer_entities"] or [])})
            res = (e["rcsb_entry_info"]["resolution_combined"] or [None])[0]
            for pe in e["polymer_entities"]:
                ann = pe["rcsb_polymer_entity_annotation"] or []
                pfams = sorted({a["annotation_id"] for a in ann if a["type"] == "Pfam"})
                sf2 = [p for p in pfams if p in SF2_PFAMS]
                if not sf2:
                    continue
                ids_ = pe["rcsb_polymer_entity_container_identifiers"]
                org = (pe["rcsb_entity_source_organism"] or [{}])[0]
                up = (pe["uniprots"] or [{}])[0]
                rows.append({
                    "pdb": e["rcsb_id"],
                    "entity": pe["rcsb_id"],
                    "auth_chains": ",".join(ids_["auth_asym_ids"] or []),
                    "label_chains": ",".join(ids_["asym_ids"] or []),
                    "description": pe["rcsb_polymer_entity"]["pdbx_description"],
                    "uniprot": up.get("rcsb_id", ""),
                    "uniprot_name": ((up.get("rcsb_uniprot_protein") or {}).get("name") or {}).get("value", ""),
                    "organism": org.get("scientific_name", ""),
                    "taxid": org.get("ncbi_taxonomy_id", ""),
                    "sf2_pfams": ",".join(sf2),
                    "all_pfams": ",".join(pfams),
                    "method": e["exptl"][0]["method"],
                    "resolution": res,
                    "released": e["rcsb_accession_info"]["initial_release_date"][:10],
                    "title": e["struct"]["title"],
                    "ligands": ",".join(ligs),
                })
        time.sleep(0.2)
    return pd.DataFrame(rows)


def download(pdb):
    out = CIF / f"{pdb.lower()}.cif.gz"
    if out.exists() and out.stat().st_size > 0:
        return
    r = requests.get(f"https://files.rcsb.org/download/{pdb}.cif.gz", timeout=300)
    r.raise_for_status()
    out.write_bytes(r.content)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metadata-only", action="store_true", help="write entries.tsv, download nothing")
    ap.add_argument("--ids", help="comma-separated PDB IDs: download only these (pilot)")
    args = ap.parse_args()
    CIF.mkdir(parents=True, exist_ok=True)
    if (DATA / "entries.tsv").exists() and args.ids:
        df = pd.read_csv(DATA / "entries.tsv", sep="\t", dtype=str)
    else:
        ids = search_entries()
        print(f"{len(ids)} entries match", file=sys.stderr)
        df = fetch_metadata(ids)
        df.to_csv(DATA / "entries.tsv", sep="\t", index=False)
        print(f"{len(df)} SF2 polymer entities in {df.pdb.nunique()} entries", file=sys.stderr)
    if args.metadata_only:
        return
    todo = [i.strip().upper() for i in args.ids.split(",")] if args.ids else list(df.pdb.unique())
    for n, pdb in enumerate(todo, 1):
        try:
            download(pdb)
        except Exception as ex:
            print(f"  download failed {pdb}: {ex}", file=sys.stderr)
        if n % 25 == 0:
            print(f"  {n} downloaded", file=sys.stderr)


if __name__ == "__main__":
    main()
