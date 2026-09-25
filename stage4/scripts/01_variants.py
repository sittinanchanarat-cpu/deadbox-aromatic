#!/usr/bin/env python
"""Human variants at the upstream anchor and its shell, from ClinVar and gnomAD v4.

Genes: every reviewed human DEAD-box protein in the Stage 1 census, plus those the
InterPro Q-motif profile misses (EXTRA), located with the same rule (01_census.calls).
Positions per protein (UniProt numbering):
  UA      the anchor (F/W)            UA+1   exposed residue after it (partner hot spot)
  Q-7/6   Q-motif stacker on adenine  Q-4    Q-motif Pro   Q  the Q-motif Gln
  loop    UA-8..UA+3 (anchor loop, excluding the above)
Source: the gnomAD GraphQL API (gene -> clinvar_variants, variants in gnomad_r4),
canonical Ensembl transcript. A variant's protein position is taken from hgvsp and
kept only if its reference residue matches the UniProt sequence at that position
(otherwise the transcript numbering differs from UniProt and it is flagged).

Output
  data/genes.tsv                 gene, accession, positions, UniProt sequence
  data/gnomad/<gene>.json        raw API responses (cached)
  results/variants_positions.tsv every ClinVar or gnomAD variant at a listed position
  results/gnomad_missense_counts.tsv  per gene: missense allele count at UA vs loop vs all
"""
import json
import re
import sys
import time
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
S1 = ROOT.parent / "stage1"
sys.path.insert(0, str(S1 / "scripts"))
from importlib import import_module  # noqa: E402

census_calls = import_module("01_census").calls
EXTRA = {"Q92499": "DDX1", "Q9NVP1": "DDX18", "Q9NR30": "DDX21", "Q9BQ39": "DDX50",
         "Q8N8A6": "DDX51", "Q86TM3": "DDX53"}
# Q-motif Gln set by hand where UniProt has no usable Q-motif span:
# DDX1 (MAAF...DIQ29: anchor F4 = Q-25); DDX51 (...SSYFPVQ225; no aromatic at Q-22..Q-30)
MANUAL_Q = {"Q92499": 29, "Q8N8A6": 225}
AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
       "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
       "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "*"}
QUERY = """{ gene(gene_symbol: "%s", reference_genome: GRCh38) { gene_id canonical_transcript_id
  clinvar_variants { variant_id hgvsp clinical_significance review_status major_consequence }
  variants(dataset: gnomad_r4) { variant_id hgvsp consequence exome { ac an } genome { ac an } } } }"""


def uniprot(acc):
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.tsv?fields=accession,gene_primary,ft_motif,sequence"
    with urllib.request.urlopen(url, timeout=60) as r:
        _, row = r.read().decode().strip().split("\n")
    a, gene, motif, seq = row.split("\t")
    return gene, motif, seq


def genes():
    c = pd.read_csv(S1 / "results" / "census_all.tsv", sep="\t", low_memory=False)
    h = c[(c.taxid == 9606) & (c.Reviewed == "reviewed")]
    rows = []
    accs = list(h.acc) + [a for a in EXTRA if a not in set(h.acc)]
    for acc in accs:
        gene, motif, seq = uniprot(acc)
        if acc in MANUAL_Q:  # a synthetic span ending at the Gln reuses the census rule
            q = MANUAL_Q[acc]
            motif = f'MOTIF {q - 27}..{q + 1}; /note="Q motif"'
        call = census_calls(pd.Series({"Motif": motif, "Sequence": seq}))
        if pd.isna(call.get("ua_pos")) or not call.get("ua_found", False):
            print(f"{gene}: no anchor (no F/W/Y at Q-22..Q-30); skipped", file=sys.stderr)
            continue
        ua, q = int(call.ua_pos), int(call.q_pos)
        stack = next((q - o for o in (7, 6) if seq[q - o - 1] in "FWYH"), q - 7)
        pos = {"UA": ua, "UA+1": ua + 1, "Q-stacker": stack, "Q-4": q - 4, "Q": q}
        for k in range(ua - 8, ua + 4):
            if k not in pos.values():
                pos[f"loop{k - ua:+d}"] = k
        rows.append({"acc": acc, "gene": gene, "ua": f"{seq[ua - 1]}{ua}", "q": q, "seq": seq,
                     "positions": json.dumps(pos)})
    g = pd.DataFrame(rows)
    g.to_csv(ROOT / "data" / "genes.tsv", sep="\t", index=False)
    return g


def gnomad(gene):
    out = ROOT / "data" / "gnomad" / f"{gene}.json"
    if out.exists():
        return json.loads(out.read_text())
    req = urllib.request.Request("https://gnomad.broadinstitute.org/api",
                                 data=json.dumps({"query": QUERY % gene}).encode(),
                                 headers={"Content-Type": "application/json"})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=180) as r:
                d = json.load(r)
            if "errors" in d:
                raise RuntimeError(d["errors"][0]["message"])
            out.write_text(json.dumps(d))
            time.sleep(6)  # stay under the API rate limit
            return d
        except Exception as ex:
            print(f"{gene}: {ex}; retry", file=sys.stderr)
            time.sleep(30 * (attempt + 1))
    return None


def parse_hgvsp(h):
    m = re.match(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|=|\*|fs|del|dup)?", h or "")
    if not m:
        return None, None, None
    return AA3.get(m.group(1)), int(m.group(2)), m.group(3)


def main():
    (ROOT / "data" / "gnomad").mkdir(parents=True, exist_ok=True)
    g = genes()
    hits, counts = [], []
    for r in g.itertuples():
        d = gnomad(r.gene)
        if d is None or d["data"]["gene"] is None:
            print(f"{r.gene}: no gnomAD data", file=sys.stderr)
            continue
        gd = d["data"]["gene"]
        pos = json.loads(r.positions)
        where = {v: k for k, v in pos.items()}
        mism = 0
        for src, key in (("ClinVar", "clinvar_variants"), ("gnomAD", "variants")):
            for v in gd[key]:
                ref, p, alt = parse_hgvsp(v.get("hgvsp"))
                if p is None or p not in where:
                    continue
                ok = p <= len(r.seq) and r.seq[p - 1] == ref
                mism += not ok
                ac = sum((v.get(x) or {}).get("ac", 0) for x in ("exome", "genome")) if src == "gnomAD" else None
                hits.append({"gene": r.gene, "acc": r.acc, "anchor": r.ua, "site": where[p], "pos": p,
                             "source": src, "hgvsp": v["hgvsp"], "consequence": v.get("major_consequence") or v.get("consequence"),
                             "clinvar": v.get("clinical_significance", ""), "review": v.get("review_status", ""),
                             "gnomad_ac": ac, "uniprot_match": ok, "variant_id": v["variant_id"]})
        # gnomAD missense burden (UniProt-matched variants only): anchor vs loop vs whole protein
        mis = []
        for v in gd["variants"]:
            if v.get("consequence") != "missense_variant" or not v.get("hgvsp"):
                continue
            ref, p, _ = parse_hgvsp(v["hgvsp"])
            if p is not None and p <= len(r.seq) and r.seq[p - 1] == ref:  # same filter as the site hits
                mis.append((p, sum((v.get(x) or {}).get("ac", 0) for x in ("exome", "genome"))))
        mis = pd.DataFrame(mis, columns=["p", "ac"])
        ua = pos["UA"]
        counts.append({"gene": r.gene, "anchor": r.ua, "transcript": gd["canonical_transcript_id"],
                       "missense_variants_total": len(mis), "length": len(r.seq),
                       "missense_sites_per_100aa": round(100 * mis.p.nunique() / len(r.seq), 1),
                       "missense_at_UA": int((mis.p == ua).sum()), "missense_AC_at_UA": int(mis.ac[mis.p == ua].sum()),
                       "missense_in_loop": int(mis.p.between(ua - 8, ua + 3).sum()),
                       "numbering_mismatches": mism})
    h = pd.DataFrame(hits)
    h.to_csv(ROOT / "results" / "variants_positions.tsv", sep="\t", index=False)
    pd.DataFrame(counts).to_csv(ROOT / "results" / "gnomad_missense_counts.tsv", sep="\t", index=False)
    print(f"{len(g)} genes, {len(h)} variants at listed positions", file=sys.stderr)


if __name__ == "__main__":
    main()
