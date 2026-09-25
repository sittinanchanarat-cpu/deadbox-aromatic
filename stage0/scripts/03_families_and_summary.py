#!/usr/bin/env python
"""Stage 0, step 3: family labels and the decision-gate table.

First run writes data/family_curation.tsv with keyword-based `family_auto`
and an empty `family_curated` column. Fill `family_curated` by hand wherever
the auto label is wrong or `?`; later runs keep your edits and only append
newly seen proteins. The summary uses family_curated when set, else family_auto.

Outputs
  data/family_curation.tsv    one row per UniProt (or description) -- EDIT THIS
  results/sites_labelled.tsv  sites.tsv + family
  results/decision_gate.tsv   per family: stacker frequency, identity, Q presence, spacers
"""
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA, RES = ROOT / "data", ROOT / "results"

# Order matters: first match wins.
RULES = [
    ("DEAH", r"\bDHX(8|15|16|32|33|35|38|40)\b|\bPRP(2|16|22|43)\b|PRPF(2|16|22|43)|DEAH|MLE\b|\bRHA\b|DHX9|DHX36|DHX29|DHX30|DHX57|DHX37|DHX34"),
    ("Ski2-like", r"\bMTR4\b|MTREX|SKIV2L|\bSKI2\b|BRR2|SNRNP200|HEL308|\bHJM\b|\bHEL(Q|S)\b|\bPOLQ\b|SUV3|SUPV3L1"),
    ("RLR/Dicer", r"RIG-?I|DDX58|RIGI|MDA5|IFIH1|LGP2|DHX58|DICER|DCR-?\d|DCL\d"),
    ("RecQ", r"RECQ|\bBLM\b|\bWRN\b|SGS1|RECQL"),
    ("viral", r"\bNS3\b|NPH-?II|NS3H|NONSTRUCTURAL PROTEIN 3|GENOME POLYPROTEIN|POLYPROTEIN"),
    ("DEAD", r"\bDDX\d+|DEAD|\bEIF-?4A|EIF4A|VASA|MSS116|DBP\d+|\bDED1\b|\bPRP(5|28)\b|PRPF28|UAP56|SUB2|RHLB|DBPA|CSHA|SRMB|HERA|\bVAD1\b|\bMJ0669\b|\bFAL1\b|\bROK1\b|\bHAS1\b|\bMAK5\b|\bSPB4\b|\bDRS1\b|RNA HELICASE"),
    ("Rad3/XPD", r"XPD|ERCC2|RAD3|DING|FANCJ|BRIP1|RTEL"),
    ("SWI2/SNF2", r"SNF2|SWI2|CHD\d|INO80|SWR1|RAD54|MOT1|BTAF1|ISWI|SMARCA"),
    ("other", r"RECG|MFD|UVRB|PRIA|HRQ|TOPOISOMERASE|REVERSE GYRASE|CAS3|ASCC3"),
]


def auto_family(text):
    t = (text or "").upper()
    for fam, pat in RULES:
        if re.search(pat, t):
            return fam
    return "?"


def main():
    ent = pd.read_csv(DATA / "entries.tsv", sep="\t", dtype=str).fillna("")
    ent["key"] = ent.uniprot.where(ent.uniprot != "", "desc:" + ent.description.str.upper())

    cur_path = DATA / "family_curation.tsv"
    prot = (ent.groupby("key")
               .agg(uniprot_name=("uniprot_name", "first"), description=("description", "first"),
                    organism=("organism", "first"), n_entries=("pdb", "nunique"),
                    example_pdbs=("pdb", lambda s: ",".join(sorted(set(s))[:5])))
               .reset_index())
    prot["family_auto"] = (prot.uniprot_name + " " + prot.description).map(auto_family)
    if cur_path.exists():
        old = pd.read_csv(cur_path, sep="\t", dtype=str).fillna("")
        prot = prot.merge(old[["key", "family_curated", "note"]], on="key", how="left").fillna("")
    else:
        prot["family_curated"], prot["note"] = "", ""
    prot.sort_values(["family_auto", "uniprot_name"]).to_csv(cur_path, sep="\t", index=False)
    fam = prot.set_index("key").apply(lambda r: r.family_curated or r.family_auto, axis=1)

    sites = pd.read_csv(RES / "sites.tsv", sep="\t")
    # map each site to the protein of the chain carrying the Walker K (fallback: stacker, Q)
    chain_key = {(r.pdb, c): r.key for r in ent.itertuples() for c in r.auth_chains.split(",")}

    def site_chain(r):
        for col in ("walker_lys", "stacker_chain", "q_gln"):
            v = r[col]
            if isinstance(v, str) and v:
                return v.split(":")[0]
        return None

    sites["chain"] = sites.apply(site_chain, axis=1)
    sites["key"] = [chain_key.get((p, c), "") for p, c in zip(sites.pdb, sites.chain)]
    sites["family"] = sites.key.map(fam).fillna("?")
    sites = sites[sites.family != "exclude"].copy()
    sites["protein"] = sites.key.map(prot.set_index("key").uniprot_name)
    sites["has_stacker"] = sites.has_recA1_stacker.astype(bool)
    # Q motif = Gln on adenine N6/N7 at the Q-motif register (WK-23 +/- 2); SNF2
    # remodellers have a Gln on adenine at WK-27, which is a different residue
    Q_REG = (21, 25)
    q_any = sites.q_gln.fillna("") != ""
    q_reg = sites.q_to_walkerK.between(*Q_REG)
    sites["has_q"] = q_any & (q_reg | sites.q_to_walkerK.isna())
    sites["stacker_aa"] = sites.stacker.fillna("").str[:3]
    ua_path = RES / "ua.tsv"
    if ua_path.exists():
        ua = pd.read_csv(ua_path, sep="\t").drop(columns="pdb")
        sites = sites.merge(ua, on="site", how="left")
    sites["ua_state"] = sites.get("ua_state", pd.Series(dtype=str)).fillna("no_walkerK")
    # the UA is the aromatic at WK-48 +/- 4 (Prp5 W257 = WK-49); inserted aromatics
    # at other registers (Hel308, RadD, Mtr4 at WK-55..60) are different residues
    UA_REG = (44, 52)
    sites["ua_inserted"] = sites.ua_state.eq("inserted") & sites.get("ua_to_walkerK", pd.Series(index=sites.index, dtype=float)).between(*UA_REG)
    # sites where the loop isn't in the model carry no evidence about the anchor
    sites["ua_informative"] = ~sites.ua_state.isin(["not_modelled", "no_walkerK"])
    sites["ua_aa_ins"] = sites.get("ua_aa", pd.Series(index=sites.index, dtype=str)).where(sites.ua_inserted, "").fillna("")
    sites.to_csv(RES / "sites_labelled.tsv", sep="\t", index=False)

    # Collapse to one row per protein so heavily crystallised proteins
    # (eIF4A, Mss116, Prp43 ...) don't dominate the counts.
    per_prot = (sites.groupby(["family", "key"])
                     .agg(protein=("protein", "first"), n_sites=("site", "size"),
                          frac_stacked=("has_stacker", "mean"), frac_q=("has_q", "mean"),
                          stacker_aa=("stacker_aa", lambda s: s[s != ""].mode().iat[0] if (s != "").any() else ""),
                          stacker_to_q=("stacker_to_q", "median"),
                          stacker_to_walkerK=("stacker_to_walkerK", "median"),
                          q_to_walkerK=("q_to_walkerK", "median"),
                          ua_informative=("ua_informative", "sum"),
                          ua_inserted=("ua_inserted", "sum"),
                          ua_aa=("ua_aa_ins", lambda s: s[s != ""].mode().iat[0] if (s != "").any() else ""),
                          ua_to_walkerK=("ua_to_walkerK", "median"),
                          ua_bfactor_z=("ua_bfactor_z", "median"),
                          loop_missing=("loop_missing", "median"))
                     .reset_index())
    per_prot.to_csv(RES / "per_protein.tsv", sep="\t", index=False)

    # a protein "has the UA" if the anchor is inserted in at least half of the
    # sites where the loop is actually modelled
    per_prot["ua_evaluable"] = per_prot.ua_informative > 0
    per_prot["has_ua"] = per_prot.ua_evaluable & (per_prot.ua_inserted >= 0.5 * per_prot.ua_informative)
    per_prot.to_csv(RES / "per_protein.tsv", sep="\t", index=False)
    gate = (per_prot.assign(stacked=per_prot.frac_stacked >= 0.5, q=per_prot.frac_q >= 0.5)
                    .groupby("family")
                    .agg(n_proteins=("key", "size"),
                         with_stacker=("stacked", "sum"),
                         with_q=("q", "sum"),
                         F=("stacker_aa", lambda s: (s == "PHE").sum()),
                         W=("stacker_aa", lambda s: (s == "TRP").sum()),
                         Y=("stacker_aa", lambda s: (s == "TYR").sum()),
                         H=("stacker_aa", lambda s: (s == "HIS").sum()),
                         med_stacker_to_q=("stacker_to_q", "median"),
                         med_stacker_to_walkerK=("stacker_to_walkerK", "median"),
                         min_stacker_to_walkerK=("stacker_to_walkerK", "min"),
                         max_stacker_to_walkerK=("stacker_to_walkerK", "max"),
                         ua_evaluable=("ua_evaluable", "sum"),
                         with_ua=("has_ua", "sum"),
                         ua_W=("ua_aa", lambda s: (s == "W").sum()),
                         ua_F=("ua_aa", lambda s: (s == "F").sum()),
                         ua_Y=("ua_aa", lambda s: (s == "Y").sum()),
                         ua_H=("ua_aa", lambda s: (s == "H").sum()),
                         med_ua_to_walkerK=("ua_to_walkerK", "median"))
                    .reset_index())
    gate.to_csv(RES / "decision_gate.tsv", sep="\t", index=False)
    print(gate.to_string(index=False))
    n_unl = (prot.family_auto.eq("?") & prot.family_curated.eq("")).sum()
    if n_unl:
        print(f"\n{n_unl} proteins still unlabelled -> edit {cur_path}")


if __name__ == "__main__":
    main()
