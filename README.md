# The DEAD-box upstream anchor

Computational study of the conserved aromatic about 25 residues before the Q-motif
Gln of DEAD-box RNA helicases (Tanner et al., Mol Cell 2003): a buried Phe that
contacts adenine N6, and its Trp variant in four RNP-assembly subfamilies
(Prp5, Prp28, Mak5, Spb4). SC Lab, Dept. of Biochemistry, Mahidol University.

## Rebuild everything

```bash
./run_all.sh          # all stages + figures (several hours; stage 3 trees ~2 h on 8 threads)
./run_all.sh 3 4      # only the listed stages (each needs the earlier ones to have run)
```

Environment: micromamba env `helicase` (`env/helicase.yml`). Data versions and
access dates: `VERSIONS.md`. Large downloads and regenerable tables are not tracked
(see `.gitignore`); `run_all.sh` fetches them again.

## Stages

| Stage | Question | Key outputs |
|---|---|---|
| `stage0` | PDB: where is the anchor, what does it touch, does it move with partners? | `per_protein.tsv`, `n6_geometry.tsv`, `partner_contacts.tsv`, `loop_summary.tsv` |
| `stage1` | UniProt reference proteomes: anchor identity per subfamily and genus | `census_balanced_assigned.tsv` (regenerated), `genus_table.tsv` |
| `stage2` | AlphaFold models: is the called aromatic inserted at adenine N6? | `af_anchor.tsv`, `shell_covariation.tsv` |
| `stage3` | Phylogeny: ancestral state, switches, independent origins of Trp | `asr_summary.tsv`, `asr_switches.tsv`, `family_tree_trp.txt` |
| `stage4` | Human variants: ClinVar, gnomAD v4, cBioPortal | `variants_positions.tsv`, `constraint.tsv` |
| `figures` | Main-figure data panels (SC Lab style) and structure renderings | `figures/out/` |

Stage 0 family labels are hand-curated in `stage0/data/family_curation.tsv`; keep
that file under version control.

## Conventions

- The anchor (UA) is located from the Q-motif Gln: the F/W at Q-25, else the nearest
  F/W in Q-22..Q-30, else a Tyr (Rok1's anchor is the Phe at Q-29, confirmed in
  AlphaFold models).
- Partner contacts use deposited coordinates only; crystal-symmetry mates are excluded.
- Human variants are kept only where the transcript residue matches UniProt.
