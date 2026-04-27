# Handoff: Local AGP QZA SILVA Transfer Rerun

Date: 2026-04-26

This handoff documents the work done to replace the PRIME-derived AGP
taxonomy source with the local AGP feature table, recover the corresponding
ASV sequences from redbiom, assign taxonomy with the same DADA2/SILVA pipeline
used for validation datasets, rebuild AGP dCST labels, rerun external
taxon-label transfer, and rerun downstream AGP-dCST/outcome association
analyses.

The short version: the local AGP data do contain md5 ASV feature IDs. I
recovered a complete md5-to-sequence map for those ASVs from the redbiom
`Deblur_2021.09-Illumina-16S-V4-100nt-50b3a2` context, assigned SILVA taxonomy
with the validation taxonomy code path, rebuilt the AGP absorb hierarchy, and
reran transfer for all validation comparisons. The transfer now has perfect
sample retention at depths 1-4 for all 15 validation comparisons.

## Why This Was Done

The previous paper-facing transfer analysis used a PRIME/project AGP taxonomy
table. That table was already taxon-level and did not expose ASV sequences, so
the AGP discovery labels and validation labels were not assigned by exactly the
same sequence-to-taxonomy procedure.

The concern was that the transfer script matched validation taxon labels into
an AGP taxon dictionary. If the AGP dictionary and validation dictionaries were
made by different taxonomic assignment workflows, then some validation samples
could lose all usable AGP-label support at deeper dCST depths. That sample loss
could bias the direct-transfer validation analyses.

The goal was therefore:

1. use the local AGP filtered count table as the abundance source;
2. recover actual ASV sequences corresponding to the AGP md5 feature IDs;
3. assign taxonomy to AGP ASVs with the same DADA2/SILVA references and label
   formatting used for validation datasets;
4. collapse the AGP ASV count table to those SILVA labels;
5. rebuild AGP dCST assignments from that collapsed table;
6. rerun taxon-label transfer for all validation datasets;
7. verify that no samples are lost at any depth;
8. rerun all AGP-dCST/outcome association analyses under the new transfer.

## Repositories And Main Roots

The paper-facing repository is:

`/Users/pgajer/current_projects/linf`

The analysis-scale repository and outputs are:

`/Users/pgajer/current_projects/gut_microbiome`

The local AGP source dataset is:

`/Users/pgajer/current_projects/AGP`

The main output roots generated or updated by this work are:

`/Users/pgajer/current_projects/gut_microbiome/outputs/agp_silva_taxonomy/2026-04-26-redbiom-md5-sequence-map`

`/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_analysis/runs/2026-04-26-agp-silva-local-qza-absorb-depth4`

`/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_validation/frozen_agp_silva_local_qza_2026-04-26`

## Source Data Findings

### Local AGP Data

The local AGP QIIME 2 artifact is:

`/Users/pgajer/current_projects/AGP/data/raw.nobloom.minfeat.mindepth.biom.qza`

It is a `FeatureTable[Frequency]` artifact with:

- 24,566 samples;
- 286,255 md5-like ASV feature IDs;
- Qiita-style sample IDs such as `10317.x`.

The matching local AGP metadata file is:

`/Users/pgajer/current_projects/AGP/data/All_good.tsv`

The preprocessing notes are:

`/Users/pgajer/current_projects/AGP/docs/AGP_preprocessing.tex`

Those notes identify the AGP redbiom context as:

`Deblur_2021.09-Illumina-16S-V4-100nt-50b3a2`

The notes also identify Qiita study `10317`, 100 nt ASVs, and bloom-filtered
minimum-depth AGP count processing.

### PRIME AGP Data

The PRIME/project AGP dataset available to the paper workflow was taxon-level.
I did not find usable ASV sequences there, which means it cannot support a true
rerun of taxonomy assignment through the same sequence-level validation
pipeline.

## Redbiom Search And Sequence Recovery

I used a temporary environment with `redbiom` and `biom-format`:

`/tmp/linf-redbiom-venv`

The exact redbiom context exists:

`Deblur_2021.09-Illumina-16S-V4-100nt-50b3a2`

`redbiom summarize contexts` reported approximately 313,690 samples and
8,625,124 features for this context.

I also confirmed that related Deblur V4 contexts exist at 90, 100, 125, 150,
200, 250, and 300 nt. The 100 nt context is the one matching the local AGP
preprocessing notes.

Important behavior:

- `redbiom fetch samples` without `--md5` emits nucleotide sequence feature IDs;
- `redbiom fetch samples --md5 True` emits a sequence-to-md5 map;
- fetching all 24,566 AGP samples directly was too slow;
- fetching all features contained in the context was fast enough and complete.

The successful strategy was:

```bash
redbiom fetch features-contained \
  --context Deblur_2021.09-Illumina-16S-V4-100nt-50b3a2 \
  > /Users/pgajer/current_projects/gut_microbiome/outputs/agp_silva_taxonomy/2026-04-26-redbiom-md5-sequence-map/agp_context_features_100nt.txt
```

I then computed md5(sequence) over all context feature sequences and matched
those hashes to the 286,255 md5 feature IDs in the local AGP QZA.

Result:

- local AGP md5 ASVs: 286,255;
- matched to redbiom 100 nt context sequences: 286,255;
- missing local md5 IDs: 0.

The complete md5-to-sequence mapping for local AGP ASVs is:

`/Users/pgajer/current_projects/gut_microbiome/outputs/agp_silva_taxonomy/2026-04-26-redbiom-md5-sequence-map/agp_context_features_local_md5_map.tsv`

The missing-ID audit file is:

`/Users/pgajer/current_projects/gut_microbiome/outputs/agp_silva_taxonomy/2026-04-26-redbiom-md5-sequence-map/agp_context_features_missing_local_md5.txt`

That file should be empty.

## Scripts Added Or Modified

All scripts below are in:

`/Users/pgajer/current_projects/gut_microbiome/scripts`

Added:

- `prepare_agp_redbiom_md5_taxonomy_inputs.py`
  - filters a pre-extracted local AGP md5 feature list against the redbiom
    md5-to-sequence map;
  - matches local md5 IDs to redbiom context sequences;
  - writes the `asv-map` used by `assign_agp_md5_taxonomy.R` plus a lightweight
    placeholder seqtab artifact.
- `assign_agp_md5_taxonomy.R`
  - assigns AGP ASV taxonomy with DADA2/SILVA;
  - mirrors the validation taxonomy assignment approach;
  - works row-wise on md5/sequence input rather than requiring a huge dense ASV
    table.
- `collapse_agp_qza_by_taxonomy.py`
  - collapses the local AGP QZA ASV count table to assigned SILVA taxon labels.
- `build_agp_silva_transfer_association_summary.py`
  - builds the association manifest;
  - audits overlap and nominal depth-specific assignment;
  - gathers top associations;
  - independently recomputes 2x2 Fisher label associations from assignments;
  - compares old and new frozen-AGP transfer outputs.

Modified:

- `run_agp_absorb_depthscan.R`
  - changed `compact_taxon_label()` so already-compact labels are preserved
    instead of converted to `Unassigned`;
  - this is needed because the collapsed AGP count table columns are already
    compact labels such as `Bacteroides dorei`.
- `run_external_case_control_validation_frozen_agp.R`
  - added `--carry-forward-terminal-depths true|false`, default false;
  - added terminal-depth carry-forward behavior;
  - added fallback to dominant mapped AGP taxon when no depth-1 root candidate
    has positive support but the sample has mapped AGP counts;
  - improved normalization for semicolon taxonomy strings and genus-only labels,
    which was needed for PRIME-style full taxonomy matrices such as HMP2 and
    PRJEB84421.

Parse verification for the modified R transfer script passed with `Rscript`.

## Taxonomy Assignment Outputs

Main taxonomy output root:

`/Users/pgajer/current_projects/gut_microbiome/outputs/agp_silva_taxonomy/2026-04-26-redbiom-md5-sequence-map`

Input and mapping files:

- `agp_local_qza_sample_ids.txt`
- `agp_local_qza_md5_feature_ids.txt`
- `agp_context_features_100nt.txt`
- `agp_context_features_local_md5_map.tsv`
- `agp_context_features_missing_local_md5.txt`
- `taxonomy_inputs/agp_md5_asv_sequences.tsv`
- `taxonomy_inputs/agp_md5_placeholder_seqtab.tsv.gz`
- `taxonomy_inputs/prepare_agp_redbiom_md5_taxonomy_inputs_summary.tsv`

Taxonomy assignment files:

- `taxonomy_assignment/asv_taxonomy.tsv`
- `taxonomy_assignment/taxon_label_summary.tsv`
- `taxonomy_assignment/taxonomy_summary.tsv`
- `taxonomy_cache/sequence_taxonomy_cache.tsv.gz`

Collapsed AGP count files:

- `agp_silva_collapsed_counts.tsv.gz`
- `agp_silva_collapsed_counts_summary.tsv`
- `agp_all_good_depthscan_metadata.tsv.gz`

Taxonomy summary:

| Metric | Value |
|---|---:|
| Input ASVs | 286,255 |
| Unique sequences | 286,255 |
| Unique SILVA taxon labels | 3,620 |
| Exact species labels | 1,881 |
| Genus-level labels | 89,798 |

Collapsed count summary:

| Metric | Value |
|---|---:|
| Samples | 24,566 |
| Input ASVs | 286,255 |
| Collapsed taxon labels | 3,620 |

The DADA2/SILVA references used were:

`/Users/pgajer/current_projects/gut_microbiome/data/reference/dada2_taxonomy/silva_nr99_v138.2_toSpecies_trainset.fa.gz`

`/Users/pgajer/current_projects/gut_microbiome/data/reference/dada2_taxonomy/silva_v138.2_assignSpecies.fa.gz`

## AGP Absorb Hierarchy Rerun

Main AGP absorb hierarchy output root:

`/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_analysis/runs/2026-04-26-agp-silva-local-qza-absorb-depth4`

The rerun used:

- count table: `agp_silva_collapsed_counts.tsv.gz`;
- metadata: `agp_all_good_depthscan_metadata.tsv.gz`;
- maximum depth: 4;
- `n0` schedule: `50,25,12,10`;
- run mode: exploratory;
- omnibus mode: standard;
- Fisher permutations: `1000`;
- minimum cases: `10`.

Representative command:

```bash
Rscript scripts/run_agp_absorb_depthscan.R \
  --output-dir /Users/pgajer/current_projects/gut_microbiome/outputs/dcst_analysis/runs/2026-04-26-agp-silva-local-qza-absorb-depth4 \
  --counts /Users/pgajer/current_projects/gut_microbiome/outputs/agp_silva_taxonomy/2026-04-26-redbiom-md5-sequence-map/agp_silva_collapsed_counts.tsv.gz \
  --metadata /Users/pgajer/current_projects/gut_microbiome/outputs/agp_silva_taxonomy/2026-04-26-redbiom-md5-sequence-map/agp_all_good_depthscan_metadata.tsv.gz \
  --linf-root /Users/pgajer/current_projects/linf \
  --agp-project PRJEB11419 \
  --n0-schedule 50,25,12,10 \
  --run-mode exploratory \
  --omnibus-mode standard \
  --fisher-b 1000 \
  --min-cases 10
```

Input and filtering summary:

- raw AGP matrix: 24,566 samples x 3,620 SILVA labels;
- samples kept after minimum library size: 24,421 / 24,566;
- features kept after prevalence filtering: 340 / 3,620;
- prevalence threshold: 5 percent, or 1,222 samples.

Depth summary:

| Depth | Labels | Samples | n0 | Largest label | Largest n |
|---:|---:|---:|---:|---|---:|
| 1 | 50 | 24,421 | 50 | `Bacteroides dorei` | 8,734 |
| 2 | 158 | 24,421 | 25 | `Bacteroides dorei__Faecalibacterium prausnitzii` | 2,320 |
| 3 | 402 | 24,421 | 12 | `Bacteroides dorei__Faecalibacterium prausnitzii__Agathobacter faecis` | 372 |
| 4 | 464 | 24,421 | 10 | `Segatella sp.__Faecalibacterium prausnitzii__Agathobacter faecis__Bacteroides dorei` | 170 |

Generated hierarchy assets:

- `agp_absorb_assignments.tsv.gz`
- `agp_absorb_label_summary_by_depth.tsv`
- `agp_absorb_object.rds`
- `agp_absorb_reassignment_summary.tsv`
- `baseline_manifest.tsv`
- `depth_summary.tsv`
- `largest_labels_by_depth.tsv`
- `taxon_reference.tsv`
- `omnibus_by_depth.tsv`
- `GO_NO_GO_depth3_depth4.md`
- `depth12_crosstab.tsv`
- `depth23_crosstab.tsv`
- `depth34_crosstab.tsv`
- `depth12_heatmap.png`
- `depth23_heatmap.png`
- `depth34_heatmap.png`
- `omnibus_q_heatmap.png`
- `omnibus_v_heatmap.png`
- `run_summary.md`

## Transfer Rerun

Main transfer output root:

`/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_validation/frozen_agp_silva_local_qza_2026-04-26`

The transfer rerun used the local-QZA SILVA AGP hierarchy:

`/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_analysis/runs/2026-04-26-agp-silva-local-qza-absorb-depth4/agp_absorb_assignments.tsv.gz`

and taxon reference:

`/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_analysis/runs/2026-04-26-agp-silva-local-qza-absorb-depth4/taxon_reference.tsv`

Transfer was run with:

- maximum depth: 4;
- terminal-depth carry-forward: enabled.

Validation comparisons rerun:

| Dataset | Comparisons |
|---|---|
| Halfvarson 2017 | `IBD_vs_Healthy`, `Crohn_vs_Healthy`, `UC_vs_Healthy` |
| Gevers 2014 | `Crohn_vs_Healthy` |
| HMP2 / IBDMDB | `IBD_vs_Healthy`, `Crohn_vs_Healthy`, `UC_vs_Healthy`, `Crohn_vs_UC` |
| PRJEB84421 | `Inflammatory_vs_Healthy`, `Crohn_vs_Healthy`, `OFG_vs_Healthy`, `Crohn_vs_OFG` |
| Jacobs 2023 IBS 250 bp | `IBS_vs_Healthy` |
| Jacobs 2023 IBS 150 bp sample-level | `IBS_vs_healthy` |
| Jacobs 2023 IBS 150 bp subject-level | `IBS_vs_healthy` |

No-sample-loss audit:

- comparisons checked: 15;
- depth rows checked: 60;
- initial AGP-overlap loss rows: 0;
- depth-loss rows: 0;
- comparisons with overlap loss: 0.

Audit files:

- `NO_SAMPLE_LOSS_AUDIT.tsv`
- `NO_SAMPLE_LOSS_AUDIT_summary.txt`

In the corrected transfer implementation, the overlap audit should be read in
two layers:

- root overlap retention: whether `samples_after_raw_qc = samples_after_agp_overlap`
- nominal depth-specific assignment: whether a sample received a valid AGP
  lineage-set assignment at depth `d` under the current nearest-lineage-set
  transfer rule

Under the finalized semantics, repeated lineage labels can be fully valid at
later nominal depths when those same repeated lineage sets are realized by AGP
at those depths. Accordingly, transfer reporting should be interpreted in terms
of valid nominal depth-specific assignment, not by imposing an extra
lineage-length-based notion of "effective depth."

## Per-Comparison Transfer Assets

For each comparison, the transfer output directory contains:

- `*_agp_taxon_mapping.csv`
- `*_agp_taxon_reference.csv`
- `*_dcst_assignments.csv`
- `*_depth1_results.csv`
- `*_depth2_results.csv`
- `*_depth3_results.csv`
- `*_depth4_results.csv`
- `*_depth1_sizes.png`
- `*_depth2_sizes.png`
- `*_depth3_sizes.png`
- `*_depth4_sizes.png`
- `*_mapping_summary.csv`
- `*_metadata_subset.csv`
- `*_run_info.csv`
- `*_validation_summary.md`

Dataset-specific directories under the transfer root:

- `gevers_2014/`
- `halfvarson_2017/`
- `hmp2/`
- `jacobs_2023_ibs_250bp/`
- `jacobs_2023_ibs_150bp/sample_level_depth1_4/`
- `jacobs_2023_ibs_150bp/subject_level_depth1_4/`
- `prjeb84421/`

There are 15 comparison assignment files, 60 depth-specific result CSVs, and
60 depth-size PNGs in the new transfer output family.

## Association Summary Outputs

Summary output root:

`/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_validation/frozen_agp_silva_local_qza_2026-04-26/summaries`

Generated files:

- `README.md`
- `association_manifest.tsv`
- `no_sample_loss_audit.tsv`
- `top_associations_by_dataset.tsv`
- `recomputed_label_associations.tsv`
- `old_vs_new_association_delta.tsv`

Rows, excluding header:

| File | Rows |
|---|---:|
| `association_manifest.tsv` | 60 |
| `no_sample_loss_audit.tsv` | 60 |
| `top_associations_by_dataset.tsv` | 468 |
| `recomputed_label_associations.tsv` | 2,889 |
| `old_vs_new_association_delta.tsv` | 4,969 |

Independent recomputation matched the R-generated transfer result tables on:

- `n_DCST`;
- `n_cases`;
- `mapped_samples`.

for all 2,889 label rows.

Quick signal count among top-10-per-depth rows:

- q <= 0.10: 58 rows;
- q <= 0.05: 46 rows;
- q <= 0.01: 24 rows.

## New Main Association Results

### Halfvarson 2017

Halfvarson remains the strongest positive validation family, but the label
story changes materially.

Depth-1 findings:

| Comparison | New best signal | OR | q | n |
|---|---|---:|---:|---:|
| IBD vs Healthy | `Segatella sp.` depletion | 0.140625 | 1.28e-08 | 533 |
| Crohn vs Healthy | `Segatella sp.` depletion | 0.1449 | 4.06e-07 | 254 |
| UC vs Healthy | `Segatella sp.` depletion | 0.1376 | 5.66e-08 | 333 |

Other depth-1 Halfvarson signals:

- `Bacteroides dorei` is enriched in IBD, Crohn, and UC comparisons.
- `Akkermansia muciniphila` is depleted in IBD, Crohn, and UC comparisons.
- `Prevotellaceae fam.` is depleted in Halfvarson IBD at q approximately 0.070.

Old-vs-new interpretation:

- old paper-facing best signals emphasized broad `Bacteroides`,
  `Faecalibacterium`, and `Akkermansia` labels;
- new local-QZA SILVA transfer sharpens the leading Halfvarson pattern to
  `Segatella sp.` depletion, `Bacteroides dorei` enrichment, and
  `Akkermansia muciniphila` depletion.

### HMP2 / IBDMDB

The old HMP2 positive transfer signal does not survive the local-QZA SILVA
transfer rerun.

Key shift:

- old best HMP2 IBD signal: depth-2
  `Bacteroides__Faecalibacterium`, q = 0.027, n = 164;
- new best HMP2 IBD signal: depth-1 `Alistipes sp.`, q = 0.742, n = 166.

Other HMP2 comparisons are also nonsignificant:

- Crohn vs Healthy: new best q approximately 0.771;
- UC vs Healthy: new best q approximately 1;
- Crohn vs UC: new best q approximately 1.

Paper implication: any current manuscript claim that HMP2 provides a corrected
direct-transfer replication signal should be removed or softened.

### Gevers 2014

Gevers remains negative under the new transfer.

Key shift:

- old best Gevers Crohn signal: depth-2
  `Bacteroides__Dialister`, q approximately 0.972, n = 308;
- new best Gevers Crohn signal: depth-4
  `Bacteroides fragilis__Faecalibacterium prausnitzii__Bacteroides ovatus__Bacteroides thetaiotaomicron`,
  q approximately 0.883, n = 396.

Paper implication: Gevers continues to act as a boundary-setting negative
validation cohort, now without sample-loss ambiguity.

### Jacobs 2023 IBS

The IBS results remain weak overall. The old 150 bp `Eukaryota` signal appears
to have been a taxonomy-transfer artifact.

Key shifts:

- Jacobs 250 bp IBS: old best q approximately 0.722, new best q approximately
  1;
- Jacobs 150 bp sample-level: old best `Eukaryota`, q approximately 1.11e-41;
  new best `Bacillota phy.`, q approximately 0.00692;
- Jacobs 150 bp subject-level: old best `Eukaryota`, q approximately 1.63e-29;
  new best `Bacillota phy.`, q approximately 0.0236.

Paper implication: the old `Eukaryota` result should not be retained as a
biological claim.

### PRJEB84421

PRJEB84421 gains some evidence for `Blautia sp.`-related contrasts, especially
in OFG comparisons.

Key new best signals:

- Inflammatory vs Healthy: `Blautia sp.`, q approximately 0.0588;
- OFG vs Healthy: `Blautia sp.`, q approximately 0.00191;
- Crohn vs OFG: `Blautia sp.`, OR approximately 0.1119, q approximately
  0.0354;
- Crohn vs Healthy: remains nonsignificant, new best q approximately 0.985.

Paper implication: PRJEB84421 may now be useful as an OFG-sensitive external
cohort rather than a simple Crohn-vs-healthy replication cohort.

## Paper Impact

This rerun changes the paper materially.

What is strengthened:

- the sample-loss/mapping-bias concern is resolved;
- all transfer analyses now retain every eligible sample at depths 1-4;
- taxonomy assignment is now harmonized between AGP discovery and validation
  datasets;
- Halfvarson remains strongly positive after this correction.

What changes:

- the main Halfvarson label story changes from broad genus/family-level labels
  to more specific local-QZA SILVA labels;
- HMP2 no longer supports a significant direct-transfer claim;
- the old Jacobs 150 bp `Eukaryota` result should be treated as an artifact;
- PRJEB84421 now has potentially interesting OFG-associated `Blautia sp.`
  signals;
- all transfer coverage figures and mapping-loss claims based on the old
  transfer run are obsolete.

Manuscript text that likely needs revision:

- the direct transfer validation Results section;
- the Methods description of the AGP source table and taxonomy assignment;
- Figure 4 or any figure that shows transfer coverage loss;
- Supplementary Table S11B transfer summary;
- Supplementary Table S13 transfer mapping balance;
- appendix transfer association tables if the local-QZA SILVA transfer becomes
  the canonical analysis.

Old manuscript numbers that should be replaced if this rerun is adopted:

- old AGP hierarchy: 40 / 153 / 479 / 744 labels;
- new AGP hierarchy: 50 / 158 / 402 / 464 labels;
- old AGP source description: PRIME-derived SILVA gut-project abundance table,
  34,711 AGP rows, 30,290 retained, 274 taxa;
- new AGP source description: local AGP QZA with 24,566 samples x 286,255 ASVs,
  redbiom md5-to-sequence recovery, DADA2/SILVA reassignment, collapse to 3,620
  labels, then hierarchy filtering to 24,421 samples x 340 labels.

## Current Manuscript Assets Affected

Existing paper assets in `linf` that are now stale relative to the local-QZA
SILVA transfer rerun:

- `/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex`
- `/Users/pgajer/current_projects/linf/papers/gut-dcst/assets/tables/TABLE_S11B_transfer_external_validation_summary.tsv`
- `/Users/pgajer/current_projects/linf/papers/gut-dcst/assets/tables/TABLE_S11B_transfer_external_validation_summary.tex`
- `/Users/pgajer/current_projects/linf/papers/gut-dcst/assets/tables/TABLE_S13_external_label_transfer_mapping_balance.tsv`
- `/Users/pgajer/current_projects/linf/papers/gut-dcst/assets/tables/TABLE_S13_external_label_transfer_mapping_balance.tex`
- `/Users/pgajer/current_projects/linf/papers/gut-dcst/assets/tables/TABLE_S14_S16_external_transfer_master.tsv`
- `/Users/pgajer/current_projects/linf/papers/gut-dcst/assets/tables/TABLE_S14_S16_external_transfer_association_tables.tex`
- `/Users/pgajer/current_projects/linf/papers/gut-dcst/assets/tables/external_validation_transfer_tables/`
- `/Users/pgajer/current_projects/linf/papers/gut-dcst/assets/tables/validation_dataset_appendix_tables/`

I have not overwritten those paper assets in this step. They remain the old
paper-facing state until a deliberate table/figure/manuscript migration is
performed.

## Recommended Next Steps

1. Decide whether the local-QZA SILVA transfer becomes the canonical transfer
   analysis for the manuscript.
2. If yes, regenerate manuscript-facing transfer summary tables from:
   `/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_validation/frozen_agp_silva_local_qza_2026-04-26/summaries`
3. Replace the mapping-balance table with an overlap/nominal-assignment audit table
   or remove that table entirely.
4. Rewrite the direct-transfer Results section around:
   - perfect transfer retention;
   - Halfvarson positive results;
   - HMP2 negative result after harmonization;
   - Gevers negative result;
   - PRJEB84421 OFG/Blautia result;
   - Jacobs artifact cleanup.
5. Update Methods to describe local AGP QZA, redbiom md5-to-sequence recovery,
   DADA2/SILVA taxonomy assignment, and collapsed AGP count generation.
6. Rebuild the manuscript PDF and visually check updated tables and figures.

## Minimal Starting Points For Future Codex Threads

To inspect the full new transfer state, start with:

`/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_validation/frozen_agp_silva_local_qza_2026-04-26/summaries/README.md`

To inspect the overlap/nominal-assignment audit, start with:

`/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_validation/frozen_agp_silva_local_qza_2026-04-26/NO_SAMPLE_LOSS_AUDIT.tsv`

To inspect the new AGP hierarchy, start with:

`/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_analysis/runs/2026-04-26-agp-silva-local-qza-absorb-depth4/depth_summary.tsv`

To inspect the taxonomy assignment, start with:

`/Users/pgajer/current_projects/gut_microbiome/outputs/agp_silva_taxonomy/2026-04-26-redbiom-md5-sequence-map/taxonomy_assignment/taxonomy_summary.tsv`

To inspect the code changes, start with:

`/Users/pgajer/current_projects/gut_microbiome/scripts/prepare_agp_redbiom_md5_taxonomy_inputs.py`

`/Users/pgajer/current_projects/gut_microbiome/scripts/assign_agp_md5_taxonomy.R`

`/Users/pgajer/current_projects/gut_microbiome/scripts/collapse_agp_qza_by_taxonomy.py`

`/Users/pgajer/current_projects/gut_microbiome/scripts/run_agp_absorb_depthscan.R`

`/Users/pgajer/current_projects/gut_microbiome/scripts/run_external_case_control_validation_frozen_agp.R`

`/Users/pgajer/current_projects/gut_microbiome/scripts/build_agp_silva_transfer_association_summary.py`
