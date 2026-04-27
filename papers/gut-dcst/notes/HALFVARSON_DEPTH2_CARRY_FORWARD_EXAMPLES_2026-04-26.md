# Halfvarson IBD examples with no true depth-2 AGP absorb mapping

Date: 2026-04-26

Criterion used here:

- `final_case_control_label = IBD_case`
- `depth1_frozen_agp_absorb` is non-empty
- `depth2_frozen_agp_absorb == depth1_frozen_agp_absorb`

Interpretation:

- these samples overlap the AGP taxon universe and receive a depth-1 frozen AGP label;
- they do **not** reach a more specific retained depth-2 AGP absorb descendant in the current local-QZA SILVA transfer;
- the depth-2 column is non-empty only because the transfer script carries the shallower label forward.

Summary for `halfvarson_2017__IBD_vs_Healthy`:

- qualifying IBD-case samples: `41`
- unique subjects represented: `35`
- phenotype split: `25` Crohn disease, `16` ulcerative colitis

Ten concrete examples:

| run_accession | phenotype | subject_id | timepoint | dominant_agp_taxon | depth1 | depth2 | agp_mapped_fraction | nonzero_features | top5 taxa |
|---|---|---|---:|---|---|---|---:|---:|---|
| `ERR1746351` | Ulcerative_colitis | `1629:21` | 8 | `Gemmiger sp.` | `Gemmiger sp.` | `Gemmiger sp.` | 0.9663 | 165 | `Gemmiger sp.=16496; Faecalibacterium prausnitzii=13549; Bifidobacterium breve=12407; Blautia sp.=11119; Ligilactobacillus sp.=9000` |
| `ERR1746365` | Crohn_disease | `1629:69` | 1 | `Intestinibacter bartlettii` | `Blautia sp.` | `Blautia sp.` | 0.0513 | 55 | `Escherichia-Shigella coli=921387; Intestinibacter bartlettii=31756; Blautia sp.=9490; Lachnospiraceae fam.=2224; Enterocloster sp.=2044` |
| `ERR1746366` | Crohn_disease | `1629:69` | 2 | `[Ruminococcus] gnavus group gnavus` | `Lachnospiraceae fam.` | `Lachnospiraceae fam.` | 0.4611 | 61 | `Escherichia-Shigella coli=123482; [Ruminococcus] gnavus group gnavus=36262; Salmonella enterica=26862; Blautia hansenii=16809; Lachnospiraceae fam.=15696` |
| `ERR1746382` | Crohn_disease | `1629:93` | 4 | `Gemmiger formicilis` | `Gemmiger formicilis` | `Gemmiger formicilis` | 0.9992 | 98 | `Gemmiger formicilis=10466; Gemmiger sp.=8872; Bacteroides dorei=5338; Acidaminococcus intestini=4013; Bacteroides ovatus=4007` |
| `ERR1746399` | Ulcerative_colitis | `1629:122` | 1 | `Ruminococcus bromii` | `Ruminococcus bromii` | `Ruminococcus bromii` | 0.9890 | 125 | `Ruminococcus bromii=36134; Blautia sp.=28374; Blautia massiliensis=16975; Christensenellaceae R-7 group sp.=13781; Erysipelotrichaceae UCG-003 bacterium=12770` |
| `ERR1746404` | Ulcerative_colitis | `1629:136` | 1 | `Bacteroides stercoris` | `Gemmiger sp.` | `Gemmiger sp.` | 0.9434 | 161 | `Bacteroides stercoris=41306; Gemmiger sp.=27705; Faecalibacterium prausnitzii=20132; Bacteroides intestinalis=20131; Anaerostipes hadrus=18933` |
| `ERR1746410` | Ulcerative_colitis | `1629:138` | 4 | `CAG-352 bromii` | `CAG-352 bromii` | `CAG-352 bromii` | 0.9936 | 181 | `CAG-352 bromii=99669; Bacteroides sp.=46783; Prevotellaceae NK3B31 group sp.=32990; Agathobacter faecis=20911; Faecalibacterium prausnitzii=18676` |
| `ERR1746429` | Ulcerative_colitis | `1629:152` | 4 | `Ruminococcus bromii` | `Ruminococcus bromii` | `Ruminococcus bromii` | 0.9231 | 161 | `Ruminococcus bromii=28439; Faecalibacterium prausnitzii=19264; Bacteroides dorei=18370; Rhodospirillales ord.=17736; Bacteroides ovatus=15674` |
| `ERR1746450` | Crohn_disease | `1629:191` | 4 | `Blautia sp.` | `Blautia sp.` | `Blautia sp.` | 0.9797 | 49 | `Blautia sp.=196970; Enterococcus faecium=9119; Lachnospiraceae fam.=8249; Blautia hansenii=7888; Escherichia-Shigella coli=4264` |
| `ERR1746455` | Crohn_disease | `1629:191` | 2 | `Blautia sp.` | `Blautia sp.` | `Blautia sp.` | 0.2561 | 43 | `Lactobacillus amylovorus=34106; Limosilactobacillus sp.=17313; Blautia sp.=7847; Bifidobacterium breve=5465; Fusobacterium nucleatum=5120` |
