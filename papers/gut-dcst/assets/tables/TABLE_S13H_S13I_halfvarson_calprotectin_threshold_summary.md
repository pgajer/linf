# Halfvarson calprotectin threshold sensitivity summary

## Sample universe
- Crohn runs with finite positive calprotectin in the Halfvarson validation matrix: 237
- Unique subjects: 50
- Repeated-sample subjects: 40
- Calprotectin threshold sensitivity split at 250 ug/g: 121 high / 116 lower visits

## Outputs
- TABLE_S13H_halfvarson_calprotectin_threshold250_sensitivity.tsv
- TABLE_S13H_halfvarson_calprotectin_threshold250_sensitivity.tex
- TABLE_S13I_halfvarson_calprotectin_threshold_sweep.tsv
- TABLE_S13I_halfvarson_calprotectin_threshold_sweep.tex
- FIGURE_S2_halfvarson_calprotectin_threshold_sweep.png

## Notes
- The manuscript-facing calprotectin follow-up is thresholded and sensitivity-oriented.
- The 250 ug/g table is retained as a clinically legible high-inflammatory-burden layer.
- Mechanistic context-aware F. prausnitzii modeling is intentionally left as a second-stage extension rather than part of this primary pipeline.

## Threshold-250 top signals

- Rebuilt depth 1: Bacteroides dorei | OR 1.715 [0.853, 3.500], q = 0.1391
- Rebuilt depth 2: Bacteroides dorei__Bacteroides uniformis | OR 2.270 [0.956, 5.739], q = 0.2014
- Rebuilt depth 3: Segatella sp.__Faecalibacterium prausnitzii__Escherichia-Shigella coli | OR 0.145 [0.015, 0.677], q = 0.0440
- Rebuilt depth 4: Segatella sp.__Faecalibacterium prausnitzii__Escherichia-Shigella coli | OR 0.145 [0.015, 0.677], q = 0.0538
- AGP transfer depth 1: Bacteroides fragilis | OR 0.242 [0.042, 0.954], q = 0.1989
- AGP transfer depth 2: Segatella sp.__Bacteroides dorei | OR 0.097 [0.002, 0.721], q = 0.1058
- AGP transfer depth 3: Segatella sp.__Bacteroides dorei__Faecalibacterium prausnitzii | OR 0.110 [0.002, 0.845], q = 0.1844
- AGP transfer depth 4: Segatella sp.__Bacteroides dorei__Faecalibacterium prausnitzii__Agathobacter faecis | OR 0.000 [0.000, 0.799], q = 0.1438

## Threshold-sweep top signals

- 100 ug/g | Rebuilt depth 1: Segatella sp. | OR 0.887 [0.419, 1.936], q = 0.7207
- 100 ug/g | Rebuilt depth 2: Bacteroides dorei__Bacteroides ovatus | OR 3.832 [1.098, 20.661], q = 0.1455
- 100 ug/g | Rebuilt depth 3: Bacteroides dorei__Bacteroides ovatus | OR 3.832 [1.098, 20.661], q = 0.2182
- 100 ug/g | Rebuilt depth 4: Bacteroides dorei__Bacteroides ovatus | OR 3.832 [1.098, 20.661], q = 0.2667
- 100 ug/g | AGP transfer depth 1: Lachnospiraceae fam. | OR 0.250 [0.050, 1.100], q = 0.2406
- 100 ug/g | AGP transfer depth 2: Lachnospiraceae fam.__Faecalibacterium prausnitzii | OR 0.222 [0.033, 1.188], q = 0.2813
- 100 ug/g | AGP transfer depth 3: Lachnospiraceae fam.__Faecalibacterium prausnitzii__Bacteroides dorei | OR 0.223 [0.033, 1.192], q = 0.3315
- 100 ug/g | AGP transfer depth 4: Lachnospiraceae fam.__Faecalibacterium prausnitzii__Bacteroides dorei | OR 0.228 [0.034, 1.219], q = 0.3208
- 150 ug/g | Rebuilt depth 1: Segatella sp. | OR 0.790 [0.387, 1.639], q = 0.4948
- 150 ug/g | Rebuilt depth 2: Bacteroides dorei__Bacteroides ovatus | OR 2.126 [0.782, 6.754], q = 0.7504
- 150 ug/g | Rebuilt depth 3: Segatella sp.__Faecalibacterium prausnitzii__Escherichia-Shigella coli | OR 0.371 [0.101, 1.280], q = 0.6145
- 150 ug/g | Rebuilt depth 4: Segatella sp.__Faecalibacterium prausnitzii__Escherichia-Shigella coli | OR 0.371 [0.101, 1.280], q = 0.7510
- 150 ug/g | AGP transfer depth 1: Lachnospiraceae fam. | OR 0.335 [0.067, 1.471], q = 0.4140
- 150 ug/g | AGP transfer depth 2: Bacteroides dorei__Bacteroides ovatus | OR 3.562 [0.770, 33.483], q = 0.4743
- 150 ug/g | AGP transfer depth 3: Lachnospiraceae fam.__Faecalibacterium prausnitzii__Bacteroides dorei | OR 0.295 [0.044, 1.572], q = 0.5483
- 150 ug/g | AGP transfer depth 4: Lachnospiraceae fam.__Faecalibacterium prausnitzii__Bacteroides dorei | OR 0.302 [0.045, 1.611], q = 0.6760
- 200 ug/g | Rebuilt depth 1: Segatella sp. | OR 0.542 [0.267, 1.088], q = 0.0698
- 200 ug/g | Rebuilt depth 2: Segatella sp.__Faecalibacterium prausnitzii | OR 0.542 [0.267, 1.088], q = 0.3056
- 200 ug/g | Rebuilt depth 3: Segatella sp.__Faecalibacterium prausnitzii__Escherichia-Shigella coli | OR 0.191 [0.033, 0.754], q = 0.0878
- 200 ug/g | Rebuilt depth 4: Segatella sp.__Faecalibacterium prausnitzii__Escherichia-Shigella coli | OR 0.191 [0.033, 0.754], q = 0.1073
- 200 ug/g | AGP transfer depth 1: Lachnospiraceae fam. | OR 0.180 [0.018, 0.935], q = 0.1263
- 200 ug/g | AGP transfer depth 2: Segatella sp.__Bacteroides dorei | OR 0.076 [0.002, 0.569], q = 0.0319
- 200 ug/g | AGP transfer depth 3: Segatella sp.__Bacteroides dorei__Faecalibacterium prausnitzii | OR 0.088 [0.002, 0.679], q = 0.1325
- 200 ug/g | AGP transfer depth 4: Segatella sp.__Bacteroides dorei__Faecalibacterium prausnitzii__Agathobacter faecis | OR 0.000 [0.000, 0.643], q = 0.0718
- 250 ug/g | Rebuilt depth 1: Bacteroides dorei | OR 1.715 [0.853, 3.500], q = 0.1391
- 250 ug/g | Rebuilt depth 2: Bacteroides dorei__Bacteroides uniformis | OR 2.270 [0.956, 5.739], q = 0.2014
- 250 ug/g | Rebuilt depth 3: Segatella sp.__Faecalibacterium prausnitzii__Escherichia-Shigella coli | OR 0.145 [0.015, 0.677], q = 0.0440
- 250 ug/g | Rebuilt depth 4: Segatella sp.__Faecalibacterium prausnitzii__Escherichia-Shigella coli | OR 0.145 [0.015, 0.677], q = 0.0538
- 250 ug/g | AGP transfer depth 1: Bacteroides fragilis | OR 0.242 [0.042, 0.954], q = 0.1989
- 250 ug/g | AGP transfer depth 2: Segatella sp.__Bacteroides dorei | OR 0.097 [0.002, 0.721], q = 0.1058
- 250 ug/g | AGP transfer depth 3: Segatella sp.__Bacteroides dorei__Faecalibacterium prausnitzii | OR 0.110 [0.002, 0.845], q = 0.1844
- 250 ug/g | AGP transfer depth 4: Segatella sp.__Bacteroides dorei__Faecalibacterium prausnitzii__Agathobacter faecis | OR 0.000 [0.000, 0.799], q = 0.1438
