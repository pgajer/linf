# Exploratory Halfvarson Crohn calprotectin context model

- Sample universe: 215 Crohn runs across 50 subjects
- Exact duplicated run accessions removed before fitting: 19
- Epsilon added before log10 transforms: 1e-06

## Actual feature map

- Faecalibacterium prausnitzii anchor: matched 1 feature(s): Faecalibacterium prausnitzii
- cross_feeding: matched 23 feature(s): Bifidobacterium breve; Bifidobacterium sp.; Bifidobacterium animalis; Bifidobacterium bifidum; Bifidobacterium longum; Bifidobacterium dentium; Bifidobacterium angulatum; Bifidobacterium pseudolongum; Bifidobacterium faecale; Bifidobacterium pseudocatenulatum; Blautia sp.; Blautia massiliensis; Blautia obeum; Blautia hansenii; Blautia caecimuris; Blautia hydrogenotrophica; Blautia glucerasea; Blautia wexlerae; Blautia coccoides; Roseburia sp.; Roseburia inulinivorans; Roseburia hominis; Roseburia intestinalis
- oxy_tolerators: matched 2 feature(s): Escherichia-Shigella coli; Enterobacteriaceae fam.
- bacteroides_bg: matched 21 feature(s): Bacteroides dorei; Bacteroides uniformis; Bacteroides ovatus; Bacteroides fragilis; Bacteroides thetaiotaomicron; Bacteroides stercoris; Bacteroides caccae; Bacteroides cellulosilyticus; Bacteroides sartorii; Bacteroides eggerthii; Bacteroides vulgatus; Bacteroides sp.; Bacteroides salyersiae; Bacteroides intestinalis; Bacteroides coprophilus; Bacteroides coprocola; Bacteroides plebeius; Bacteroides clarus; Bacteroides mediterraneensis; Bacteroides ndongoniae; [Bacteroides] pectinophilus group pectinophilus
- mucin_context: matched 2 feature(s): Akkermansia muciniphila; Akkermansia sp.

## Interaction-term summary

- log10_Fp:cross_feeding: estimate -0.005 [-0.021, 0.011], approx. p = 0.5231
- log10_Fp:oxy_tolerators: estimate 0.011 [-0.003, 0.024], approx. p = 0.1374
- log10_Fp:bacteroides_bg: estimate -0.015 [-0.037, 0.007], approx. p = 0.1907
- log10_Fp:mucin_context: estimate 0.002 [-0.015, 0.019], approx. p = 0.8121

## Simple slopes

- cross_feeding | low (Q1): slope 0.006 [-0.031, 0.042]
- cross_feeding | median: slope 0.003 [-0.035, 0.041]
- cross_feeding | high (Q3): slope 0.001 [-0.039, 0.041]
- oxy_tolerators | low (Q1): slope -0.009 [-0.047, 0.030]
- oxy_tolerators | median: slope 0.003 [-0.035, 0.041]
- oxy_tolerators | high (Q3): slope 0.015 [-0.029, 0.058]
- bacteroides_bg | low (Q1): slope 0.008 [-0.029, 0.044]
- bacteroides_bg | median: slope 0.003 [-0.035, 0.041]
- bacteroides_bg | high (Q3): slope 0.001 [-0.039, 0.040]
- mucin_context | low (Q1): slope 0.000 [-0.043, 0.043]
- mucin_context | median: slope 0.003 [-0.035, 0.041]
- mucin_context | high (Q3): slope 0.006 [-0.040, 0.051]
