# dCSTs for Vaginal Microbiome Data

## Introduction

Community state types (CSTs) are the standard way to describe the
structure of vaginal microbial communities. They were originally defined
by hierarchical clustering of 16S rRNA data (Ravel et al. 2011) and
later standardized by the VALENCIA nearest-centroid classifier (France
et al. 2020).

**Dominant Community State Types (dCSTs)** offer an alternative that
requires no clustering and no training data. A dCST is defined by
$`L^\infty`$-normalization: divide each sample’s abundance vector by its
maximum value, then group samples by which taxon achieves that maximum.
This is a purely geometric operation on compositional data — it
identifies the dominant species in each sample and groups samples
accordingly.

This vignette demonstrates the dCST workflow on a real vaginal 16S
dataset and shows that dCSTs recover the established CST structure with
high concordance.

**Reference:** Gajer, P. & Ravel, J. (2025). A New Approach to
Compositional Data Analysis using $`L^\infty`$-normalization with
Applications to Vaginal Microbiome. arXiv:2503.21543.

## The Valencia 2k dataset

The `valencia2k` dataset is a stratified subsample of 2,000 vaginal
samples from the Valencia CST-classifier training set (France et
al. 2020). It contains species-level relative abundances for 178 taxa
along with the original Valencia CST and sub-CST assignments.

``` r

data(valencia2k)

## Compositional matrix: 2000 samples x 178 taxa
dim(valencia2k$rel)
#> [1] 2000  178

## Valencia CST assignments
table(valencia2k$cst$Val_CST)
#> 
#>    I   II  III IV-A IV-B IV-C    V 
#>  557   69  582  141  454  117   80

## Top 10 most abundant taxa (by mean relative abundance)
means <- sort(colMeans(valencia2k$rel), decreasing = TRUE)
head(round(means, 4), 10)
#>     Lactobacillus_iners Lactobacillus_crispatus   Gardnerella_vaginalis 
#>                  0.2609                  0.2539                  0.1203 
#>  Lactobacillus_jensenii       Atopobium_vaginae                   BVAB1 
#>                  0.0418                  0.0370                  0.0297 
#>   Lactobacillus_gasseri   Sneathia_sanguinegens         g_Streptococcus 
#>                  0.0283                  0.0198                  0.0194 
#>           g_Megasphaera 
#>                  0.0158
```

The dominant taxa are familiar vaginal species: *Lactobacillus iners*,
*L. crispatus*, *Gardnerella vaginalis*, and others. This is exactly the
kind of data where dominant-feature assignments are most intuitive: most
samples are dominated by a single species, and the identity of that
species is biologically meaningful.

## Reconstructing counts

The bundled matrix is compositional (rows sum to 1). For workflows that
require count-like input (e.g.,
[`filter.asv()`](https://pgajer.github.io/linf/reference/filter.asv.md)),
we can reconstruct approximate counts using the per-sample read depths:

``` r

count_mat <- sweep(valencia2k$rel, 1, valencia2k$reads, "*")
count_mat <- round(count_mat)
storage.mode(count_mat) <- "integer"

## Check: library sizes
summary(rowSums(count_mat))
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>    3026   22281   47958   57111   74849  409411
```

Since the Valencia 13k data was already filtered (min 3,000 reads, taxa
present in $`\ge`$ 1% of samples), we skip
[`filter.asv()`](https://pgajer.github.io/linf/reference/filter.asv.md)
here and work directly with the compositional matrix.

## $`L^\infty`$ normalization

$`L^\infty`$-normalization divides each row by its maximum value. After
normalization, every sample has at least one coordinate equal to 1, and
all other coordinates lie in $`[0, 1]`$.

``` r

Z <- normalize.linf(valencia2k$rel)

## Every row max should be 1
summary(apply(Z, 1, max))
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>       1       1       1       1       1       1
```

## Dominance sample sets: dominant species assignment

The dominant species is the column that achieves the maximum in each
row. Samples with the same dominant species form a depth-1 dominance
sample set. The assignment is deterministic and sample-independent: it
depends only on the sample’s own abundance profile, never on other
samples.

``` r

dominant.features <- linf.dominant.features(Z)

## How many distinct dominant taxa?
cat("Distinct dominant taxa:", length(dominant.features$observed.levels), "\n")
#> Distinct dominant taxa: 35

## Dominance sample-set size distribution (top 15)
dominance_tab <- sort(table(dominant.features$label), decreasing = TRUE)
head(dominance_tab, 15)
#> 
#>     Lactobacillus_iners Lactobacillus_crispatus   Gardnerella_vaginalis 
#>                     606                     564                     360 
#>                   BVAB1  Lactobacillus_jensenii   Lactobacillus_gasseri 
#>                      91                      67                      62 
#>       Atopobium_vaginae   Sneathia_sanguinegens         g_Streptococcus 
#>                      51                      40                      36 
#>       g_Bifidobacterium          g_Anaerococcus          g_Enterococcus 
#>                      26                      14                      13 
#>     g_Corynebacterium_1                BVAB_TM7        g_Staphylococcus 
#>                      11                       8                       7
```

We see the expected vaginal community structure: a few well-supported
dominance sample sets (the *Lactobacillus* species, *Gardnerella*,
BVAB1) and a long tail of rarer dominant taxa.

## Depth-1 dCSTs

In practice, many provisional dominance sample sets contain very few
samples. **Dominant Community State Types (dCSTs)** retain only states
above a minimum support threshold ($`n_0`$) and handle samples from
low-support sets according to the selected policy.

The package supports two policies for handling low-support provisional
sets:

- **“pure”**: named dCSTs stay pure, and samples from low-support sets
  are labeled `RARE_DOMINANT`
- **“absorb”**: samples from low-support sets are reassigned to a
  retained state by restricted argmax

``` r

dcst1 <- linf.csts(Z, n0 = 30, low.freq.policy = "pure")

## Retained depth-1 states
dcst1$retained.feature.labels
#> [1] "Lactobacillus_iners"     "Lactobacillus_crispatus"
#> [3] "Gardnerella_vaginalis"   "BVAB1"                  
#> [5] "Lactobacillus_jensenii"  "Lactobacillus_gasseri"  
#> [7] "Atopobium_vaginae"       "Sneathia_sanguinegens"  
#> [9] "g_Streptococcus"

## dCST assignment (pure policy)
dcst1_tab <- sort(table(dcst1$lineage.label), decreasing = TRUE)
dcst1_tab
#> 
#>     Lactobacillus_iners Lactobacillus_crispatus   Gardnerella_vaginalis 
#>                     606                     564                     360 
#>           RARE_DOMINANT                   BVAB1  Lactobacillus_jensenii 
#>                     123                      91                      67 
#>   Lactobacillus_gasseri       Atopobium_vaginae   Sneathia_sanguinegens 
#>                      62                      51                      40 
#>         g_Streptococcus 
#>                      36
```

### Concordance with Valencia CSTs

The key test: do dCSTs correspond to the established Valencia CSTs?
Let’s build a cross-tabulation:

``` r

## Build concordance table: rows = dCSTs, columns = Valencia CSTs
dcst_labels <- dcst1$lineage.label.absorb  # use absorb view for clean comparison
val_cst     <- valencia2k$cst$Val_CST

concordance <- table(dCST = dcst_labels, Valencia = val_cst)

## Show as percentage: what fraction of each dCST falls in each Valencia CST?
pct <- round(100 * prop.table(concordance, margin = 1), 1)
pct
#>                          Valencia
#> dCST                          I    II   III  IV-A  IV-B  IV-C     V
#>   Atopobium_vaginae         0.0   0.0   0.0   4.9  88.5   3.3   3.3
#>   BVAB1                     0.0   0.0   0.0 100.0   0.0   0.0   0.0
#>   g_Streptococcus           0.0   0.0   0.0   0.0   0.0 100.0   0.0
#>   Gardnerella_vaginalis     0.0   1.3   0.3   8.2  86.9   3.3   0.0
#>   Lactobacillus_crispatus  95.9   0.0   1.0   0.2   0.0   2.9   0.0
#>   Lactobacillus_gasseri     0.0  98.4   0.0   0.0   0.0   1.6   0.0
#>   Lactobacillus_iners       0.0   0.2  91.4   0.6   4.0   2.1   1.8
#>   Lactobacillus_jensenii    0.0   0.0   0.0   0.0   0.0   4.3  95.7
#>   Sneathia_sanguinegens     0.0   0.0   2.3  13.6  81.8   2.3   0.0
```

The concordance is striking. Each dCST maps predominantly to a single
Valencia CST:

- *Lactobacillus crispatus* dCST $`\to`$ CST I (dominated by *L.
  crispatus*)
- *Lactobacillus iners* dCST $`\to`$ CST III (dominated by *L. iners*)
- *Gardnerella vaginalis* dCST $`\to`$ CST IV-B
- BVAB1 dCST $`\to`$ CST IV-A
- *Lactobacillus gasseri* dCST $`\to`$ CST II
- *Lactobacillus jensenii* dCST $`\to`$ CST V

``` r

## Reverse view: what fraction of each Valencia CST falls in each dCST?
pct_rev <- round(100 * prop.table(concordance, margin = 2), 1)
pct_rev
#>                          Valencia
#> dCST                          I    II   III  IV-A  IV-B  IV-C     V
#>   Atopobium_vaginae         0.0   0.0   0.0   2.1  11.9   1.7   2.5
#>   BVAB1                     0.0   0.0   0.0  67.4   0.0   0.0   0.0
#>   g_Streptococcus           0.0   0.0   0.0   0.0   0.0  57.3   0.0
#>   Gardnerella_vaginalis     0.0   7.2   0.2  22.7  74.7  11.1   0.0
#>   Lactobacillus_crispatus 100.0   0.0   1.0   0.7   0.0  14.5   0.0
#>   Lactobacillus_gasseri     0.0  91.3   0.0   0.0   0.0   0.9   0.0
#>   Lactobacillus_iners       0.0   1.4  98.6   2.8   5.5  11.1  13.8
#>   Lactobacillus_jensenii    0.0   0.0   0.0   0.0   0.0   2.6  83.8
#>   Sneathia_sanguinegens     0.0   0.0   0.2   4.3   7.9   0.9   0.0
```

This two-way concordance confirms that dCSTs are an alternative
construct of CSTs: they recover the same groupings without any
clustering algorithm or training data.

## Depth-2 refinement

Large dCSTs can be refined by dropping the dominant species and
repeating the procedure on the remaining taxa. This produces a two-level
hierarchy that reveals sub-community structure.

``` r

dcst2 <- refine.linf.csts(Z, dcst1, n0 = 15, refinement.factor = 2,
                           low.freq.policy = "pure", verbose = FALSE)

## Depth-2 dCST table
dcst2_tab <- sort(table(dcst2$lineage.label), decreasing = TRUE)
dcst2_tab
#> 
#>             Lactobacillus_iners__Gardnerella_vaginalis 
#>                                                    204 
#>        Lactobacillus_crispatus__Lactobacillus_jensenii 
#>                                                    151 
#>           Lactobacillus_crispatus__Lactobacillus_iners 
#>                                                    131 
#>                                          RARE_DOMINANT 
#>                                                    123 
#>                     Lactobacillus_iners__RARE_DOMINANT 
#>                                                    105 
#>            Lactobacillus_iners__Lactobacillus_jensenii 
#>                                                    102 
#>           Lactobacillus_iners__Lactobacillus_crispatus 
#>                                                    101 
#>                   Gardnerella_vaginalis__RARE_DOMINANT 
#>                                                     95 
#>               Gardnerella_vaginalis__Atopobium_vaginae 
#>                                                     94 
#>               Lactobacillus_crispatus__g_Lactobacillus 
#>                                                     88 
#>             Gardnerella_vaginalis__Lactobacillus_iners 
#>                                                     85 
#>                           BVAB1__Gardnerella_vaginalis 
#>                                                     63 
#>                 Lactobacillus_crispatus__RARE_DOMINANT 
#>                                                     63 
#>                                  Lactobacillus_gasseri 
#>                                                     62 
#>           Lactobacillus_crispatus__g_Corynebacterium_1 
#>                                                     42 
#>               Atopobium_vaginae__Gardnerella_vaginalis 
#>                                                     37 
#>                                        g_Streptococcus 
#>                                                     36 
#>           Gardnerella_vaginalis__Sneathia_sanguinegens 
#>                                                     33 
#>                           Gardnerella_vaginalis__BVAB1 
#>                                                     29 
#>                                   BVAB1__RARE_DOMINANT 
#>                                                     28 
#>               Lactobacillus_iners__g_Corynebacterium_1 
#>                                                     27 
#>                   Lactobacillus_iners__g_Streptococcus 
#>                                                     26 
#>                Gardnerella_vaginalis__g_Fastidiosipila 
#>                                                     24 
#>            Lactobacillus_jensenii__Lactobacillus_iners 
#>                                                     24 
#>                  Lactobacillus_jensenii__RARE_DOMINANT 
#>                                                     24 
#>                   Sneathia_sanguinegens__RARE_DOMINANT 
#>                                                     24 
#>                      Lactobacillus_iners__g_Finegoldia 
#>                                                     22 
#>         Lactobacillus_crispatus__Gardnerella_vaginalis 
#>                                                     21 
#>         Lactobacillus_crispatus__Lactobacillus_reuteri 
#>                                                     21 
#>             Lactobacillus_iners__Lactobacillus_gasseri 
#>                                                     19 
#>          Lactobacillus_jensenii__Gardnerella_vaginalis 
#>                                                     19 
#>         Lactobacillus_crispatus__Lactobacillus_gasseri 
#>                                                     17 
#>           Sneathia_sanguinegens__Gardnerella_vaginalis 
#>                                                     16 
#> Lactobacillus_crispatus__g_Clostridium_sensu_stricto_1 
#>                                                     15 
#>                  Lactobacillus_crispatus__g_Finegoldia 
#>                                                     15 
#>                       Atopobium_vaginae__RARE_DOMINANT 
#>                                                     14
```

### Concordance with Valencia sub-CSTs

``` r

dcst2_labels <- dcst2$lineage.label.absorb
val_subcst   <- valencia2k$cst$Val_subCST

concordance2 <- table(dCST = dcst2_labels, Valencia = val_subcst)

## Percentage by dCST (row-wise)
pct2 <- round(100 * prop.table(concordance2, margin = 1), 1)

## Show only dCSTs with at least 20 samples for readability
dcst2_sizes <- rowSums(concordance2)
pct2[dcst2_sizes >= 20, , drop = FALSE]
#>                                                  Valencia
#> dCST                                                I-A   I-B    II III-A III-B
#>   Atopobium_vaginae__Gardnerella_vaginalis          0.0   0.0   0.0   0.0   0.0
#>   BVAB1__Gardnerella_vaginalis                      0.0   0.0   0.0   0.0   0.0
#>   g_Streptococcus                                   0.0   0.0   0.0   0.0   0.0
#>   Gardnerella_vaginalis                             0.0   0.0   0.0   0.0   0.0
#>   Gardnerella_vaginalis__Atopobium_vaginae          0.0   0.0   0.7   0.0   0.0
#>   Gardnerella_vaginalis__BVAB1                      0.0   0.0   0.0   0.0   0.0
#>   Gardnerella_vaginalis__g_Fastidiosipila           0.0   0.0   0.0   0.0   0.0
#>   Gardnerella_vaginalis__Lactobacillus_iners        0.0   0.0   3.3   0.0   0.8
#>   Gardnerella_vaginalis__Sneathia_sanguinegens      0.0   0.0   0.0   0.0   0.0
#>   Lactobacillus_crispatus__g_Corynebacterium_1     58.0  42.0   0.0   0.0   0.0
#>   Lactobacillus_crispatus__g_Finegoldia            60.0  40.0   0.0   0.0   0.0
#>   Lactobacillus_crispatus__g_Lactobacillus         90.2   9.8   0.0   0.0   0.0
#>   Lactobacillus_crispatus__Gardnerella_vaginalis   41.7  54.2   0.0   0.0   0.0
#>   Lactobacillus_crispatus__Lactobacillus_iners     44.4  51.1   0.0   0.0   4.5
#>   Lactobacillus_crispatus__Lactobacillus_jensenii  82.0  18.0   0.0   0.0   0.0
#>   Lactobacillus_crispatus__Lactobacillus_reuteri   74.2  25.8   0.0   0.0   0.0
#>   Lactobacillus_gasseri                             0.0   0.0  98.4   0.0   0.0
#>   Lactobacillus_iners                               0.0   0.0   0.0   4.5  18.2
#>   Lactobacillus_iners__g_Corynebacterium_1          0.0   0.0   0.0  91.3   8.7
#>   Lactobacillus_iners__g_Finegoldia                 0.0   0.0   0.0  81.6  18.4
#>   Lactobacillus_iners__g_Streptococcus              0.0   0.0   0.0  67.7  29.0
#>   Lactobacillus_iners__Gardnerella_vaginalis        0.0   0.0   0.0  39.3  50.8
#>   Lactobacillus_iners__Lactobacillus_crispatus      0.0   0.0   0.0  48.2  51.8
#>   Lactobacillus_iners__Lactobacillus_gasseri        0.0   0.0   5.0  80.0  15.0
#>   Lactobacillus_iners__Lactobacillus_jensenii       0.0   0.0   0.0  77.6  13.1
#>   Lactobacillus_jensenii__Gardnerella_vaginalis     0.0   0.0   0.0   0.0   0.0
#>   Lactobacillus_jensenii__Lactobacillus_iners       0.0   0.0   0.0   0.0   0.0
#>   Sneathia_sanguinegens__Gardnerella_vaginalis      0.0   0.0   0.0   0.0   2.5
#>                                                  Valencia
#> dCST                                               IV-A  IV-B IV-C0 IV-C1 IV-C2
#>   Atopobium_vaginae__Gardnerella_vaginalis          3.9  92.2   0.0   0.0   0.0
#>   BVAB1__Gardnerella_vaginalis                    100.0   0.0   0.0   0.0   0.0
#>   g_Streptococcus                                   0.0   0.0   9.0  58.2   3.0
#>   Gardnerella_vaginalis                            10.0  46.7  23.3   0.0   0.0
#>   Gardnerella_vaginalis__Atopobium_vaginae          0.0  99.3   0.0   0.0   0.0
#>   Gardnerella_vaginalis__BVAB1                     75.8  24.2   0.0   0.0   0.0
#>   Gardnerella_vaginalis__g_Fastidiosipila           3.4  96.6   0.0   0.0   0.0
#>   Gardnerella_vaginalis__Lactobacillus_iners        1.7  94.2   0.0   0.0   0.0
#>   Gardnerella_vaginalis__Sneathia_sanguinegens      2.4  97.6   0.0   0.0   0.0
#>   Lactobacillus_crispatus__g_Corynebacterium_1      0.0   0.0   0.0   0.0   0.0
#>   Lactobacillus_crispatus__g_Finegoldia             0.0   0.0   0.0   0.0   0.0
#>   Lactobacillus_crispatus__g_Lactobacillus          0.0   0.0   0.0   0.0   0.0
#>   Lactobacillus_crispatus__Gardnerella_vaginalis    4.2   0.0   0.0   0.0   0.0
#>   Lactobacillus_crispatus__Lactobacillus_iners      0.0   0.0   0.0   0.0   0.0
#>   Lactobacillus_crispatus__Lactobacillus_jensenii   0.0   0.0   0.0   0.0   0.0
#>   Lactobacillus_crispatus__Lactobacillus_reuteri    0.0   0.0   0.0   0.0   0.0
#>   Lactobacillus_gasseri                             0.0   0.0   0.0   0.0   1.6
#>   Lactobacillus_iners                               4.5  13.6  27.3   0.0   9.1
#>   Lactobacillus_iners__g_Corynebacterium_1          0.0   0.0   0.0   0.0   0.0
#>   Lactobacillus_iners__g_Finegoldia                 0.0   0.0   0.0   0.0   0.0
#>   Lactobacillus_iners__g_Streptococcus              0.0   0.0   0.0   0.0   0.0
#>   Lactobacillus_iners__Gardnerella_vaginalis        1.2   8.7   0.0   0.0   0.0
#>   Lactobacillus_iners__Lactobacillus_crispatus      0.0   0.0   0.0   0.0   0.0
#>   Lactobacillus_iners__Lactobacillus_gasseri        0.0   0.0   0.0   0.0   0.0
#>   Lactobacillus_iners__Lactobacillus_jensenii       0.0   0.0   0.0   0.0   0.0
#>   Lactobacillus_jensenii__Gardnerella_vaginalis     0.0   0.0   0.0   0.0   0.0
#>   Lactobacillus_jensenii__Lactobacillus_iners       0.0   0.0   0.0   0.0   0.0
#>   Sneathia_sanguinegens__Gardnerella_vaginalis     15.0  82.5   0.0   0.0   0.0
#>                                                  Valencia
#> dCST                                              IV-C3 IV-C4     V
#>   Atopobium_vaginae__Gardnerella_vaginalis          0.0   0.0   3.9
#>   BVAB1__Gardnerella_vaginalis                      0.0   0.0   0.0
#>   g_Streptococcus                                  28.4   1.5   0.0
#>   Gardnerella_vaginalis                             3.3  16.7   0.0
#>   Gardnerella_vaginalis__Atopobium_vaginae          0.0   0.0   0.0
#>   Gardnerella_vaginalis__BVAB1                      0.0   0.0   0.0
#>   Gardnerella_vaginalis__g_Fastidiosipila           0.0   0.0   0.0
#>   Gardnerella_vaginalis__Lactobacillus_iners        0.0   0.0   0.0
#>   Gardnerella_vaginalis__Sneathia_sanguinegens      0.0   0.0   0.0
#>   Lactobacillus_crispatus__g_Corynebacterium_1      0.0   0.0   0.0
#>   Lactobacillus_crispatus__g_Finegoldia             0.0   0.0   0.0
#>   Lactobacillus_crispatus__g_Lactobacillus          0.0   0.0   0.0
#>   Lactobacillus_crispatus__Gardnerella_vaginalis    0.0   0.0   0.0
#>   Lactobacillus_crispatus__Lactobacillus_iners      0.0   0.0   0.0
#>   Lactobacillus_crispatus__Lactobacillus_jensenii   0.0   0.0   0.0
#>   Lactobacillus_crispatus__Lactobacillus_reuteri    0.0   0.0   0.0
#>   Lactobacillus_gasseri                             0.0   0.0   0.0
#>   Lactobacillus_iners                              18.2   4.5   0.0
#>   Lactobacillus_iners__g_Corynebacterium_1          0.0   0.0   0.0
#>   Lactobacillus_iners__g_Finegoldia                 0.0   0.0   0.0
#>   Lactobacillus_iners__g_Streptococcus              0.0   0.0   3.2
#>   Lactobacillus_iners__Gardnerella_vaginalis        0.0   0.0   0.0
#>   Lactobacillus_iners__Lactobacillus_crispatus      0.0   0.0   0.0
#>   Lactobacillus_iners__Lactobacillus_gasseri        0.0   0.0   0.0
#>   Lactobacillus_iners__Lactobacillus_jensenii       0.0   0.0   9.3
#>   Lactobacillus_jensenii__Gardnerella_vaginalis     0.0   0.0 100.0
#>   Lactobacillus_jensenii__Lactobacillus_iners       0.0   0.0 100.0
#>   Sneathia_sanguinegens__Gardnerella_vaginalis      0.0   0.0   0.0
```

The depth-2 refinement recovers sub-CST structure. For example, the
*Lactobacillus crispatus* dCST splits based on the second-most-abundant
taxon, and these sub-dCSTs align with Valencia sub-CSTs I-A and I-B.

## The full pipeline

For convenience,
[`linf.dcst.landmark.pipeline()`](https://pgajer.github.io/linf/reference/linf.dcst.landmark.pipeline.md)
runs the entire workflow in one call: normalization, depth-1 dCSTs,
depth-2 refinement, and landmark computation.

``` r

out <- linf.dcst.landmark.pipeline(
  count_mat,
  feature.ids    = colnames(count_mat),
  feature.labels = colnames(count_mat),
  n0.depth1      = 30,
  n0.depth2      = 15,
  refinement.factor = 2,
  low.freq.policy   = "pure",
  landmark.view     = "absorb",
  verbose           = FALSE
)

names(out)
#> [1] "linf.rel"         "dcst.depth1"      "dcst.depth2"      "landmarks.depth1"
#> [5] "landmarks.depth2" "params"
```

## Landmark points

Each dCST dominance-lineage has characteristic landmark points: the most
extreme sample (highest value of the dominant taxon), the least extreme
(lowest value while still dominant), and mean/median representatives.
These help characterize what a “typical” or “extreme” member of each
community type looks like.

``` r

## Depth-1 landmarks
lm1 <- out$landmarks.depth1$landmarks

## Show landmarks for a few key dCST dominance-lineages
key_lineages <- c("Lactobacillus_crispatus", "Lactobacillus_iners",
                  "Gardnerella_vaginalis")
key_lm <- lm1[lm1$lineage.id %in% key_lineages, ]

## For each dCST, show the endpoint.max landmark's target value
## (how dominant is the dominant species in the most extreme sample?)
ep_max <- key_lm[key_lm$landmark.type == "endpoint.max",
                 c("lineage.label", "point.name", "target.value")]
ep_max
#>             lineage.label point.name target.value
#> 1 Lactobacillus_crispatus         s1            1
#> 5   Gardnerella_vaginalis         s6            1
#> 9     Lactobacillus_iners         s7            1

## And the endpoint.min (least dominant while still assigned to this lineage)
ep_min <- key_lm[key_lm$landmark.type == "endpoint.min",
                 c("lineage.label", "point.name", "target.value")]
ep_min
#>              lineage.label point.name target.value
#> 2  Lactobacillus_crispatus       s658   0.00504737
#> 6    Gardnerella_vaginalis       s577   0.01120848
#> 10     Lactobacillus_iners       s312   0.00000000
```

The endpoint.max landmarks show samples where the dominant species
comprises nearly the entire community, while endpoint.min landmarks show
the boundary cases — samples where dominance is marginal.

## Summary

This vignette demonstrated the dCST workflow on real vaginal 16S data:

1.  **$`L^\infty`$-normalization** maps compositional data to the
    $`L^\infty`$-simplex
2.  **Dominant-feature assignment** assigns each sample to its dominant
    taxon
3.  **Truncated dCSTs** retain well-supported states and handle rare
    dominance patterns
4.  **Depth-2 refinement** reveals sub-community structure
5.  **Landmark points** characterize extreme and typical community
    members

The concordance with established Valencia CSTs is high, validating dCSTs
as a principled, deterministic alternative to clustering-based community
typing. Key advantages:

- **No clustering algorithm:** assignment depends only on the sample
  itself
- **No training data:** new samples can be assigned without a reference
  database
- **Interpretable:** dCST labels directly name the dominant species
- **Hierarchical:** refinement naturally produces sub-types

## References

- France, M. T., Ma, B., Gajer, P., et al. (2020). VALENCIA: a nearest
  centroid classification method for vaginal microbial communities based
  on composition. *Microbiome*, 8(1), 166.
- Gajer, P. & Ravel, J. (2025). A New Approach to Compositional Data
  Analysis using $`L^\infty`$-normalization with Applications to Vaginal
  Microbiome. arXiv:2503.21543.
- Ravel, J., Gajer, P., Abdo, Z., et al. (2011). Vaginal microbiome of
  reproductive-age women. *PNAS*, 108(Suppl 1), 4680–4687.
