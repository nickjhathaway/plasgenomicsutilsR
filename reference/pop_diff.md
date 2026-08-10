# Population differentiation between metadata groups, per SNP

Computes a per-SNP differentiation statistic for every pair of `group`
levels, from alt-allele dosages. Run it on the **full, unpruned**
genotypes – LD-pruning removes the differentiating SNPs you want here.
Feed the result to
[`pop_diff_matrix()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff_matrix.md)
/
[`plot_diff_heatmap()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_diff_heatmap.md)
for a group summary, or
[`top_differentiating_snps()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/top_differentiating_snps.md)
for markers.

## Usage

``` r
pop_diff(
  x,
  group = NULL,
  statistic = c("jost_d", "gst_hedrick", "fst"),
  meta = NULL,
  clamp = TRUE,
  genotype = NULL
)
```

## Arguments

- x:

  A
  [PopStructure](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PopStructure.md)
  (uses its full genotype matrix for the active samples) or a genotype
  matrix (samples x SNPs, 0/1/2, `NA`).

- group:

  Metadata column (for a `PopStructure`) or a per-sample grouping vector
  (for a matrix).

- statistic:

  `"jost_d"` (Jost's D, default), `"gst_hedrick"` (Hedrick's
  standardized G'st), or `"fst"` (Hudson's Fst); all in \[0, 1\]. Plain
  Nei's Gst is not offered – it is strongly deflated when within-group
  diversity is high (typical genome-wide in *P. falciparum*); its
  standardized form G'st is provided instead.

- meta:

  When `x` is a matrix, an optional data frame with a `sample` column
  plus `group`; otherwise `group` is a vector aligned to the matrix
  rows.

- clamp:

  Clamp small negative estimates to 0 (default `TRUE`).

- genotype:

  Optional genotype matrix (samples x SNPs) or
  [`load_genotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/load_genotypes.md)
  list to use **instead** of a `PopStructure`'s stored matrix – pass the
  **full, unpruned** set here (e.g.
  `load_genotypes(vcf, prune = FALSE)`) so differentiation is measured
  on every SNP while PCA/UMAP keep using the pruned matrix. Ignored when
  `x` is a matrix.

## Value

A `pop_diff` object: a list with `D` (a SNP x pair matrix of the
statistic), `snp`, `groups`, `pairs`, `statistic`, and the group
`freqs`.

## References

Jost, L. (2008) G_ST and its relatives do not measure differentiation.
*Molecular Ecology* 17, 4015-4026.
[doi:10.1111/j.1365-294X.2008.03887.x](https://doi.org/10.1111/j.1365-294X.2008.03887.x)

Nei, M. & Chesser, R. K. (1983) Estimation of fixation indices and gene
diversities. *Annals of Human Genetics* 47, 253-259.
[doi:10.1111/j.1469-1809.1983.tb00993.x](https://doi.org/10.1111/j.1469-1809.1983.tb00993.x)

Hedrick, P. W. (2005) A standardized genetic differentiation measure.
*Evolution* 59, 1633-1638.
[doi:10.1111/j.0014-3820.2005.tb01814.x](https://doi.org/10.1111/j.0014-3820.2005.tb01814.x)

Hudson, R. R., Slatkin, M. & Maddison, W. P. (1992) Estimation of levels
of gene flow from DNA sequence data. *Genetics* 132, 583-589.
[doi:10.1093/genetics/132.2.583](https://doi.org/10.1093/genetics/132.2.583)

Bhatia, G., Patterson, N., Sankararaman, S. & Price, A. L. (2013)
Estimating and interpreting F_ST: the impact of rare variants. *Genome
Research* 23, 1514-1521.
[doi:10.1101/gr.154831.113](https://doi.org/10.1101/gr.154831.113)
