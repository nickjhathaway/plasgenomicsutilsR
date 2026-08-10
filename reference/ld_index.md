# Multilocus linkage disequilibrium: the index of association

`Ia` and its sample-size-standardised form `rbarD` measure whether
alleles at *unlinked* loci travel together – the genome-wide signature
of clonal propagation. Both are 0 under free recombination and rise as
reproduction becomes more clonal. `rbarD` is the one to compare between
datasets, since `Ia` grows with the number of loci. Both are also
upward-biased in small samples, so a group of ten will score above a
group of a hundred drawn from the same population – read them beside
`n_samples`.

## Usage

``` r
ld_index(
  x,
  group = NULL,
  maf = 0,
  max_snps = LD_MAX_SNPS,
  het = c("missing", "dosage"),
  min_samples = 4,
  meta = NULL,
  genotype = NULL
)
```

## Arguments

- x:

  A
  [PopStructure](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PopStructure.md)
  or a genotype matrix (samples x SNPs, alt dosage, `chr:pos` column
  names).

- group:

  Metadata column name, or a vector aligned to the rows; `NULL` pools
  every sample.

- maf:

  Skip loci below this minor-allele frequency within the group.

- max_snps:

  Loci sampled before the pairwise distances, evenly along the genome.

- het:

  How a heterozygous call is read; see
  [`pop_diversity()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diversity.md).

- min_samples:

  Skip groups smaller than this.

- meta, genotype:

  As in
  [`pop_diversity()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diversity.md).

## Value

A tibble with `group`, `n_samples`, `n_loci`, `ia`, `rbar_d`.

## Details

A multilocus genotype is only defined where every locus is called, so
the loci carrying a gap are dropped and every sample is kept – over
thousands of loci even slight missingness would otherwise leave almost
no complete sample, and it is the samples the statistic is about.
`n_loci` reports what survived. The estimator is the pairwise-distance
decomposition of Brown et al. (1980), standardised by Agapow & Burt
(2001).

## References

Brown, A. H. D., Feldman, M. W. & Nevo, E. (1980) Multilocus structure
of natural populations of *Hordeum spontaneum*. *Genetics* 96, 523-536.
[doi:10.1093/genetics/96.2.523](https://doi.org/10.1093/genetics/96.2.523)

Agapow, P.-M. & Burt, A. (2001) Indices of multilocus linkage
disequilibrium. *Molecular Ecology Notes* 1, 101-102.
[doi:10.1046/j.1471-8278.2000.00014.x](https://doi.org/10.1046/j.1471-8278.2000.00014.x)

## Examples

``` r
ps <- example_pop_structure(umap = FALSE)
ld_index(ps, group = "country")
#> # A tibble: 2 × 5
#>   group    n_samples n_loci    ia  rbar_d
#>   <fct>        <int>  <int> <dbl>   <dbl>
#> 1 Cambodia        30     92 11.8  0.141  
#> 2 Ghana           30    425  2.74 0.00672
```
