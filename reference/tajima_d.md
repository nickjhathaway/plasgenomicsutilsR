# Tajima's D from segregating sites and per-site heterozygosity

Tajima's D from segregating sites and per-site heterozygosity

## Usage

``` r
tajima_d(h, n)
```

## Arguments

- h:

  Per-site heterozygosity at each usable site (see
  [`pop_diversity()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diversity.md)).

- n:

  Number of sampled gene copies (haploid: the sample count). A mean over
  sites with unequal missingness is fine; it is rounded for the harmonic
  sums.

## Value

Tajima's D, or `NA` when there are too few segregating sites or samples.

## References

Tajima, F. (1989) Statistical method for testing the neutral mutation
hypothesis by DNA polymorphism. *Genetics* 123, 585-595.
[doi:10.1093/genetics/123.3.585](https://doi.org/10.1093/genetics/123.3.585)
