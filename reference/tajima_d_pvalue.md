# p-value for a Tajima's D

Two-sided significance against the standard neutral model, either by
Tajima's own beta approximation (the default, and what pegas reports as
`Pval.beta`) or by treating D as standard normal (`Pval.normal`). Both
agree with pegas to floating point.

## Usage

``` r
tajima_d_pvalue(D, n, S, method = c("beta", "normal"))
```

## Arguments

- D:

  Tajima's D.

- n:

  Number of sampled gene copies.

- S:

  Number of segregating sites.

- method:

  `"beta"` (Tajima's approximation over D's attainable range) or
  `"normal"`.

## Value

A two-sided p-value, or `NA` when D is undefined.

## Details

**Read it as conservative, and low-powered.** The variance term behind D
assumes *no recombination*, which is the most variance the statistic can
have; a heavily recombining organism like *P. falciparum* therefore
under-rejects. Measured on a real 249-sample cohort, per-gene D reached
`p < 0.05` for only 1-2.5% of genes – less than chance would give –
while 72% of genes had negative D. Per gene there are usually only a
handful of segregating sites, which leaves almost no power.

What the test asks is also not quite what a scan wants. It tests the
standard neutral model as a whole, so a population that has expanded
gives negative D genome-wide and a small p-value says "not a
constant-size neutral population", not "this gene is under selection".
For shortlisting loci, rank on `tajima_d` itself, or use the
`tajima_percentile` that
[`pop_diversity()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diversity.md)
reports, and keep the p-value as context.

## References

Tajima, F. (1989) Statistical method for testing the neutral mutation
hypothesis by DNA polymorphism. *Genetics* 123, 585-595.
[doi:10.1093/genetics/123.3.585](https://doi.org/10.1093/genetics/123.3.585)

## Examples

``` r
tajima_d_pvalue(-2.1, n = 40, S = 25)
#> [1] 0.01500924
```
