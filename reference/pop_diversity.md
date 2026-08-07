# Within-population genetic diversity

Nucleotide diversity, expected heterozygosity, Watterson's theta,
Tajima's D and haplotype diversity, computed for each metadata group and
reported genome-wide, per gene, or in sliding windows.

## Usage

``` r
pop_diversity(
  x,
  group = NULL,
  by = c("genome", "gene", "window"),
  genes = NULL,
  window = 10000,
  step = NULL,
  accessible = NULL,
  het = c("missing", "dosage"),
  min_snps = DIVERSITY_MIN_SNPS,
  max_missing = 0.1,
  min_samples = 4,
  genotype = NULL,
  meta = NULL
)
```

## Arguments

- x:

  A
  [PopStructure](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PopStructure.md),
  or a genotype matrix (samples x SNPs, alt dosage 0/1/2, `NA` for
  missing) with `chr:pos` column names.

- group:

  Metadata column naming the grouping (for a `PopStructure`) or a vector
  aligned to the matrix rows. `NULL` treats every sample as one
  population.

- by:

  `"genome"` (default), `"gene"`, or `"window"`.

- genes:

  Gene table for `by = "gene"` (columns `name`, `chr`, `start`, `end`);
  defaults to
  [PF3D7_GENES](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF3D7_GENES.md).
  Coordinates are 0-based half-open CDS spans.

- window, step:

  Window size and step in bp for `by = "window"`. `step` defaults to
  `window`, giving abutting windows; set it smaller for a **sliding**
  scan – `window = 5000, step = 2500` steps a 5 kb window along in 2.5
  kb hops, so consecutive windows share half their span. Sliding is the
  usual way to scan Tajima's D: a fixed grid can split a signal across
  two windows and dilute it in both, and overlapping windows also make
  the track read more smoothly. The cost is that neighbouring windows
  are no longer independent – do not count them as separate findings
  (see
  [`selection_peaks()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/selection_peaks.md)).

- accessible:

  Callable regions as a data frame with `chrom`, `start`, `end` – e.g.
  [PF3D7_CORE_REGIONS](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF3D7_CORE_REGIONS.md).
  Sets the denominator of `pi` and `theta_w`.

- het:

  How to read a heterozygous call: `"missing"` or `"dosage"`.

- min_snps:

  Fewest SNPs for Tajima's D (default 3).

- max_missing:

  Per-SNP missingness allowed into the Tajima's D block.

- min_samples:

  Skip a group with fewer samples than this.

- genotype:

  Genotype matrix to use instead of a `PopStructure`'s stored one – pass
  the **full, unpruned** set, since LD-pruning removes the very sites
  diversity is measured over.

- meta:

  When `x` is a matrix, a data frame with `sample` plus `group`.

## Value

A tibble with one row per group x unit: `group`, the unit's identity,
then `n_samples`, `n_snps`, `n_sites`, `seg_sites`, `he`, `pi`,
`theta_w`, `tajima_d`, `n_hap`, `hap_div`, `shannon_h`,
`simpson_lambda`, `evenness`. `tajima_p` is the two-sided test against
the standard neutral model and `tajima_percentile` places each unit's D
within its own group's distribution – read
[`tajima_d_pvalue()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/tajima_d_pvalue.md)
before using the former, which is conservative under recombination.

## Details

**`pi` is per accessible base pair** – the sum of per-site
heterozygosity over the number of callable sites (`n_sites`), so windows
of different SNP density stay comparable. Pass `accessible` to say which
bases were callable; without it every base of the unit is assumed
callable, which inflates nothing but does assume your VCF was called
across the whole span. `he` is the same per-site quantity averaged over
the **SNPs** instead, which is what papers usually label expected
heterozygosity (Hs/He); it is not nucleotide diversity and the two are
never comparable.

The parasite is haploid: `het = "missing"` (the default) treats a
heterozygous call as a mixed infection and drops it at that site, so an
allele count is a count of samples. `het = "dosage"` instead splits the
call evenly between the two alleles.

Tajima's D needs a single sample count, so it is computed on a
complete-case block – SNPs missing in more than `max_missing` of the
group are dropped, then any sample still missing a call. `n_taj_snps` /
`n_taj_samples` report what survived, and the statistic is `NA` below
`min_snps` SNPs or four samples.

## References

Nei, M. (1987) *Molecular Evolutionary Genetics*. Columbia University
Press.

Watterson, G. A. (1975) On the number of segregating sites in genetical
models without recombination. *Theoretical Population Biology* 7,
256-276.
[doi:10.1016/0040-5809(75)90020-9](https://doi.org/10.1016/0040-5809%2875%2990020-9)

Korunes, K. L. & Samuk, K. (2021) pixy: Unbiased estimation of
nucleotide diversity and divergence in the presence of missing data.
*Molecular Ecology Resources* 21, 1359-1368.
[doi:10.1111/1755-0998.13326](https://doi.org/10.1111/1755-0998.13326)

## See also

[`pop_diff()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff.md)
for between-group differentiation,
[`ld_index()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ld_index.md)
and
[`read_ld_decay()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/read_ld_decay.md)
for linkage disequilibrium,
[`run_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_ihs.md)
for selection.

## Examples

``` r
ps <- example_pop_structure(umap = FALSE)
pop_diversity(ps, group = "country")
#> # A tibble: 2 × 24
#>   group chr   start   end unit  n_samples n_snps n_sites seg_sites    he      pi
#>   <fct> <chr> <dbl> <dbl> <chr>     <int>  <int>   <dbl>     <int> <dbl>   <dbl>
#> 1 Camb… NA       NA    NA geno…        30     49  2.33e7        24 0.159 3.35e-7
#> 2 Ghana NA       NA    NA geno…        30     49  2.33e7        36 0.194 4.08e-7
#> # ℹ 13 more variables: theta_w <dbl>, tajima_d <dbl>, tajima_p <dbl>,
#> #   n_taj_snps <int>, n_taj_samples <dbl>, n_hap_snps <int>,
#> #   n_hap_samples <int>, n_hap <int>, hap_div <dbl>, shannon_h <dbl>,
#> #   simpson_lambda <dbl>, evenness <dbl>, tajima_percentile <dbl>
```
