# Do two populations carry the same haplotype at a locus? An allele-split IBS test

Measures whether the carriers of a variant in two populations are drawn
from a common haplotype pool, using pairwise identity (IBS) rather than
called IBD, so it never depends on the caller reaching the
cross-population pairs the question turns on. The statistic is a
difference of differences – cross-population carrier-vs-carrier identity
minus reference-vs-reference identity – scored against a null of
SNP-count-matched windows tiled off the focal chromosome. It answers
"one haplotype spreading across the boundary" versus "independent
origins", which the IBD block-extension family (the iHS-shaped, within-
population question; see
[`ibd_block_extension_by_allele()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_block_extension_by_allele.md))
cannot settle on its own.

## Usage

``` r
locus_ibs_by_allele(
  x,
  locus,
  allele,
  carrier = NULL,
  reference = NULL,
  group_a,
  group_b,
  group = NULL,
  meta = NULL,
  n_snps = 53L,
  span_tol = 2,
  exclude_focal_chr = TRUE,
  min_sites = 10L,
  drop_ibd = NULL,
  map = NULL
)
```

## Arguments

- x:

  Genotypes: a
  [parasite_haplotypes](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/parasite_haplotypes.md)
  object, a
  [`load_genotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/load_genotypes.md)
  list (its `allele-index` encoding, which keeps missingness), or a bare
  integer matrix (samples x SNPs) with `map`. Missing calls must be
  `NA`; an imputed panel is warned about, since imputation would be read
  as agreement. For IBS and IBD to be comparable the calls should be the
  ones the IBD run used (hmmibd-rs dominant-allele from `FORMAT/AD`);
  the allele-index panel is `GT`-derived and agrees with it wherever
  calls are unmixed.

- locus:

  The focal locus: a one-row data frame or list with `chr`, `start`,
  `end` (0-based half-open, used for the `drop_ibd` overlap) and
  optional `centre` (the focal window is centred here; defaults to the
  interval midpoint). A multi-row interval frame is collapsed to its
  span.

- allele:

  The variant to split on: a metadata column name or a named
  `sample -> state` vector (as
  [`allele_states()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/allele_states.md)
  returns).

- carrier, reference:

  State sets (matched with `%in%`), as in
  [`carrier_genome_wide_relatedness()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/carrier_genome_wide_relatedness.md).

- group_a, group_b:

  The two populations, each a set of `group` values. The contrast is
  between them; a sample in neither is unused.

- group:

  Metadata column defining the populations (default `"region"`, or the
  object's own group column when it carries one).

- meta:

  Sample metadata (`sample` plus `group` and the `allele` column). Taken
  from a `parasite_haplotypes`' `$meta` when not given.

- n_snps:

  Window size in SNPs (default `53`).

- span_tol:

  Keep null windows whose physical span is within this factor of the
  focal window's span, either side (default `2`; `NULL` or `0` to not
  bound).

- exclude_focal_chr:

  Drop the focal chromosome from the null (default `TRUE`).

- min_sites:

  Comparable sites a pair needs in a window to contribute (default
  `10`).

- drop_ibd:

  Optionally also report the statistic with pairs already IBD across the
  locus removed from both strata, which separates "a lineage moved
  across the border" from "the carrier populations share a background".
  Pass an
  [IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
  (its blocks name the IBD pairs), or `TRUE` to use `x` if it is one.
  Adds a second result row.

- map:

  When `x` is a bare genotype matrix, a data frame of `chr`/`pos` with
  one row per SNP column giving each marker's position; ignored when `x`
  carries its own `$map`.

## Value

A tibble, one row per pair set (`"all pairs"`, and `"IBD pairs dropped"`
when `drop_ibd` is given), carrying the per-window null on the `null`
attribute:

- `locus`, `pairs`:

  the locus label and which pair set the row is.

- `n_a_car`, `n_b_car`, `n_a_ref`, `n_b_ref`:

  samples in each population x stratum.

- `cross_car`, `cross_ref`:

  cross-population mean IBS in each stratum.

- `a_car`, `b_car`:

  within-population carrier cohesion, each population on its own.

- `d_ind`:

  `cross_car - cross_ref`, the statistic.

- `null_n`, `null_mean`, `null_sd`:

  the matched-window null.

- `z`:

  `(d_ind - null_mean) / null_sd`.

- `p_two`:

  two-sided empirical p against the null (a floor when small).

- `focal_span_kb`, `span_tol`, `n_snps`:

  the window the null was matched to.

## What it does and does not say

A positive, significant result means the carrier chromosomes in the two
populations are drawn from a common pool rather than independently
derived. It does **not** count how many backgrounds that pool holds;
pair it with a tabulation of the other coding variants each carrier also
carries. A negative difference means the carriers are *less* alike
across the boundary than the reference samples are – independent
backgrounds – even where each population's carriers cluster tightly on
their own (`a_car`, `b_car`).

## Matched windows

Windows are matched on **SNP count** (`n_snps`), not physical width,
because IBS is a mean over a window's SNPs and its variance is set by
that count. A physical-span bound (`span_tol`) is layered on top to stop
a count-matched window that happens to span far more bp from inflating
the null variance. The focal chromosome is excluded from the null so the
locus cannot seed its own baseline, and a pair needs `min_sites`
comparable sites in a window to contribute. With a few hundred windows
the smallest achievable two-sided p is a resolution floor, so read a p
at the floor as "no matched window reached the observed value" rather
than as an exact probability.

## See also

[`ibd_block_extension_by_allele()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_block_extension_by_allele.md)
for the within-population iHS-shaped contrast,
[`carrier_genome_wide_relatedness()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/carrier_genome_wide_relatedness.md)
for the genome-wide relatedness control,
[`parasite_haplotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/parasite_haplotypes.md)
and
[`load_genotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/load_genotypes.md)
for the genotype inputs.
