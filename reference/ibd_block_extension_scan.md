# Genome-wide scan of IBD block extension, and loci placed against it

Runs
[`ibd_block_extension_test()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_block_extension_test.md)
over fixed-width windows tiling the genome, then reports where each
window – and each locus you name – sits in its own group's distribution.
This is the interface to reach for first: a single locus's
`paired_ratio` has no interpretable scale on its own, because the
statistic is inflated at *every* locus by length-biased sampling, and
the scan is what measures that inflation for the cohort in hand.

## Usage

``` r
ibd_block_extension_scan(
  x,
  loci = NULL,
  regions = NULL,
  width = 500,
  step = 2000,
  group = NULL,
  within = 0,
  sharing = c("overlap", "complete"),
  min_ref_blocks = 1L,
  min_pairs = 5L,
  min_windows = 50L,
  meta = NULL
)
```

## Arguments

- x:

  An
  [IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
  built with `blocks =` and `meta =`.

- loci:

  Optional loci to place against the scan: gene names from the object's
  track, or an interval data frame as
  [`ibd_block_extension_test()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_block_extension_test.md)
  takes. Each comes back with its `percentile` in its group's window
  distribution. `NULL` (default) returns the windows alone.

- regions:

  Regions to tile (`chr`/`chrom`, `start`, `end`). Defaults to
  [PF3D7_CORE_REGIONS](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF3D7_CORE_REGIONS.md)
  for a `pf3d7` object, which keeps the hypervariable subtelomeres out
  of the background; pass your own for any other reference.

- width:

  Window width in bp (default `500`).

- step:

  Distance between window starts in bp (default `2000`). See the
  resolution note: finer buys nothing.

- group, within, sharing, min_ref_blocks, min_pairs, meta:

  Passed to
  [`ibd_block_extension_test()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_block_extension_test.md).

- min_windows:

  Groups with fewer than this many tested windows get `NA` percentiles
  and a warning, since a distribution over a handful of windows is not
  one (default `50`).

## Value

A tibble with
[`ibd_block_extension_test()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_block_extension_test.md)'s
columns plus:

- `source`:

  `"window"` for the background rows, `"locus"` for anything named in
  `loci`.

- `percentile`:

  where this row's `paired_ratio` falls among the **windows** of the
  same group, as a percentage at or below it. The number to quote.

- `n_windows`:

  tested windows behind that group's distribution.

Loci come first, then windows by descending `paired_ratio`. The per-pair
ratios are not retained; call
[`ibd_block_extension_test()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_block_extension_test.md)
directly on a locus when you want them.

## Details

Windows rather than genes as the background, for two reasons. They are
all the same width, which matters under `sharing = "complete"`, where a
segment cannot be shorter than the locus is wide and so wider loci score
higher for no biological reason. And they owe nothing to annotation, so
a target between genes still shows up.

## What the resolution is, and is not

Do not read a window as localising selection to its own width. The value
at a point is set by the segments covering it, and those run tens of kb
either side, so neighbouring windows repeat each other almost exactly:
on one real cohort the correlation of `log2(paired_ratio)` between
windows 250 bp apart was 0.997, still 0.973 at 2 kb, and only fell to
0.65 at the median segment length of 33 kb. Every window inside *pfcrt*,
*pfgch1* and *pfkelch13* returned an identical ratio. The resolution
floor is the segment-length scale, so a sub-genic target is out of reach
here – use
[`run_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_ihs.md)
and
[`plot_ehh()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ehh.md),
which work at SNP resolution.

That is also why `step` defaults to 2 kb rather than something finer: a
denser grid costs time and memory and returns the same numbers.

## Why there is no q-value

Windows overlap the same segments, so they are not independent tests.
The core genome over the median segment length is on the order of a few
hundred independent positions against tens of thousands of windows,
which makes a Benjamini-Hochberg q across windows both invalid and
wildly conservative. `percentile` is the quantity to quote; it does not
care how finely the genome was cut. Per-window p-values are still
reported, unadjusted, so you can see which windows are tested at all.

## See also

[`ibd_block_extension_test()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_block_extension_test.md)
for one locus and what the columns mean,
[`ibd_block_extension_by_allele()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_block_extension_by_allele.md)
for splitting the pairs by carriage of a variant.

## Examples

``` r
if (FALSE) { # \dontrun{
ibd <- ibd_results(blocks = "hmm.txt", meta = meta, group_col_in_meta = "region")
scan <- ibd_block_extension_scan(ibd, loci = c("pfcrt", "pfgch1"))
subset(scan, source == "locus", c(locus, group, paired_ratio, percentile))
} # }
```
