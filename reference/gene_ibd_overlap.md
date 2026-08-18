# Per-gene IBD-block overlap between groups

For each gene and each pair of groups, the fraction of sample pairs that
share an IBD **block overlapping the gene** – a pair counts when any of
its IBD segments overlaps the (optionally padded) gene interval, so a
segment spanning the gene counts even with no genotyped SNP inside it.
The denominator is *all* pairs compared in that group-pair (from the
analyzed-sample set), so pairs that are never IBD still count against
the fraction.

## Usage

``` r
gene_ibd_overlap(x, genes = NULL, group = NULL, within = 0)
```

## Arguments

- x:

  An
  [IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
  with IBD `blocks` and `meta`.

- genes:

  Gene names from the object's track (default all), or a gene-interval
  data frame (`name`, `chr`/`chrom`, `start`, `end`) to use instead of
  the track.

- group:

  Metadata column defining the groups (default: the first non-`sample`
  column of `meta`).

- within:

  Pad each gene interval by this many bp on both sides (default `0`).

## Value

A tibble, one row per gene x group-pair: `gene` (a unique display label
– the `name`, disambiguated as `name (gene_id)` where GFF Names repeat
across gene families), `name`, `gene_id`, `chr`, `start`, `end`,
`group_a`, `group_b`, `n_pairs_ibd`, `n_pairs_total`, `frac_pairs_ibd`.

## Details

Needs an
[IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
built with `blocks =` and `meta =` (or use the Python `ibd_gene_overlap`
tool and pass its table as `gene_overlap =`).

## Examples

``` r
# needs IBD segments, so build a small object that has them
# segments shorter than 15 kb or carrying under 15 SNPs are discarded as spurious,
# so these are long enough to survive that
blocks <- data.frame(
  sample1 = c("s1", "s1", "s2"), sample2 = c("s2", "s3", "s3"),
  chr = "7", start = c(395000, 100000, 395000), end = c(415000, 130000, 415000),
  different = 0, Nsnp = 40)
meta <- data.frame(sample = paste0("s", 1:3), region = c("north", "north", "south"))
ibd <- ibd_results(blocks = blocks, meta = meta, group_col_in_meta = "region",
                   genes = PF_EXAMPLE_DRUG_GENES)

# the fraction of pairs sharing the gene by IBD block, per group pair
gene_ibd_overlap(ibd, genes = "pfcrt")
#> # A tibble: 3 × 11
#>   gene  name  gene_id       chr    start    end group_a group_b n_pairs_ibd
#>   <fct> <chr> <chr>         <chr>  <dbl>  <dbl> <fct>   <fct>         <int>
#> 1 pfcrt pfcrt PF3D7_0709000 7     403221 406317 north   north             1
#> 2 pfcrt pfcrt PF3D7_0709000 7     403221 406317 north   south             1
#> 3 pfcrt pfcrt PF3D7_0709000 7     403221 406317 south   south             0
#> # ℹ 2 more variables: n_pairs_total <dbl>, frac_pairs_ibd <dbl>
```
