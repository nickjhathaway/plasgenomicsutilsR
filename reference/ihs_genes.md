# Summarise a haplotype scan per gene

The strongest signal inside each gene, which is how selection scans are
usually reported: one row per gene rather than per SNP.

## Usage

``` r
ihs_genes(scan, genes = NULL, within = 0, min_snps = 1)
```

## Arguments

- scan:

  The tibble from
  [`run_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_ihs.md),
  [`run_rsb()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_rsb.md)
  or
  [`run_xpehh()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_xpehh.md).

- genes:

  Gene table (`name`, `chr` or `chrom`, `start`, `end`); defaults to
  [PF3D7_GENES](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF3D7_GENES.md).
  Coordinates are 0-based half-open.

- within:

  Widen each gene by this many bp on both sides before matching SNPs.
  Note what this does to a table of "top genes": one strong SNP is then
  credited to every gene within `within` bp of it, so a single sweep can
  fill several rows that share a `peak_pos`. `peak_in_gene` says whether
  the peak is actually inside the gene or was pulled in from the flank,
  and `n_snps` counts the widened window.

- min_snps:

  Genes with fewer scored SNPs than this are dropped.

## Value

A tibble with the grouping column, `gene`, `chr`, `start`, `end`,
`n_snps`, `max_neg_log10_p`, `max_abs_value`, `peak_pos` and
`peak_in_gene`.

## Examples

``` r
ps <- example_pop_structure(umap = FALSE)
hap <- parasite_haplotypes(ps, maf = 0.05)
ihs_genes(run_ihs(hap, group = "country"), genes = PF_EXAMPLE_DRUG_GENES, min_snps = 1)
#> # A tibble: 0 × 0
```
