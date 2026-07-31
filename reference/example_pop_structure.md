# Public population-structure example datasets

Builds a
[PopStructure](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PopStructure.md)
object from a bundled **public** genotype matrix. Two datasets ship:

- `"ghana_cambodia"` (default) – 60 Pf7 samples (30 Ghana, 30 Cambodia)
  at 49 biallelic SNPs; two well-separated populations, a clean minimal
  demo.

- `"africa"` – 258 published East/Central-African samples (`country`:
  DRC, Kenya, Tanzania, Uganda; finer `site` sub-regions; a macro
  `region`) at the 2,000 SNPs that most differentiate the sites (top
  Jost's D; see
  [`top_differentiating_snps()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/top_differentiating_snps.md)),
  so the regional structure is clear – a richer, multi-region demo for
  the combined UMAP + admixture figure
  ([`plot_structure_figure()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_structure_figure.md)).

## Usage

``` r
example_pop_structure(
  dataset = c("ghana_cambodia", "africa"),
  umap = TRUE,
  seed = 42
)
```

## Arguments

- dataset:

  Which bundled dataset to load.

- umap:

  Also compute a UMAP embedding (needs uwot); skipped with a message if
  the package is missing.

- seed:

  Random seed for the UMAP.

## Value

A
[PopStructure](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PopStructure.md)
object with a `meta` data frame.

## Examples

``` r
ps <- example_pop_structure(umap = FALSE)
ps
#> <PopStructure> 60 of 60 samples, 49 PCs
#>   UMAP: -   sNMF: -   meta: country 
```
