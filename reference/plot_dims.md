# Suggested output dimensions for an IBD plot

Returns a `c(width, height)` (inches) that scales with how much a plot
will draw: the genome fraction shown across (chromosomes kept), and the
number of group panels / triangle features down. The `plot_*()`
functions already attach this to the plot they return, so
[`save_plot()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/save_plot.md)
uses it automatically; call `plot_dims()` yourself only to inspect or
override the numbers.

## Usage

``` r
plot_dims(
  x,
  type = c("manhattan", "selection", "tugofwar", "heatmap", "triangles"),
  groups = NULL,
  chroms = NULL,
  skip_chr = NULL,
  genes = NULL,
  snps = NULL,
  ncol = NULL,
  label_genes = FALSE
)
```

## Arguments

- x:

  An
  [IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
  object.

- type:

  One of `"manhattan"`, `"selection"`, `"tugofwar"`, `"heatmap"`,
  `"triangles"`.

- groups, chroms, skip_chr, genes, snps, ncol, label_genes:

  The same selection arguments you pass to the plot, so the counts match
  what will be drawn.

## Value

A named numeric vector `c(width = , height = )` in inches.

## Examples

``` r
plot_dims(example_ibd_results(), "tugofwar")
#>  width height 
#>   16.0    9.4 
```
