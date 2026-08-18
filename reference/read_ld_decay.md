# Read an LD-decay table

Reads what `plasgenomicsutils ld_decay` writes, ready for
[`plot_ld_decay()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ld_decay.md).
The half-decay distance is attached as an `ld_half_decay` attribute –
recomputed from the binned means if the separate table is not supplied –
and the scan's `max_dist`, `max_snps` and `maf` are attached from the
file header, so the thinning behind a curve travels with it.

## Usage

``` r
read_ld_decay(path, half_decay = NULL, levels = NULL)
```

## Arguments

- path:

  TSV(.gz) written by `plasgenomicsutils ld_decay`.

- half_decay:

  Optional TSV from `--half-decay-output`; recomputed when absent.

- levels:

  Optional group order; defaults to a natural sort.

## Value

A tibble of `group`, `bin_start`, `bin_end`, `bin_mid`, `n_pairs`,
`mean_r2`, `median_r2`.

## Details

The decay scan lives in the Python package because it is quadratic in
the number of SNPs inside each window; on a 249-sample, 28k-SNP callset
it runs in about 7 seconds there.
[`ld_index()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ld_index.md)
stays here, being linear and quick.

## See also

[`ld_index()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ld_index.md),
[`plot_ld_decay()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ld_decay.md)

## Examples

``` r
f <- tempfile(fileext = ".tsv")
write.table(data.frame(group = "north", bin_mid = c(500, 1500, 2500),
                       mean_r2 = c(0.40, 0.22, 0.14)),
            f, sep = "\t", quote = FALSE, row.names = FALSE)
read_ld_decay(f)
#> # A tibble: 3 × 3
#>   group bin_mid mean_r2
#>   <fct>   <dbl>   <dbl>
#> 1 north     500    0.4 
#> 2 north    1500    0.22
#> 3 north    2500    0.14
```
