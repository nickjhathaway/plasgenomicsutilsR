# Normalise a chromosome name to a bare number string

`"Pf3D7_07_v3"`, `"chr7"`, `"07"`, `7` all become `"7"`. The organelles
follow the bundled datasets: the apicoplast (`"Pf3D7_API_v3"`) becomes
`"API"` and the mitochondrion (`"Pf_M76611"`, `"Pf3D7_MIT_v3"`) `"MIT"`,
so a table read from a FASTA-named file joins
[PF3D7_GENES](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF3D7_GENES.md)
on every sequence, not just the fourteen chromosomes.

## Usage

``` r
normalise_chr(c)
```

## Arguments

- c:

  A chromosome name (character or numeric), scalar or vector.

## Value

Character vector of normalised chromosome ids.

## Examples

``` r
# every spelling of a chromosome collapses to the same key, so tables from
# different tools join
normalise_chr(c("Pf3D7_07_v3", "chr7", "7", "Pf3D7_API_v3", "Pf_M76611"))
#> [1] "7"   "7"   "7"   "API" "MIT"
```
