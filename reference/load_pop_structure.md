# Load a saved PopStructure workspace

Reads an `.rds` written by `PopStructure$save()` (or plain
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html)).

## Usage

``` r
load_pop_structure(file)
```

## Arguments

- file:

  Path to the `.rds` file.

## Value

The
[PopStructure](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PopStructure.md)
object.
