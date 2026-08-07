# Load a saved PopStructure workspace

Reads an `.rds` written by `PopStructure$save()` (or plain
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html)). The workspace is
re-bound to the installed version of the class, so a file written by an
older version of the package gains the methods added since.

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
