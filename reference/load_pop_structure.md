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

## Examples

``` r
ps <- example_pop_structure(umap = FALSE)
f <- tempfile(fileext = ".rds")
ps$save(f)
load_pop_structure(f)
#> <PopStructure> 60 of 60 samples, 49 PCs
#>   UMAP: -   sNMF: -   meta: country 
```
