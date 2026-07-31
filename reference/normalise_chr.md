# Normalise a chromosome name to a bare number string

`"Pf3D7_07_v3"`, `"chr7"`, `"07"`, `7` all become `"7"`.

## Usage

``` r
normalise_chr(c)
```

## Arguments

- c:

  A chromosome name (character or numeric), scalar or vector.

## Value

Character vector of normalised chromosome ids.
