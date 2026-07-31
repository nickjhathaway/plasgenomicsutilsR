# Pf3D7 drug-resistance gene coordinates

A small `data.frame` of well-known *P. falciparum* drug-resistance genes
(`name`, `chr`, `start`, `end`; 1-based, Pf3D7 assembly) for use as the
`genes` track of an
[IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
object (gene reference lines on Manhattan plots, and the drug-gene
triangles). These are public reference-genome coordinates.

## Usage

``` r
EXAMPLE_DRUG_GENES
```

## Format

A data frame with columns `name`, `chr`, `start`, `end`.
