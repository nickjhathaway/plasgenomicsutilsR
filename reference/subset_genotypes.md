# Restrict a genotype panel to a set of samples

Narrows the rows and leaves everything else alone, returning the same
kind of object it was given – a
[`load_genotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/load_genotypes.md)
list comes back a list, with `sample.id` narrowed to match and `snp.id`
/ `allele` / `pruned` / `positions` / `variants` carried through.

## Usage

``` r
subset_genotypes(x, samples, strict = FALSE)
```

## Arguments

- x:

  A
  [`load_genotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/load_genotypes.md)
  list, a genotype matrix (samples x SNPs), or a
  [PopStructure](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PopStructure.md)
  (which is narrowed with its own `$subset()`).

- samples:

  The samples to keep: a character vector of ids, a
  [`parasite_haplotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/parasite_haplotypes.md)
  object (its kept samples), a `PopStructure`, another
  [`load_genotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/load_genotypes.md)
  list, or a metadata data frame with a `sample` column.

- strict:

  Error when `samples` names ids the panel does not hold (default
  `FALSE`, which warns and keeps the intersection).

## Value

The same kind of object as `x`, holding only those samples, in the
panel's own order.

## Details

The use it was written for is running a diversity or differentiation
analysis on the monoclonals a haplotype set kept:

    mono <- parasite_haplotypes(geno$genotype, fws = fws, min_fws = 0.92)
    pop_diversity(subset_genotypes(geno, mono)$genotype, group = "region", meta = meta)

Note it subsets **samples**, not SNPs: the full SNP panel is kept, which
is what a diversity estimate wants. Handing `hap$hap` to such an
analysis instead would silently use the MAF-filtered, imputed 0/1 matrix
the scans need, which is a different panel.

## See also

[`haplotype_samples()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/haplotype_samples.md),
[PopStructure](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PopStructure.md)'s
`$subset()`.

## Examples

``` r
ps <- example_pop_structure(umap = FALSE)
hap <- parasite_haplotypes(ps, maf = 0.05)

# a bare matrix
dim(subset_genotypes(ps$genotype(), hap))
#> [1] 60 49

# a load_genotypes()-shaped list keeps its other slots
geno <- list(genotype = ps$genotype(), sample.id = ps$get_samples(),
             snp.id = colnames(ps$genotype()), pruned = TRUE, positions = "0-based")
str(subset_genotypes(geno, hap)[c("sample.id", "pruned", "positions")], max.level = 1)
#> List of 3
#>  $ sample.id: chr [1:60] "PH1708-C" "PH1729-C" "PH1747-C" "RCN07826" ...
#>  $ pruned   : logi TRUE
#>  $ positions: chr "0-based"
```
