# Merge loaded genotype sets into one

Column-binds the genotype matrices from two or more
[`load_genotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/load_genotypes.md)
results, keyed on the `chr:pos` names the loader assigns.

## Usage

``` r
merge_genotypes(
  ...,
  samples = c("identical", "common"),
  on_overlap = c("error", "first", "last")
)
```

## Arguments

- ...:

  Two or more
  [`load_genotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/load_genotypes.md)
  results.

- samples:

  `"identical"` (default) requires the same sample set in every input;
  `"common"` keeps the intersection.

- on_overlap:

  What to do when the same `chr:pos` appears in more than one input:
  `"error"` (default), or `"first"` / `"last"` to keep that input's
  calls.

## Value

A list shaped like
[`load_genotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/load_genotypes.md)'s:
`genotype`, `sample.id`, `snp.id`, `allele`, `pruned`, `positions`,
`variants`. `pruned` is `TRUE` only when every input was pruned – a SNP
added to a pruned set never faced pruning itself, which is usually the
reason for adding it.

## Details

Merging here rather than upstream is the point. Two callers write
different INFO and FORMAT fields, and reconciling them so
`bcftools merge` will accept both is real work that changes nothing
about the answer – by the time a callset is a dosage matrix, all that
survives is the calls themselves, and those combine by position. So a
variant that had to be re-called separately, because the original
callset arrived already filtered and the caller's own output was gone,
can be added without rebuilding the VCF it came from.

What is checked, because these are the ways a merge is silently wrong:

- **`allele`** must agree. One set counting alternate alleles and
  another counting reference alleles are indistinguishable after the
  fact, and mixing them inverts the dosages of whichever half disagrees.

- **Samples** must be the same set, in whatever order; the columns are
  realigned to the first set's order. `samples = "common"` intersects
  instead, saying how many it dropped.

- **Positions** must not collide. The same `chr:pos` from two callers is
  two answers to one question, and picking silently would hide the
  disagreement – name which one wins with `on_overlap`.

The result is sorted by chromosome and position, so it reads like a
callset rather than like the order the files happened to be given in.
`snp.id` is renumbered: SNPRelate's ids are per-file integers starting
at 1, so they collide across files and mean nothing once merged –
`chr:pos` is the identity that survives.

## See also

[`load_genotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/load_genotypes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
main <- load_genotypes("cohort.vcf.gz", gds = "cohort.gds", prune = TRUE)
extra <- load_genotypes("one_recalled_snp.vcf.gz", gds = "extra.gds", prune = FALSE)
both <- merge_genotypes(main, extra)
} # }
```
