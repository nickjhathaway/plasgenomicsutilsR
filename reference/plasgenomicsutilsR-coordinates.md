# Genomic coordinate conventions

**Everything is 0-based.** There is one rule, so nothing has to be
converted per call and there is no part of the package you have to
remember an exception for:

- **Intervals are half-open `[start, end)`** – the BED convention. This
  covers the bundled datasets
  ([PF3D7_GENES](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF3D7_GENES.md),
  [PF_EXAMPLE_DRUG_GENES](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF_EXAMPLE_DRUG_GENES.md),
  [PF3D7_CORE_REGIONS](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF3D7_CORE_REGIONS.md),
  [PF3D7_PARALOG_GENES](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF3D7_PARALOG_GENES.md)),
  any `genes` track you supply, the IBD `blocks`, `locus =` arguments,
  and
  [`bed_intersect()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/bed_intersect.md).
  `end - start` is the width in bp, and intervals that merely touch
  (`end1 == start2`) do not overlap.

- **Variant positions are 0-based too**, including the `pos` column and
  the `chr:pos` `snp_id` of the per-SNP tables, and the `snps =`
  arguments that select them. A variant sits inside a gene when
  `pos >= start & pos < end`.

Formats that number differently are converted once, at the boundary, and
never leak inward: the PlasmoDB GFF (1-based inclusive) is shifted where
the gene datasets are built in `data-raw/`; `hmmibd-rs` blocks (0-based,
both endpoints inclusive) get `end + 1` when an
[IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
reads them; and VCF `POS` becomes `POS - 1` in the companion Python
package before any table reaches R.

A `snp_id` is therefore `chr:pos0`, one less than the position the VCF
or a genome browser shows. Pass `--with-pos-vcf` to the Python tools to
carry the 1-based position alongside as `pos_vcf` when you want to look
variants up by eye. Ids already present in an input file are never
trusted as keys, since `bcftools annotate --set-id` may have written
either `%POS` or `%POS0` and the file does not record which.
