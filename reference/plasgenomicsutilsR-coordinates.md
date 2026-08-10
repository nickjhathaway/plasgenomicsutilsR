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
reads them; and VCF `POS` becomes `POS - 1` both in the companion Python
package before any table reaches R and in
[`load_genotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/load_genotypes.md),
which is where positions enter on the R side (SNPRelate reports 1-based
`POS`). That last one matters more than it looks: a scan built from
those genotypes –
[`run_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_ihs.md),
[`beta_score()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/beta_score.md),
[`pop_diff_snps()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff_snps.md)
– would otherwise sit one base off every IBD table and every interval,
which breaks exact joins between tracks and mis-assigns a SNP sitting on
a gene boundary. A
[PopStructure](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PopStructure.md)
records the convention (`$positions()`); a matrix built elsewhere, or by
a version before this, can be declared with
`PopStructure$new(..., one_based = TRUE)` and is shifted on the way in.

**One deliberate exception: amino-acid positions are 1-based**, counting
the initiator methionine as residue 1, both where you supply them
([`aa_intervals()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/aa_intervals.md))
and where they are reported
([`snp_aa_positions()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/snp_aa_positions.md))
– as is `codon_base`, which runs 1/2/3 through the codon in transcript
orientation. Residue numbering is too settled in the literature to
renumber: *pfcrt* K76T is codon 76 everywhere it is written down, and a
package that called it 75 would be wrong in the only sense that matters.
Verified against the sequence – codon 1 of *pfcrt* is `403,222-403,224`,
which reads `ATG`, and codon 1 of the minus-strand *pfkelch13* is
`1,726,995-1,726,997`, `CAT` on the forward strand and so `ATG` read
along the transcript. Genomic coordinates in those same tables (`start`,
`end`) stay 0-based half-open like everything else; `codon_positions` is
the one place a 1-based *position* is reported, because those are the
numbers people quote.

A `snp_id` is therefore `chr:pos0`, one less than the position the VCF
or a genome browser shows. Pass `--with-pos-vcf` to the Python tools to
carry the 1-based position alongside as `pos_vcf` when you want to look
variants up by eye. Ids already present in an input file are never
trusted as keys, since `bcftools annotate --set-id` may have written
either `%POS` or `%POS0` and the file does not record which.
