# The amino acid a SNP falls in

The reverse of
[`aa_intervals()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/aa_intervals.md):
given SNP positions and a GFF, which codon of which transcript each one
sits in. That is how a variant gets talked about – "*pfdhps* A437G" –
and it is the step between "these SNPs are in this gene" and "these are
the amino acids they change".

## Usage

``` r
snp_aa_positions(snps, gff, keep = c("all", "hits"), one_based_snps = FALSE)
```

## Arguments

- snps:

  A data frame with `snp_id` (`"chr:pos"`) or `chr` and `pos` columns;
  any other columns are carried through.

- gff:

  A GFF path, or the result of
  [`read_gff_cds()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/read_gff_cds.md)
  (parse once, reuse).

- keep:

  `"all"` (default) keeps non-coding SNPs with `NA` annotations;
  `"hits"` keeps only those in a CDS.

- one_based_snps:

  The positions in `snps` are 1-based (VCF `POS`). Genotype-matrix
  column names from
  [`load_genotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/load_genotypes.md)
  are, so set this when feeding those in; the package's own tables are
  0-based, which is the default.

## Value

`snps` with `transcript_id`, `gene_id`, `aa_position` (1-based),
`codon_base` (1/2/3 in transcript orientation), `strand` and `coding`
added.

## Details

`aa_position` is **1-based**, counting the initiator methionine as 1,
because that is how the literature numbers residues. `codon_base` says
which of the codon's three bases the SNP is, in transcript orientation,
so on a minus-strand gene `codon_base == 1` is the *highest* genomic
coordinate of the three.

A SNP outside any CDS gets `NA` (or is dropped by `keep = "hits"`) –
introns, UTRs and intergenic space are all simply non-coding here. A SNP
inside overlapping isoforms yields one row per transcript, since the
codon it hits can differ between them.

Being in a codon says nothing about whether the residue actually
changes: that needs the alleles and the reading frame's other two bases,
which are not part of this.

## See also

[`aa_intervals()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/aa_intervals.md)
for the other direction,
[`annotate_snps()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/annotate_snps.md)
to first ask which gene a SNP is in.

## Examples

``` r
if (FALSE) { # \dontrun{
cds <- read_gff_cds("Pf3D7.gff")
# every genotyped SNP in pfdhps, and the residues they sit on
snps <- data.frame(snp_id = colnames(ps$genotype("full")))
snp_aa_positions(snps, cds, keep = "hits", one_based_snps = TRUE)
} # }
```
