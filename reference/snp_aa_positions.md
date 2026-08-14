# The amino acid a SNP falls in

The reverse of
[`aa_intervals()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/aa_intervals.md):
given SNP positions and a GFF, which codon of which transcript each one
sits in. That is how a variant gets talked about – "*pfdhps* A437G" –
and it is the step between "these SNPs are in this gene" and "these are
the amino acids they change".

## Usage

``` r
snp_aa_positions(
  snps,
  gff,
  keep = c("all", "hits"),
  one_based_snps = FALSE,
  fasta = NULL
)
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

- fasta:

  Optional reference sequence, for the `ref_codon` / `ref_aa` columns: a
  path or URL to a FASTA, or a named character vector of sequences.
  Defaults to whatever `gff` carried; `NULL` with a GFF holding no
  sequence leaves the two columns off.

## Value

`snps` with `transcript_id`, `gene_id`, `aa_position` (1-based),
`codon_base` (1/2/3 in transcript orientation), `strand` and `coding`
added, plus `ref_codon` and `ref_aa` when there is sequence to read them
from.

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
changes: that needs the alleles, which are not part of this. What
reference sequence does buy you is the residue the codon currently codes
for – see below.

## The reference residue

Given sequence, the result also carries `ref_codon` (the codon's three
bases, in transcript orientation, complemented on the minus strand) and
`ref_aa` (its residue, one letter, `*` for a stop). Sequence can come
from either place:

- `fasta =` – a path or URL to a genome FASTA (plain or gzipped), or a
  named vector of sequences you already have in hand.

- the GFF itself, when it ends with a `##FASTA` section.
  [`read_gff_cds()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/read_gff_cds.md)
  keeps those sequences, so nothing extra needs passing.

Without either, both columns are simply absent: the reference base
cannot be inferred from an annotation alone. Some GFFs do carry a
translated protein in a feature attribute, but too few agree on how to
make it worth reading, so it is not used.

Names are matched through
[`normalise_chr()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/normalise_chr.md),
so a FASTA headed `Pf3D7_07_v3` lines up with a GFF spelling it the same
way or differently. If nothing matches, both columns come back `NA` with
a warning naming what was compared – silence there would look like a
genome with no coding SNPs.

This is worth doing even when you think you know the answer. `ref_aa` is
what the *reference* carries, and a reference is one isolate's genome –
not a consensus, and not the ancestral or wild-type sequence. Pf3D7
reads `G` at *pfdhps* 437, the residue A437G is named for changing *to*,
so a callset aligned to 3D7 shows no variant at that position precisely
because the reference is already the mutant. This is not a quirk of one
locus: across *Plasmodium* species the reference is sometimes the
non-wild-type allele, so "REF" means "what this isolate has", never
"what came first".

Two conventions here follow from that.
[`load_genotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/load_genotypes.md)
records which allele its dosages count rather than letting REF stand in
for a baseline, and
[`run_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_ihs.md)
and
[`plot_ehh()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ehh.md)
default to `polarized = FALSE`, treating the two states as simply the
two states, since ancestral versus derived cannot be read off REF and
ALT.

## See also

[`aa_intervals()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/aa_intervals.md)
for the other direction,
[`annotate_snps()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/annotate_snps.md)
to first ask which gene a SNP is in.

## Examples

``` r
cds <- read_gff_cds(system.file("extdata", "pf3d7_drug_gene_cds.gff",
                                package = "plasgenomicsutilsR"))

# the three bases of pfkelch13 codon 580, on the minus strand: base 1 is the highest
snp_aa_positions(data.frame(chr = "Pf3D7_13_v3", pos = c(1725260, 1725259, 1725258)),
                 cds, keep = "hits", one_based_snps = TRUE)[
  , c("pos", "aa_position", "codon_base", "strand")]
#> # A tibble: 3 × 4
#>       pos aa_position codon_base strand
#>     <dbl>       <int>      <int> <chr> 
#> 1 1725260         580          1 -     
#> 2 1725259         580          2 -     
#> 3 1725258         580          3 -     

# `ref_codon` / `ref_aa` need sequence. A GFF ending in a `##FASTA` section carries its own,
# so nothing extra is passed:
gff <- tempfile(fileext = ".gff")
writeLines(c("##gff-version 3",
             "demo\t.\tCDS\t11\t25\t.\t+\t0\tID=c1;Parent=T.1;gene_id=T",
             "##FASTA", ">demo", "CCCCCCCCCCATGAAATTTGGGTAAC"), gff)
snp_aa_positions(data.frame(chr = "demo", pos = 11:22), read_gff_cds(gff),
                 keep = "hits", one_based_snps = TRUE)[
  , c("pos", "aa_position", "codon_base", "ref_codon", "ref_aa")]
#> 1 sequence(s) read from the GFF's own ##FASTA section
#> # A tibble: 12 × 5
#>      pos aa_position codon_base ref_codon ref_aa
#>    <int>       <int>      <int> <chr>     <chr> 
#>  1    11           1          1 ATG       M     
#>  2    12           1          2 ATG       M     
#>  3    13           1          3 ATG       M     
#>  4    14           2          1 AAA       K     
#>  5    15           2          2 AAA       K     
#>  6    16           2          3 AAA       K     
#>  7    17           3          1 TTT       F     
#>  8    18           3          2 TTT       F     
#>  9    19           3          3 TTT       F     
#> 10    20           4          1 GGG       G     
#> 11    21           4          2 GGG       G     
#> 12    22           4          3 GGG       G     

if (FALSE) { # \dontrun{
# every genotyped SNP that is coding, and the residue it sits on. Genotype-matrix ids
# are 0-based, so no `one_based_snps` here.
snp_aa_positions(data.frame(snp_id = colnames(ps$genotype("full"))), cds, keep = "hits")

# or point `fasta` at the released genome, read straight from the web like the GFF is.
# pfcrt codon 76 comes back "AAA" / "K".
genome <- paste0("https://plasmodb.org/common/downloads/Current_Release/",
                 "Pfalciparum3D7/fasta/data/PlasmoDB-68_Pfalciparum3D7_Genome.fasta")
snp_aa_positions(data.frame(chr = "Pf3D7_07_v3", pos = 403625), cds, keep = "hits",
                 one_based_snps = TRUE, fasta = genome)
} # }
```
