# Genomic interval of an amino-acid position

Turns "codon 76 of *pfcrt*" into a genomic interval, by walking the
protein back through the transcript's coding exons. Amino-acid positions
are how resistance markers are named and reported, while every plot and
interval tool here works in genomic coordinates, and the conversion is
not something you can do by eye: it depends on the exon structure, the
strand and the CDS phase.

## Usage

``` r
aa_intervals(positions, gff, genes = PF3D7_GENES, one_based_output = FALSE)
```

## Arguments

- positions:

  A data frame with `transcript_id` and `aa_position`. `aa_position` is
  **1-based**, counting the initiator methionine as 1, matching how
  residues are numbered in the literature – so `76` is the residue
  everyone calls 76. `transcript_id` may be a transcript id
  (`"PF3D7_0709000.1"`), a gene id (`"PF3D7_0709000"` – every transcript
  of it is returned), or a gene symbol when `genes` is given
  (`"pfcrt"`).

- gff:

  A GFF path, or the result of
  [`read_gff_cds()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/read_gff_cds.md)
  (parse once, reuse).

- genes:

  Optional gene table with `name` and `gene_id` columns, so
  `transcript_id` can be a symbol; defaults to
  [PF3D7_GENES](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF3D7_GENES.md).

- one_based_output:

  Return 1-based inclusive coordinates instead of the package's 0-based
  half-open convention. `FALSE` (default) keeps `start` 0-based so the
  output can be used as an interval table directly.

## Value

A tibble with `chr` (normalised), `chrom` (as the GFF spells it),
`start`, `end`, `name`, `transcript_id`, `gene_id`, `aa_position`,
`strand`, `codon_positions` (a comma-separated list of the three base
positions, always 1-based as coordinates are usually quoted), `n_exons`
and `spans_intron`. Positions past the end of a protein, and ids that
match nothing, are dropped with a warning naming them.

## Details

The result is shaped like the package's other interval tables – `chr`,
`start`, `end`, `name` – so it drops straight into `genes =`,
`mark_snps =`,
[`annotate_snps()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/annotate_snps.md)
or
[`bed_intersect()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/bed_intersect.md).
`name` is `<transcript_id>-AA<aa_position>`. Those genomic bounds are
0-based half-open like every other interval here, while the amino-acid
position itself is 1-based and `codon_positions` lists 1-based bases –
the deliberate exception described in
[`?"plasgenomicsutilsR-coordinates"`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plasgenomicsutilsR-coordinates.md).

A codon can straddle an intron, in which case its three bases are not
contiguous: `start` and `end` then span the intron as well, and
`spans_intron` flags it so the width is not mistaken for three bases.
`codon_positions` always lists the three base positions themselves.

## See also

[`read_gff_cds()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/read_gff_cds.md),
[`annotate_snps()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/annotate_snps.md)

## Examples

``` r
if (FALSE) { # \dontrun{
cds <- read_gff_cds("Pf3D7.gff")
aa <- data.frame(transcript_id = c("pfcrt", "pfkelch13"), aa_position = c(76, 580))
aa_intervals(aa, cds)
} # }
```
