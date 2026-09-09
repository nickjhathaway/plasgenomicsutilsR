# The amino acids a genomic range covers

For each interval, every transcript whose CDS it overlaps and the range
of codons of that transcript the interval covers – the residues a target
window, an amplicon, or a piece of gene left after
[`bed_subtract()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/bed_subtract.md)
will actually read. Covered means the codon's bases lie inside the
interval, so a window that starts 100 bp upstream of a gene reports from
codon 1, and one that starts inside the gene reports from the first
codon it holds whole; nothing is inferred from the interval's first and
last positions.

## Usage

``` r
genomic_range_aa_positions(
  ranges,
  gff,
  keep = c("all", "hits"),
  partial = FALSE,
  fasta = NULL,
  genes = PF3D7_GENES,
  chrom = "chr",
  start = "start",
  end = "end"
)
```

## Arguments

- ranges:

  An interval table (`chr` or `chrom`, `start`, `end`), such as
  [`bed_subtract()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/bed_subtract.md)
  or
  [`read_gff_features()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/read_gff_features.md)
  return, or a gene table.

- gff:

  Path or URL to a GFF3, or the table
  [`read_gff_cds()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/read_gff_cds.md)
  returns – parse once and reuse it.

- keep:

  `"all"` (the default) keeps every input row, filling the result
  columns with `NA` for an interval that covers no codon; `"hits"` drops
  those rows.

- partial:

  Count a codon the interval covers only partly (default `FALSE`).

- fasta:

  Optional sequence, as for
  [`snp_aa_positions()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/snp_aa_positions.md):
  a path or URL to the genome (gzipped is fine), or a named vector of
  sequences. When given, `aa_seq` holds the reference residues of the
  covered codons, `aa_start` to `aa_end`, in transcript order; a stop is
  `*` and an unreadable codon `X`. A GFF that ends in a `##FASTA`
  section supplies its own sequence.

- genes:

  A gene table with `name` and `gene_id`, used to fill `gene_name`.
  Defaults to
  [PF3D7_GENES](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF3D7_GENES.md);
  `NULL` leaves `gene_name` `NA`.

- chrom, start, end:

  Column names in `ranges` (defaults `"chr"`, `"start"`, `"end"`;
  `"chrom"` is accepted as a chromosome alias).

## Value

A tibble with one row per interval and covered transcript: every column
of `ranges`, then `transcript_id`, `gene_id`, `gene_name`, `aa_start`,
`aa_end`, `n_aa`, `partial_start`, `partial_end`, `strand`, `coding`
(whether any codon was covered), and `aa_seq` when sequence is
available. Rows follow the input order; an interval covering several
transcripts has several rows.

## Details

A codon is covered when all three of its bases are inside the interval.
`partial = TRUE` also counts a codon the interval holds only one or two
bases of – the codon it starts or ends in – and `partial_start` /
`partial_end` say when that happened. Because an interval is contiguous
and a transcript's coding bases are in order along it, the covered
codons of one transcript are always a contiguous run, `aa_start` to
`aa_end`.

Codon numbers are **1-based**, counting the initiator methionine as 1,
as everywhere in the package
([plasgenomicsutilsR-coordinates](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plasgenomicsutilsR-coordinates.md)).
A CDS that includes its stop codon reports it as the last codon (425 for
*pfcrt*, whose protein is 424 residues), as
[`snp_aa_positions()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/snp_aa_positions.md)
does; in `aa_seq` it is `*`. Genomic coordinates in `ranges` are 0-based
half-open, like every interval table.

## See also

[`snp_aa_positions()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/snp_aa_positions.md)
for single positions,
[`aa_intervals()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/aa_intervals.md)
for the other direction,
[`bed_subtract()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/bed_subtract.md)
and
[`read_gff_features()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/read_gff_features.md)
for where the intervals come from.

## Examples

``` r
gff <- system.file("extdata", "pf3d7_drug_gene_cds.gff", package = "plasgenomicsutilsR")
cds <- read_gff_cds(gff)

# pfcrt with its tandem repeats removed: which residues does each piece still read?
crt <- PF3D7_GENES[PF3D7_GENES$name == "pfcrt", ]
pieces <- bed_subtract(crt, tandem_repeats_to_avoid(pf3d7_tandem_repeats()), pad = 10)
genomic_range_aa_positions(pieces, cds)[, c("piece", "start", "end", "aa_start", "aa_end",
                                            "n_aa")]
#> # A tibble: 15 × 6
#>    piece  start    end aa_start aa_end  n_aa
#>    <int>  <dbl>  <dbl>    <int>  <int> <int>
#>  1     1 403221 403363        1     30    30
#>  2     2 403480 403823       32    120    89
#>  3     3 403899 404123      121    177    57
#>  4     4 404281 404418      179    222    44
#>  5     5 404478 404524       NA     NA    NA
#>  6     6 404561 404692      223    246    24
#>  7     7 404725 404877      247    271    25
#>  8     8 404929 405103      273    299    27
#>  9     9 405146 405263      301    316    16
#> 10    10 405342 405454      320    335    16
#> 11    11 405510 405635      336    366    31
#> 12    12 405736 405883      367    381    15
#> 13    13 405921 405961       NA     NA    NA
#> 14    14 406022 406129      384    399    16
#> 15    15 406175 406317      401    425    25

# the K76T codon itself, and a window holding only its middle base
k76 <- data.frame(chr = "Pf3D7_07_v3", start = c(403623, 403624), end = c(403626, 403625))
genomic_range_aa_positions(k76, cds)[, c("start", "end", "aa_start", "aa_end", "coding")]
#> # A tibble: 2 × 5
#>    start    end aa_start aa_end coding
#>    <dbl>  <dbl>    <int>  <int> <lgl> 
#> 1 403623 403626       76     76 TRUE  
#> 2 403624 403625       NA     NA FALSE 
genomic_range_aa_positions(k76, cds, partial = TRUE)[, c("start", "end", "aa_start",
                                                         "partial_start", "partial_end")]
#> # A tibble: 2 × 5
#>    start    end aa_start partial_start partial_end
#>    <dbl>  <dbl>    <int> <lgl>         <lgl>      
#> 1 403623 403626       76 FALSE         FALSE      
#> 2 403624 403625       76 TRUE          TRUE       

# with sequence, the residues themselves
fa <- system.file("extdata", "pf3d7_drug_gene_regions.fasta.gz", package = "plasgenomicsutilsR")
genomic_range_aa_positions(k76[1, ], cds, fasta = fa)$aa_seq
#> [1] "K"
```
