# Exons, introns, CDS or spans of genes, from a GFF

The intervals a target design starts from. Most of a gene's troublesome
tandem repeats sit in its introns, so a design will often want the exons
alone; the rest are in the exons, so the flow is gene, then exons, then
[`bed_subtract()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/bed_subtract.md)
with the repeats from
[`tandem_repeats_to_avoid()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/tandem_repeats_to_avoid.md).
This reads a GFF3 and returns the pieces you ask for as interval rows
that drop straight into
[`bed_subtract()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/bed_subtract.md),
[`bed_intersect()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/bed_intersect.md)
and
[`write_bed()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/write_bed.md).

## Usage

``` r
read_gff_features(
  gff,
  feature = c("exon", "intron", "cds", "gene"),
  ids = NULL,
  per = c("gene", "transcript"),
  genes = PF3D7_GENES
)
```

## Arguments

- gff:

  Path or URL to a GFF3 file, plain or gzipped, as for
  [`read_gff_cds()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/read_gff_cds.md).

- feature:

  One of `"exon"`, `"intron"`, `"cds"`, `"gene"`.

- ids:

  Which genes: gene ids (`"PF3D7_0709000"`), the GFF's `Name` attribute
  (`"CRT"`, case-insensitive), a `name` from `genes` (`"pfcrt"`), or
  transcript ids, which restrict that gene to the transcript named.
  `NULL` (the default) returns every gene.

- per:

  `"gene"` merges isoforms; `"transcript"` keeps them apart.

- genes:

  A gene table with `name` and `gene_id` columns, so `ids` can use the
  package's friendly names. Defaults to
  [PF3D7_GENES](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF3D7_GENES.md);
  `NULL` to use the GFF's `Name` only.

## Value

A tibble with one row per interval: `gene_id`, `name` (the `genes` name,
else the GFF `Name`, else `NA`), `transcript_id` when
`per = "transcript"`, `feature`, `index`, `chrom` (as the GFF spells
it), `chr` (normalised), `start`, `end`, `width`, `strand`. Genes in GFF
order, intervals by position.

## Details

`feature` is what to return:

- `"exon"` – the exons (the GFF's `exon` features; a GFF with none falls
  back to its `CDS` features, with a message);

- `"intron"` – the gaps between consecutive exons;

- `"cds"` – the coding parts of the exons only;

- `"gene"` – each gene's whole span, one row per gene.

`per = "gene"` (the default) gives one set of intervals per gene,
merging the exons of every isoform, so an intron is a stretch no isoform
transcribes; `per = "transcript"` keeps each transcript's own exons and
introns. `index` numbers them in transcript orientation, so exon 1 of a
minus-strand gene is the one with the highest coordinates.

## Coordinates

**0-based half-open**, like every interval table in the package (see
[plasgenomicsutilsR-coordinates](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plasgenomicsutilsR-coordinates.md))
– converted from the GFF's 1-based inclusive coordinates here, at the
boundary.
[`read_gff_cds()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/read_gff_cds.md)
is the one reader that keeps the GFF's numbering, because
[`aa_intervals()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/aa_intervals.md)
walks codons on it; do not mix the two.

## See also

[`tandem_repeats_to_avoid()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/tandem_repeats_to_avoid.md)
and
[`bed_subtract()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/bed_subtract.md)
for what comes next,
[`read_gff_cds()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/read_gff_cds.md)
for the codon-walking reader.

## Examples

``` r
gff <- system.file("extdata", "pf3d7_drug_gene_cds.gff", package = "plasgenomicsutilsR")
read_gff_features(gff, "exon", ids = "pfcrt")
#> # A tibble: 13 × 10
#>    gene_id       name  feature index chrom      chr    start    end width strand
#>    <chr>         <chr> <chr>   <int> <chr>      <chr>  <dbl>  <dbl> <dbl> <chr> 
#>  1 PF3D7_0709000 pfcrt exon        1 Pf3D7_07_… 7     402384 403312   928 +     
#>  2 PF3D7_0709000 pfcrt exon        2 Pf3D7_07_… 7     403489 403758   269 +     
#>  3 PF3D7_0709000 pfcrt exon        3 Pf3D7_07_… 7     403937 404110   173 +     
#>  4 PF3D7_0709000 pfcrt exon        4 Pf3D7_07_… 7     404282 404415   133 +     
#>  5 PF3D7_0709000 pfcrt exon        5 Pf3D7_07_… 7     404568 404640    72 +     
#>  6 PF3D7_0709000 pfcrt exon        6 Pf3D7_07_… 7     404763 404839    76 +     
#>  7 PF3D7_0709000 pfcrt exon        7 Pf3D7_07_… 7     404935 405018    83 +     
#>  8 PF3D7_0709000 pfcrt exon        8 Pf3D7_07_… 7     405145 405196    51 +     
#>  9 PF3D7_0709000 pfcrt exon        9 Pf3D7_07_… 7     405333 405390    57 +     
#> 10 PF3D7_0709000 pfcrt exon       10 Pf3D7_07_… 7     405538 405631    93 +     
#> 11 PF3D7_0709000 pfcrt exon       11 Pf3D7_07_… 7     405824 405869    45 +     
#> 12 PF3D7_0709000 pfcrt exon       12 Pf3D7_07_… 7     406016 406071    55 +     
#> 13 PF3D7_0709000 pfcrt exon       13 Pf3D7_07_… 7     406240 406341   101 +     
read_gff_features(gff, "intron", ids = "pfcrt")[, c("index", "start", "end", "width")]
#> # A tibble: 12 × 4
#>    index  start    end width
#>    <int>  <dbl>  <dbl> <dbl>
#>  1     1 403312 403489   177
#>  2     2 403758 403937   179
#>  3     3 404110 404282   172
#>  4     4 404415 404568   153
#>  5     5 404640 404763   123
#>  6     6 404839 404935    96
#>  7     7 405018 405145   127
#>  8     8 405196 405333   137
#>  9     9 405390 405538   148
#> 10    10 405631 405824   193
#> 11    11 405869 406016   147
#> 12    12 406071 406240   169
read_gff_features(gff, "gene")
#> # A tibble: 6 × 10
#>   gene_id       name      feature index chrom   chr    start    end width strand
#>   <chr>         <chr>     <chr>   <int> <chr>   <chr>  <dbl>  <dbl> <dbl> <chr> 
#> 1 PF3D7_0417200 pfdhfr    gene        1 Pf3D7_… 4     7.48e5 7.50e5  2169 +     
#> 2 PF3D7_0523000 pfmdr1    gene        1 Pf3D7_… 5     9.56e5 9.63e5  7141 +     
#> 3 PF3D7_0629500 pfaat1    gene        1 Pf3D7_… 6     1.21e6 1.22e6  4212 -     
#> 4 PF3D7_0709000 pfcrt     gene        1 Pf3D7_… 7     4.02e5 4.06e5  3957 +     
#> 5 PF3D7_0810800 pfdhps    gene        1 Pf3D7_… 8     5.48e5 5.51e5  3162 +     
#> 6 PF3D7_1343700 pfkelch13 gene        1 Pf3D7_… 13    1.72e6 1.73e6  3278 -     

# the exons of two genes, minus the tandem repeats and 10 bp of clearance around each
exons <- read_gff_features(gff, "exon", ids = c("pfcrt", "pfkelch13"))
targets <- bed_subtract(exons, tandem_repeats_to_avoid(pf3d7_tandem_repeats()), pad = 10)
targets[, c("name", "index", "start", "end", "piece", "width")]
#> # A tibble: 25 × 6
#>    name  index  start    end piece width
#>    <chr> <int>  <dbl>  <dbl> <int> <dbl>
#>  1 pfcrt     1 402393 402401     1     8
#>  2 pfcrt     1 402478 402592     2   114
#>  3 pfcrt     1 402699 402878     3   179
#>  4 pfcrt     1 402927 402942     4    15
#>  5 pfcrt     1 403026 403312     5   286
#>  6 pfcrt     2 403489 403758     1   269
#>  7 pfcrt     3 403937 404110     1   173
#>  8 pfcrt     4 404282 404415     1   133
#>  9 pfcrt     5 404568 404640     1    72
#> 10 pfcrt     6 404763 404839     1    76
#> # ℹ 15 more rows
```
