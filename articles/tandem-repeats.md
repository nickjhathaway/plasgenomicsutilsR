# Avoiding tandem repeats

``` r

library(plasgenomicsutilsR)
```

A polymerase slips on a short tandem repeat, and the shorter the
repeating unit the shorter a run it takes: a dozen `A`s is enough, a
dozen bases of `AT` too, while a trinucleotide has to run to about
twenty before it matters and anything at all is trouble from fifty.
Designing targets – amplicons, probes, windows to call in – means
knowing where those runs are and keeping a few bases clear of them. This
article goes from a repeat finder’s output to a BED of targets with the
repeats cut out, two ways: from whole genes, and from their exons.

## The repeats

[`pf3d7_tandem_repeats()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pf3d7_tandem_repeats.md)
is every short tandem repeat in the Pf3D7 reference, found once by
Tandem Repeats Finder and an exhaustive search for perfect repeats of
units up to 10 bp
([`?pf3d7_tandem_repeats`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pf3d7_tandem_repeats.md)
has the exact commands), and bundled so it never has to be run again.
[`tandem_repeats()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/tandem_repeats.md)
reads the same kind of BED for any other genome: only the first four
fields are used, and the name has to carry the unit and copy number as
`UNIT_xCOPIES`, with or without a `chrom-start-end__` in front.

``` r

tr <- pf3d7_tandem_repeats()
nrow(tr)
#> [1] 271936
tr[1:4, c("chrom", "start", "end", "width", "repeat_unit", "period", "copies")]
#> # A tibble: 4 × 7
#>   chrom       start   end width repeat_unit period copies
#>   <chr>       <dbl> <dbl> <dbl> <chr>        <int>  <dbl>
#> 1 Pf3D7_01_v3     2   356   354 AACCCTA          7  50.6 
#> 2 Pf3D7_01_v3    12    36    24 CCTAAAC          7   3.43
#> 3 Pf3D7_01_v3    51    71    20 AACCCTA          7   2.86
#> 4 Pf3D7_01_v3    72    92    20 AACCCTA          7   2.86
```

`period` is the length of the *true* repeating unit, whatever spelling
the finder chose – `AT`, `ATAT` and `ATATAT` all describe a dinucleotide
run and all get `period = 2`. That is the number the slippage rules key
on:

``` r

table(pmin(tr$period, 5))     # 1..4 and "5" for anything longer
#> 
#>      1      2      3      4      5 
#>  37476  56813  16337  15388 145922
```

## What to avoid

### Judging each repeat

[`flag_tandem_repeats()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/flag_tandem_repeats.md)
applies the thresholds: a run is flagged at 11 bp for a homopolymer, 12
for a dinucleotide, 21 for a trinucleotide, and 50 for any unit. Those
are the values the HEOME Pf3D7 design used, and each is an argument –
`min_width_by_period` is named by period, so
`c("1" = 21, "2" = 21, "3" = 21)` holds homopolymers and dinucleotides
to the trinucleotide length, and `NULL` drops the per-period rules to
flag on total length alone.

``` r

flagged <- flag_tandem_repeats(tr)
table(period = pmin(flagged$period, 4), avoid = flagged$avoid)
#>       avoid
#> period  FALSE   TRUE
#>      1      0  37476
#>      2      0  56813
#>      3  10892   5445
#>      4 133645  27665
```

### Merging runs that flow into one another

Repeats flow into one another – an `A` homopolymer into an `AT` run, one
trinucleotide into the next – and the finders report each as its own
record.
[`merge_tandem_repeats()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/merge_tandem_repeats.md)
joins records that overlap, abut, or lie within `gap` bases, and says
what went into each block:

``` r

blocks <- merge_tandem_repeats(tr, gap = 10)
blocks[blocks$n_repeats >= 4, c("name", "width", "n_repeats", "repeat_units", "periods")][1:4, ]
#> # A tibble: 4 × 5
#>   name                  width n_repeats repeat_units                     periods
#>   <chr>                 <dbl>     <int> <chr>                            <chr>  
#> 1 Pf3D7_01_v3-2-356       354        15 AACCCTA,CCTAAAC,AACCCTG,AACCCT,… 6,7    
#> 2 Pf3D7_01_v3-370-1183    813        63 TCTTCTTA,TCTTACTTACTTACTCTTATCT… 4,6,8,…
#> 3 Pf3D7_01_v3-1311-1373    62         5 TCTTATCTTCTTACTTCTCATTACTTAC,TC… 4,6,13…
#> 4 Pf3D7_01_v3-2325-2913   588         5 TAACTCATTTATTACTACTATCATCACTTAT… 3,6,268
```

The two finders often report one run twice under different spellings
(`ATA` here, `TAT` there). `dedupe = TRUE` drops a record lying inside
another with the same unit up to rotation or reverse complement; the
blocks come out the same, and the unit list stops naming both spellings.

### Judging the merged block

[`tandem_repeats_to_avoid()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/tandem_repeats_to_avoid.md)
does merging and flagging together, and `rule` says which width the
thresholds see:

- `rule = "record"` (the default) judges each repeat on its own width,
  then keeps a merged block when any repeat in it was flagged, or when
  the block as a whole reaches 50 bp – a chain of individually harmless
  short repeats slips like one long one. This is the HEOME mask.
- `rule = "block"` merges first and judges each block on its
  **combined** width, held to the lowest threshold among the periods it
  contains. Any block with a homopolymer in it is held to 11 bp,
  whatever else is in it.

The difference is visible on two short repeats near each other. An 8 bp
`A` run and a 9 bp `ATA` run 3 bp apart are each under their threshold,
so the record rule keeps them; merged within 10 bp they are one 20 bp
stretch holding a homopolymer, and the block rule avoids it:

``` r

close <- data.frame(chrom = "Pf3D7_07_v3", start = c(100, 111), end = c(108, 120),
                    name = c("A_x8", "ATA_x3"))
tandem_repeats_to_avoid(close, gap = 10)                 # nothing
#> # A tibble: 0 × 10
#> # ℹ 10 variables: chrom <chr>, chr <chr>, start <dbl>, end <dbl>, width <dbl>,
#> #   name <chr>, n_repeats <int>, repeat_units <chr>, periods <chr>,
#> #   min_period <int>
tandem_repeats_to_avoid(close, gap = 10, rule = "block")[, c("start", "end", "width", "periods")]
#> # A tibble: 1 × 4
#>   start   end width periods
#>   <dbl> <dbl> <dbl> <chr>  
#> 1   100   120    20 1,3
```

The masks these give for the genome:

``` r

avoid <- tandem_repeats_to_avoid(tr)
avoid_block <- tandem_repeats_to_avoid(tr, gap = 10, rule = "block")
data.frame(rule = c("record", "block, gap 10"),
           blocks = c(nrow(avoid), nrow(avoid_block)),
           bases = c(sum(avoid$width), sum(avoid_block$width)))
#>            rule blocks   bases
#> 1        record  75931 4156065
#> 2 block, gap 10  72685 4423471
```

The rest of this article uses `avoid`, the HEOME mask. Swap in
`avoid_block`, or a stricter `min_width_by_period`, and everything below
runs the same.

## Workflow 1: genes, minus the repeats, to a BED

Take the genes you are after from `PF3D7_GENES`, subtract the mask with
[`bed_subtract()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/bed_subtract.md),
and pad it so nothing you design ends right at a repeat’s edge. What
comes back is the pieces of each gene that are clear of repeats, still
carrying the gene’s own columns:

``` r

genes <- PF3D7_GENES[PF3D7_GENES$name %in% c("pfcrt", "pfdhfr", "pfkelch13"), ]
targets <- bed_subtract(genes, avoid, pad = 10)
targets[, c("name", "Pf3D7_chrom", "start", "end", "piece", "width")]
#> # A tibble: 18 × 6
#>    name      Pf3D7_chrom   start     end piece width
#>    <chr>     <chr>         <dbl>   <dbl> <int> <dbl>
#>  1 pfdhfr    Pf3D7_04_v3  748087  749914     1  1827
#>  2 pfcrt     Pf3D7_07_v3  403221  403363     1   142
#>  3 pfcrt     Pf3D7_07_v3  403480  403823     2   343
#>  4 pfcrt     Pf3D7_07_v3  403899  404123     3   224
#>  5 pfcrt     Pf3D7_07_v3  404281  404418     4   137
#>  6 pfcrt     Pf3D7_07_v3  404478  404524     5    46
#>  7 pfcrt     Pf3D7_07_v3  404561  404692     6   131
#>  8 pfcrt     Pf3D7_07_v3  404725  404877     7   152
#>  9 pfcrt     Pf3D7_07_v3  404929  405103     8   174
#> 10 pfcrt     Pf3D7_07_v3  405146  405263     9   117
#> 11 pfcrt     Pf3D7_07_v3  405342  405454    10   112
#> 12 pfcrt     Pf3D7_07_v3  405510  405635    11   125
#> 13 pfcrt     Pf3D7_07_v3  405736  405883    12   147
#> 14 pfcrt     Pf3D7_07_v3  405921  405961    13    40
#> 15 pfcrt     Pf3D7_07_v3  406022  406129    14   107
#> 16 pfcrt     Pf3D7_07_v3  406175  406317    15   142
#> 17 pfkelch13 Pf3D7_13_v3 1724816 1726561     1  1745
#> 18 pfkelch13 Pf3D7_13_v3 1726617 1726997     2   380
```

*pfcrt* comes back in many pieces, which is the point: it is a
repeat-rich gene, and a design that treats it as one interval will keep
landing primers on the runs. `min_width` drops pieces too short to be
useful.

[`write_bed6()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/write_bed6.md)
writes the six-column BED – chromosome, start, end, name, score, strand
– using the reference’s own chromosome names (`Pf3D7_chrom`, not the
short `chrom` key) so the file matches the FASTA. The score is the width
unless the table has a `score` column, and `meta_cols =` carries any
other columns in a seventh `[field=value;]` column:

``` r

bed <- file.path(tempdir(), "genes_no_repeats.bed")
write_bed6(targets, bed, meta_cols = c("gene_id", "piece"))
readLines(bed)[1:3]
#> [1] "Pf3D7_04_v3\t748087\t749914\tpfdhfr\t1827\t+\t[gene_id=PF3D7_0417200;piece=1;]"
#> [2] "Pf3D7_07_v3\t403221\t403363\tpfcrt\t142\t+\t[gene_id=PF3D7_0709000;piece=1;]"  
#> [3] "Pf3D7_07_v3\t403480\t403823\tpfcrt\t343\t+\t[gene_id=PF3D7_0709000;piece=2;]"
```

Every piece of *pfcrt* is named `pfcrt` there. Tools that key on the
name want it unique within the file: `make_names_unique = TRUE` numbers
repeats of a name in file order, zero-padded to the count (`_01` to
`_15`), and `name_is_coords = TRUE` names each interval by its
coordinates instead. The two combine, so identical coordinates written
twice still get distinct names:

``` r

write_bed6(targets, bed, make_names_unique = TRUE)
readLines(bed)[2:3]
#> [1] "Pf3D7_07_v3\t403221\t403363\tpfcrt_01\t142\t+"
#> [2] "Pf3D7_07_v3\t403480\t403823\tpfcrt_02\t343\t+"
write_bed6(targets, bed, name_is_coords = TRUE)
readLines(bed)[2]
#> [1] "Pf3D7_07_v3\t403221\t403363\tPf3D7_07_v3-403221-403363\t142\t+"
```

[`write_bed()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/write_bed.md)
is the three- or four-column form, for `bedtools` and `bcftools -R`.

## Workflow 2: genes, to exons, minus the repeats

Most of a gene’s worst repeats sit in its introns, so a design will
often want the exons alone – and the repeats that remain are in the
exons, so the subtraction still follows.
[`read_gff_features()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/read_gff_features.md)
reads exons (or introns, CDS, or gene spans) out of a GFF3 as 0-based
intervals.
[`ensembl_gff_url()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ensembl_gff_url.md)
points it at a whole released annotation:

``` r

gff <- ensembl_gff_url("falciparum")         # Ensembl Protists, read straight from the web
exons <- read_gff_features(gff, "exon", ids = c("pfcrt", "pfdhfr", "pfkelch13"))
```

The package ships the annotation of six drug-resistance genes, verbatim
from PlasmoDB, which is what runs here; the coordinates are the same:

``` r

gff <- system.file("extdata", "pf3d7_drug_gene_cds.gff", package = "plasgenomicsutilsR")
exons <- read_gff_features(gff, "exon", ids = c("pfcrt", "pfdhfr", "pfkelch13"))
exons[exons$name == "pfcrt", c("index", "start", "end", "width", "strand")]
#> # A tibble: 13 × 5
#>    index  start    end width strand
#>    <int>  <dbl>  <dbl> <dbl> <chr> 
#>  1     1 402384 403312   928 +     
#>  2     2 403489 403758   269 +     
#>  3     3 403937 404110   173 +     
#>  4     4 404282 404415   133 +     
#>  5     5 404568 404640    72 +     
#>  6     6 404763 404839    76 +     
#>  7     7 404935 405018    83 +     
#>  8     8 405145 405196    51 +     
#>  9     9 405333 405390    57 +     
#> 10    10 405538 405631    93 +     
#> 11    11 405824 405869    45 +     
#> 12    12 406016 406071    55 +     
#> 13    13 406240 406341   101 +
```

`ids` takes gene ids, the GFF’s own `Name` (`"CRT"`), the package’s
names (`"pfcrt"`, via `genes = PF3D7_GENES`), or transcript ids. `index`
numbers exons in transcript orientation, so exon 1 of the minus-strand
*pfkelch13* is the one with the highest coordinates;
`per = "transcript"` keeps isoforms apart instead of merging them.

Then the same subtraction, and the same BED:

``` r

exon_targets <- bed_subtract(exons, avoid, pad = 10)
exon_targets[exon_targets$name == "pfcrt", c("index", "start", "end", "piece", "width")]
#> # A tibble: 17 × 5
#>    index  start    end piece width
#>    <int>  <dbl>  <dbl> <int> <dbl>
#>  1     1 402393 402401     1     8
#>  2     1 402478 402592     2   114
#>  3     1 402699 402878     3   179
#>  4     1 402927 402942     4    15
#>  5     1 403026 403312     5   286
#>  6     2 403489 403758     1   269
#>  7     3 403937 404110     1   173
#>  8     4 404282 404415     1   133
#>  9     5 404568 404640     1    72
#> 10     6 404763 404839     1    76
#> 11     7 404935 405018     1    83
#> 12     8 405146 405196     1    50
#> 13     9 405342 405390     1    48
#> 14    10 405538 405631     1    93
#> 15    11 405824 405869     1    45
#> 16    12 406022 406071     1    49
#> 17    13 406240 406317     1    77

bed <- file.path(tempdir(), "exons_no_repeats.bed")
write_bed6(exon_targets, bed, meta_cols = c("gene_id", "index", "piece"))
readLines(bed)[1:3]
#> [1] "Pf3D7_04_v3\t747918\t748013\tpfdhfr\t95\t+\t[gene_id=PF3D7_0417200;index=1;piece=1;]"  
#> [2] "Pf3D7_04_v3\t748087\t750009\tpfdhfr\t1922\t+\t[gene_id=PF3D7_0417200;index=1;piece=2;]"
#> [3] "Pf3D7_07_v3\t402393\t402401\tpfcrt\t8\t+\t[gene_id=PF3D7_0709000;index=1;piece=1;]"
```

How much each route keeps:

``` r

data.frame(route = c("gene minus repeats", "exons minus repeats"),
           pieces = c(nrow(targets), nrow(exon_targets)),
           bases = c(sum(targets$width), sum(exon_targets$width)))
#>                 route pieces bases
#> 1  gene minus repeats     18  6091
#> 2 exons minus repeats     27  6517
```

`feature = "intron"` gives the gaps between exons if the introns are
what you want to look at, and `feature = "cds"` the coding parts only.
