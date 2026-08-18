# Amino acids and genomic coordinates

Resistance markers are named by residue — *pfcrt* K76T, *pfdhps* A437G,
*pfkelch13* C580Y — while callsets, scans and every plot in this package
work in genomic coordinates. Going between the two is not arithmetic you
can do in your head: it depends on the transcript’s exon structure,
which strand it is on, and the CDS phase.

Two functions do it, in opposite directions:

- [`aa_intervals()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/aa_intervals.md)
  — a residue → the genomic interval of its codon
- [`snp_aa_positions()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/snp_aa_positions.md)
  — a SNP → the residue it sits on, or `NA` if it is not coding

Both read coding exons from a GFF via
[`read_gff_cds()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/read_gff_cds.md).

``` r

library(plasgenomicsutilsR)
```

## Getting an annotation

[`read_gff_cds()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/read_gff_cds.md)
takes a path or a URL, and will read a gzipped file straight from the
web, so a released annotation can be used without keeping a copy:

``` r

# VEuPathDB / PlasmoDB — the source the bundled gene datasets were built from
cds <- read_gff_cds(paste0("https://plasmodb.org/common/downloads/Current_Release/",
                           "Pfalciparum3D7/gff/data/PlasmoDB-68_Pfalciparum3D7.gff"))

# Ensembl Protists, gzipped
cds <- read_gff_cds(paste0("https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists/current/",
                           "gff3/plasmodium_falciparum/",
                           "Plasmodium_falciparum.GCA000002765v3.63.gff3.gz"))
```

The examples here use a small file that ships with the package instead:
the real CDS features of six 3D7 drug-resistance genes, taken verbatim
from PlasmoDB. Real coordinates, so every number below is the published
one.

``` r

cds <- read_gff_cds(system.file("extdata", "pf3d7_drug_gene_cds.gff",
                                package = "plasgenomicsutilsR"))
cds
#> # A tibble: 21 × 7
#>    transcript_id   gene_id       chrom         start     end strand phase
#>    <chr>           <chr>         <chr>         <int>   <int> <chr>  <int>
#>  1 PF3D7_0417200.1 PF3D7_0417200 Pf3D7_04_v3  748088  749914 +          0
#>  2 PF3D7_0523000.1 PF3D7_0523000 Pf3D7_05_v3  957890  962149 +          0
#>  3 PF3D7_0629500.1 PF3D7_0629500 Pf3D7_06_v3 1213948 1214092 -          1
#>  4 PF3D7_0629500.1 PF3D7_0629500 Pf3D7_06_v3 1214330 1216005 -          0
#>  5 PF3D7_0709000.1 PF3D7_0709000 Pf3D7_07_v3  403222  403312 +          0
#>  6 PF3D7_0709000.1 PF3D7_0709000 Pf3D7_07_v3  403490  403758 +          2
#>  7 PF3D7_0709000.1 PF3D7_0709000 Pf3D7_07_v3  403938  404110 +          0
#>  8 PF3D7_0709000.1 PF3D7_0709000 Pf3D7_07_v3  404283  404415 +          1
#>  9 PF3D7_0709000.1 PF3D7_0709000 Pf3D7_07_v3  404569  404640 +          0
#> 10 PF3D7_0709000.1 PF3D7_0709000 Pf3D7_07_v3  404764  404839 +          0
#> # ℹ 11 more rows
```

## A residue to its codon

Six markers you would recognise from a resistance table:

``` r

markers <- data.frame(
  transcript_id = c("pfcrt", "pfdhps", "pfdhps", "pfdhps", "pfdhfr", "pfmdr1", "pfkelch13"),
  aa_position   = c(     76,      437,      540,      581,      108,       86,         580))

aa_intervals(markers, cds, one_based_output = TRUE)[
  , c("name", "chr", "start", "end", "codon_positions", "strand")]
#> # A tibble: 7 × 6
#>   name                  chr     start     end codon_positions         strand
#>   <chr>                 <chr>   <dbl>   <dbl> <chr>                   <chr> 
#> 1 PF3D7_0709000.1-AA76  7      403624  403626 403624,403625,403626    +     
#> 2 PF3D7_0810800.1-AA437 8      549684  549686 549684,549685,549686    +     
#> 3 PF3D7_0810800.1-AA540 8      549993  549995 549993,549994,549995    +     
#> 4 PF3D7_0810800.1-AA581 8      550116  550118 550116,550117,550118    +     
#> 5 PF3D7_0417200.1-AA108 4      748409  748411 748409,748410,748411    +     
#> 6 PF3D7_0523000.1-AA86  5      958145  958147 958145,958146,958147    +     
#> 7 PF3D7_1343700.1-AA580 13    1725258 1725260 1725258,1725259,1725260 -
```

*pfcrt* codon 76 lands at 403,624–403,626 and *pfkelch13* codon 580 at
1,725,258–1,725,260 — the positions these mutations are reported at.

`transcript_id` took gene symbols above because `genes` defaults to
`PF3D7_GENES`. A transcript id (`"PF3D7_0709000.1"`) or a gene id
(`"PF3D7_0709000"`, returning every transcript) works the same way.

Coordinates come back **0-based half-open** by default, matching every
other interval in the package
([`?"plasgenomicsutilsR-coordinates"`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plasgenomicsutilsR-coordinates.md)),
which is what lets the result be used as an interval table directly.
`one_based_output = TRUE`, used above, gives the 1-based numbers you
would quote in text.

### When a codon straddles an intron

*pfcrt* has 13 coding exons, and four of its codons are split across two
of them:

``` r

pfcrt_all <- aa_intervals(data.frame(transcript_id = "pfcrt", aa_position = 1:424), cds,
                          one_based_output = TRUE)
subset(pfcrt_all, spans_intron)[, c("name", "start", "end", "codon_positions")]
#> # A tibble: 4 × 4
#>   name                   start    end codon_positions     
#>   <chr>                  <dbl>  <dbl> <chr>               
#> 1 PF3D7_0709000.1-AA31  403312 403491 403312,403490,403491
#> 2 PF3D7_0709000.1-AA178 404109 404283 404109,404110,404283
#> 3 PF3D7_0709000.1-AA272 404839 404937 404839,404936,404937
#> 4 PF3D7_0709000.1-AA400 406071 406242 406071,406241,406242
```

Codon 178 is bases 404,109 and 404,110 at the end of one exon, then
404,283 at the start of the next. `start` and `end` span the intron with
them, so that interval is 175 bp wide rather than 3 — `spans_intron`
flags it, and `codon_positions` always lists the three bases themselves.

### On the minus strand

*pfkelch13* is transcribed right to left, so its transcript is read from
the highest coordinate down. Codon 580’s first base is its *highest*
genomic position:

``` r

snp_aa_positions(data.frame(chr = "Pf3D7_13_v3", pos = c(1725260, 1725259, 1725258)),
                 cds, keep = "hits", one_based_snps = TRUE)[
  , c("pos", "aa_position", "codon_base", "strand")]
#> # A tibble: 3 × 4
#>       pos aa_position codon_base strand
#>     <dbl>       <int>      <int> <chr> 
#> 1 1725260         580          1 -     
#> 2 1725259         580          2 -     
#> 3 1725258         580          3 -
```

## A SNP to its residue

The other direction, on positions around *pfdhps* codon 437:

``` r

snp_aa_positions(data.frame(chr = "Pf3D7_08_v3",
                            pos = c(549684, 549685, 549686, 548400, 550700)),
                 cds, one_based_snps = TRUE)[
  , c("pos", "gene_id", "aa_position", "codon_base", "coding")]
#> # A tibble: 5 × 5
#>      pos gene_id       aa_position codon_base coding
#>    <dbl> <chr>               <int>      <int> <lgl> 
#> 1 549684 PF3D7_0810800         437          1 TRUE  
#> 2 549685 PF3D7_0810800         437          2 TRUE  
#> 3 549686 PF3D7_0810800         437          3 TRUE  
#> 4 548400 NA                     NA         NA FALSE 
#> 5 550700 NA                     NA         NA FALSE
```

`aa_position` is **1-based**, counting the initiator methionine as 1 —
the one deliberate exception to this package’s 0-based rule, at both
ends: where you supply a residue and where one is reported. Residue
numbering is too settled to renumber (*pfcrt* K76T is codon 76 wherever
it is written down, and calling it 75 would be wrong in the only sense
that matters). The genomic `start`/`end` in the same table stay 0-based
half-open; `codon_positions` is the one place a 1-based *position* is
reported, since those are the numbers people quote.

`codon_base` says which of the codon’s three bases the SNP is, in
transcript orientation. `548400` is intronic and `550700` past the CDS,
so both are non-coding and come back `NA`; `keep = "hits"` drops them
instead.

The two functions are exact inverses:

``` r

iv <- aa_intervals(markers, cds, one_based_output = TRUE)
mid <- as.integer(vapply(strsplit(iv$codon_positions, ","), `[`, character(1), 2))
back <- snp_aa_positions(data.frame(chr = iv$chrom, pos = mid), cds, keep = "hits",
                         one_based_snps = TRUE)
all(back$aa_position == markers$aa_position)
#> [1] TRUE
```

Being inside a codon is not the same as changing the residue: that needs
the alleles, which neither function has.

## The reference residue

What the reference *does* say is which residue the codon codes for
today. Hand
[`snp_aa_positions()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/snp_aa_positions.md)
some sequence and it adds `ref_codon` and `ref_aa`. `fasta =` takes the
same path-or-URL as
[`read_gff_cds()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/read_gff_cds.md),
so point it at the released genome:

``` r

cds <- read_gff_cds(paste0("https://plasmodb.org/common/downloads/Current_Release/",
                           "Pfalciparum3D7/gff/data/PlasmoDB-68_Pfalciparum3D7.gff"))
genome <- paste0("https://plasmodb.org/common/downloads/Current_Release/",
                 "Pfalciparum3D7/fasta/data/PlasmoDB-68_Pfalciparum3D7_Genome.fasta")

markers <- data.frame(
  transcript_id = c("pfcrt", "pfdhps", "pfdhps", "pfdhfr", "pfmdr1", "pfkelch13"),
  aa_position   = c(     76,      437,      581,      108,       86,         580))
codons <- aa_intervals(markers, cds, one_based_output = TRUE)
middle <- as.integer(vapply(strsplit(codons$codon_positions, ","), `[`, character(1), 2))

snp_aa_positions(data.frame(chr = codons$chrom, pos = middle), cds, keep = "hits",
                 one_based_snps = TRUE, fasta = genome)[
  , c("gene_id", "aa_position", "ref_codon", "ref_aa", "strand")]
```

    #> # A tibble: 6 × 5
    #>   gene_id       aa_position ref_codon ref_aa strand
    #>   <chr>               <int> <chr>     <chr>  <chr>
    #> 1 PF3D7_0709000          76 AAA       K      +
    #> 2 PF3D7_0810800         437 GGT       G      +
    #> 3 PF3D7_0810800         581 GCG       A      +
    #> 4 PF3D7_0417200         108 AGC       S      +
    #> 5 PF3D7_0523000          86 AAT       N      +
    #> 6 PF3D7_1343700         580 TGT       C      -

Read down `ref_aa`: K76, A581, S108, N86, C580 — the residue each of
those markers is named for. **Except *pfdhps* 437, which comes back
`G`.** That is not an error: 3D7 itself carries the 437G allele, so a
callset aligned to 3D7 shows no variant at the position everyone calls
A437G.

That is worth internalising rather than filing as a curiosity. A
reference is one isolate’s genome — not a consensus, and not the
ancestral or wild-type sequence — and across *Plasmodium* species the
reference is sometimes the non-wild-type allele. “REF” means “what this
isolate has”, never “what came first”, so a marker absent from your
callset may be absent because the reference already carries it. Checking
`ref_aa` is how you catch that before it becomes a conclusion.

Two defaults in this package follow from the same fact:
[`load_genotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/load_genotypes.md)
records which allele its dosages count instead of letting REF stand in
for a baseline, and
[`run_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_ihs.md)
/
[`plot_ehh()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ehh.md)
use `polarized = FALSE`, treating the two states as simply the two
states, because ancestral versus derived cannot be read off REF and ALT.

Reading a genome over the network costs a few seconds and it is all held
in memory, so do it once and pass the result around — `fasta =` also
takes a named vector of sequences, which is what
[`read_gff_cds()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/read_gff_cds.md)’s
own reader returns.

Sequence can also come from the GFF itself, when it ends with a
`##FASTA` section —
[`read_gff_cds()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/read_gff_cds.md)
keeps it and nothing extra needs passing:

``` r

gff <- tempfile(fileext = ".gff")
writeLines(c("##gff-version 3",
             "demo\t.\tCDS\t11\t25\t.\t+\t0\tID=c1;Parent=T.1;gene_id=T",
             "##FASTA", ">demo", "CCCCCCCCCCATGAAATTTGGGTAAC"), gff)

snp_aa_positions(data.frame(chr = "demo", pos = c(11, 14, 17, 20, 23)), read_gff_cds(gff),
                 keep = "hits", one_based_snps = TRUE)[
  , c("pos", "aa_position", "ref_codon", "ref_aa")]
#> 1 sequence(s) read from the GFF's own ##FASTA section
#> # A tibble: 5 × 4
#>     pos aa_position ref_codon ref_aa
#>   <dbl>       <int> <chr>     <chr> 
#> 1    11           1 ATG       M     
#> 2    14           2 AAA       K     
#> 3    17           3 TTT       F     
#> 4    20           4 GGG       G     
#> 5    23           5 TAA       *
```

Otherwise pass `fasta =` a path or URL to a genome — gzipped is fine,
and it reads straight from the web like
[`read_gff_cds()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/read_gff_cds.md)
does. Without sequence from either source the two columns are simply
absent: a reference base cannot be inferred from an annotation. Some
GFFs do carry a translated protein in a feature attribute, but too few
agree on how for it to be worth reading.

## Either attribute convention

PlasmoDB writes `Parent=PF3D7_0709000.1` alongside a `gene_id`; Ensembl
writes `Parent=transcript:PF3D7_0709000.1` and no gene id at all.
[`read_gff_cds()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/read_gff_cds.md)
normalises both, deriving the gene id from the transcript when it has
to. Rewriting the bundled records the way Ensembl writes them lands on
the same coordinates:

``` r

rows <- grep("^#", readLines(system.file("extdata", "pf3d7_drug_gene_cds.gff",
                                         package = "plasgenomicsutilsR")),
             invert = TRUE, value = TRUE)
ensembl_style <- sub("ID=[^;]*;Parent=([^;]*);gene_id=[^;]*.*", "Parent=transcript:\\1", rows)

sub(".*\t", "", rows[1])            # the attribute column as PlasmoDB writes it
#> [1] "ID=PF3D7_0810800.1-p1-CDS1;Parent=PF3D7_0810800.1;gene_id=PF3D7_0810800;protein_source_id=PF3D7_0810800.1-p1"
sub(".*\t", "", ensembl_style[1])   # and as Ensembl writes it
#> [1] "Parent=transcript:PF3D7_0810800.1"
```

``` r

ens <- tempfile(fileext = ".gff")
writeLines(ensembl_style, ens)

codon76 <- function(x) aa_intervals(data.frame(transcript_id = "PF3D7_0709000.1",
                                               aa_position = 76), x, genes = NULL,
                                    one_based_output = TRUE)
both <- rbind(codon76(cds), codon76(read_gff_cds(ens)))
both$read_as <- c("PlasmoDB", "Ensembl")
both[, c("read_as", "chrom", "start", "end", "gene_id")]
#> # A tibble: 2 × 5
#>   read_as  chrom        start    end gene_id      
#>   <chr>    <chr>        <dbl>  <dbl> <chr>        
#> 1 PlasmoDB Pf3D7_07_v3 403624 403626 PF3D7_0709000
#> 2 Ensembl  Pf3D7_07_v3 403624 403626 PF3D7_0709000
```

## Watch the position base

Everything in the package is 0-based, genotype matrices included:
[`load_genotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/load_genotypes.md)
shifts SNPRelate’s 1-based `POS` on the way in, so a scan built from
those genotypes lines up exactly with the IBD tables and with these
codon intervals. Nothing extra to do.

The argument below exists for tables that carry raw VCF positions —
something you read with `bcftools query`, a marker list typed out of a
paper, or a genotype matrix from a version before that shift. Getting it
wrong is a silent one-base slide:

``` r

one_snp <- data.frame(snp_id = "Pf3D7_07_v3:403624")

# 0-based, the package convention: the second base of pfcrt codon 76
snp_aa_positions(one_snp, cds)[, c("snp_id", "aa_position", "codon_base")]
#> # A tibble: 1 × 3
#>   snp_id             aa_position codon_base
#>   <chr>                    <int>      <int>
#> 1 Pf3D7_07_v3:403624          76          2

# the same number read as a VCF POS: the first base instead
snp_aa_positions(one_snp, cds, one_based_snps = TRUE)[, c("snp_id", "aa_position",
                                                          "codon_base")]
#> # A tibble: 1 × 3
#>   snp_id             aa_position codon_base
#>   <chr>                    <int>      <int>
#> 1 Pf3D7_07_v3:403624          76          1
```

[`annotate_snps()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/annotate_snps.md)
takes the same argument, for the same reason. A `PopStructure` reports
which convention its own positions follow with `$positions()`.

## Putting it to work

Because
[`aa_intervals()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/aa_intervals.md)
returns an ordinary interval table, it goes straight into the plots and
the interval tools. A few things this makes easy:

``` r

codons <- aa_intervals(data.frame(transcript_id = "pfdhps",
                                  aa_position   = c(437, 540, 581)), cds)

# 1. which genotyped SNPs actually are those codons?
snps <- data.frame(snp_id = colnames(ps$genotype("full")))
annotate_snps(snps, codons, keep = "hits")

# 2. mark them on a haplotype plot
plot_region_haplotypes(ps, "pfdhps", split = "region", pad = 20000, mark_snps = codons)

# 3. label a zoomed panel by residue instead of by gene
plot_ibd_locus(ibd, "pfdhps", pad = 20000, genes_for_track = codons,
               gene_label_angle = 90)

# 4. and the other direction: every coding SNP in a gene, with its residue
snp_aa_positions(snps, cds, keep = "hits") |>
  subset(gene_id == "PF3D7_0810800")
```

On a real callset that last step resolves ~28,000 genotyped SNPs to
~17,000 coding ones across ~3,400 genes in about a second.

## Calling the residues a panel is missing

Filtering removes SNPs, and some of what it removes is a residue you
meant to keep. Going the other way – codons minus what the panel already
holds – gives a BED to re-call from the BAMs. See
[`vignette("filling-in-missed-calls")`](https://nickjhathaway.github.io/plasgenomicsutilsR/articles/filling-in-missed-calls.md).
