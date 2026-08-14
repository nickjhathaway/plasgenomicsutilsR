# Read the CDS features of a GFF

The coding exons only, with the transcript each belongs to – what
[`aa_intervals()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/aa_intervals.md)
needs to walk a protein back onto the genome. Reading the GFF is the
slow part, so parse once and reuse the result across calls.

## Usage

``` r
read_gff_cds(gff)
```

## Arguments

- gff:

  Path to a GFF3 file, or a URL – a plain or gzipped file is read
  straight from the web, so a released annotation can be used without
  keeping a copy.

## Value

A tibble of `transcript_id`, `gene_id`, `chrom`, `start`, `end` (1-based
inclusive, as the GFF gives them), `strand` and `phase`, one row per CDS
exon.

## See also

[`aa_intervals()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/aa_intervals.md)

## Examples

``` r
# the CDS of six 3D7 drug-resistance genes ship with the package
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

if (FALSE) { # \dontrun{
# a whole released annotation, read straight from the web
cds <- read_gff_cds(paste0("https://plasmodb.org/common/downloads/Current_Release/",
                           "Pfalciparum3D7/gff/data/PlasmoDB-68_Pfalciparum3D7.gff"))
# Ensembl Protists works too, despite naming its attributes differently
cds <- read_gff_cds(paste0("https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists/current/",
                           "gff3/plasmodium_falciparum/",
                           "Plasmodium_falciparum.GCA000002765v3.63.gff3.gz"))
} # }
```
