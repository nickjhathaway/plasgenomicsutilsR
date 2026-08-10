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
if (FALSE) { # \dontrun{
# VEuPathDB / PlasmoDB, the source the bundled gene datasets were built from
cds <- read_gff_cds(paste0("https://plasmodb.org/common/downloads/Current_Release/",
                           "Pfalciparum3D7/gff/data/PlasmoDB-68_Pfalciparum3D7.gff"))
# Ensembl Protists works too, despite naming its attributes differently
cds <- read_gff_cds(paste0("https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists/current/",
                           "gff3/plasmodium_falciparum/",
                           "Plasmodium_falciparum.GCA000002765v3.63.gff3.gz"))
} # }
```
