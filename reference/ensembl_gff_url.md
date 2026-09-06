# Ensembl Protists download URLs

The GFF3 annotation and the genome FASTA of a *Plasmodium* genome, so a
released annotation can be handed straight to
[`read_gff_cds()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/read_gff_cds.md)
or
[`snp_aa_positions()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/snp_aa_positions.md)
without looking the path up by hand. Nothing is downloaded – these
return the URL.

## Usage

``` r
ensembl_gff_url(
  species,
  release = ENSEMBL_PROTISTS_RELEASE,
  assembly = NULL,
  collection = NULL
)

ensembl_genome_url(
  species,
  release = ENSEMBL_PROTISTS_RELEASE,
  masking = c("none", "soft", "hard"),
  assembly = NULL,
  collection = NULL
)
```

## Arguments

- species:

  Species name; see details.

- release:

  Ensembl Protists release number. Defaults to
  [ENSEMBL_PROTISTS_RELEASE](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ENSEMBL_PROTISTS_RELEASE.md).

- assembly:

  Assembly token to use, skipping the species lookup. Optional.

- collection:

  FTP collection directory, or `""` for none. Optional; only needed
  alongside `assembly` for a genome outside the table.

- masking:

  Repeat masking of the genome FASTA: `"none"` (the default), `"soft"`
  (repeats in lower case) or `"hard"` (repeats as `N`).

## Value

One URL, as a string.

## Details

`species` may be the Ensembl name (`"plasmodium_falciparum"`), the
species alone (`"falciparum"`), `"P. falciparum"`, or a two-letter short
name (`"pf"`, `"pv"`, `"pk"`, `"pb"`, `"pc"`, `"py"`). The short forms
name the genome Ensembl treats as that species' reference – the one
Ensembl gives a top-level directory. Species with no such genome, *P.
malariae* and *P. ovale* among them, have several competing assemblies
and so must be named in full; the error lists them, and
[`ensembl_species()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ensembl_species.md)
lists them all.

The assembly token in every Ensembl filename follows whatever assembly
the release carried, and does change (Pf3D7 is `ASM276v2` through
release 60 and `GCA000002765v3` from release 61), so it is looked up for
the release asked for rather than assumed. Pass `assembly` to skip that
lookup, which also builds URLs for a genome outside the table.

## See also

[`ensembl_species()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ensembl_species.md),
[`read_gff_cds()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/read_gff_cds.md),
[`snp_aa_positions()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/snp_aa_positions.md)

## Examples

``` r
ensembl_gff_url("falciparum")
#> [1] "https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists/release-63/gff3/plasmodium_falciparum/Plasmodium_falciparum.GCA000002765v3.63.gff3.gz"
ensembl_genome_url("falciparum")
#> [1] "https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists/release-63/fasta/plasmodium_falciparum/dna/Plasmodium_falciparum.GCA000002765v3.dna.toplevel.fa.gz"

# any of these name the same genome
ensembl_gff_url("pf")
#> [1] "https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists/release-63/gff3/plasmodium_falciparum/Plasmodium_falciparum.GCA000002765v3.63.gff3.gz"
ensembl_gff_url("P. falciparum")
#> [1] "https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists/release-63/gff3/plasmodium_falciparum/Plasmodium_falciparum.GCA000002765v3.63.gff3.gz"

# other species, and the assembly of an earlier release
ensembl_gff_url("vivax")
#> [1] "https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists/release-63/gff3/plasmodium_vivax_gca900093555/Plasmodium_vivax_gca900093555.GCA900093555v2.63.gff3.gz"
ensembl_genome_url("plasmodium_malariae_gca_900090045", masking = "soft")
#> [1] "https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists/release-63/fasta/protists_alveolata1_collection/plasmodium_malariae_gca_900090045/dna/Plasmodium_malariae_gca_900090045.PmUG01.dna_sm.toplevel.fa.gz"

if (FALSE) { # \dontrun{
# a whole released annotation, read straight from the web
cds <- read_gff_cds(ensembl_gff_url("falciparum"))
snp_aa_positions(snps, cds, fasta = ensembl_genome_url("falciparum"))
} # }
```
