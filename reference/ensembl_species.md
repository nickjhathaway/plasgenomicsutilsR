# Plasmodium genomes in an Ensembl Protists release

The names
[`ensembl_gff_url()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ensembl_gff_url.md)
and
[`ensembl_genome_url()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ensembl_gff_url.md)
accept, with the assembly each carries in that release. The default
release's table ships with the package; any other release is read from
the Ensembl FTP site once per session.

## Usage

``` r
ensembl_species(release = ENSEMBL_PROTISTS_RELEASE, pattern = NULL)
```

## Arguments

- release:

  Ensembl Protists release number. Defaults to
  [ENSEMBL_PROTISTS_RELEASE](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ENSEMBL_PROTISTS_RELEASE.md).

- pattern:

  Optional regular expression, matched against the species name – e.g.
  `"vivax"` for the *P. vivax* genomes.

## Value

A tibble of `species` (the name to pass on), `assembly` (the token that
appears in the filenames) and `collection` (the FTP collection directory
the genome sits in, `""` for the genomes Ensembl gives a top-level
directory).

## See also

[`ensembl_gff_url()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ensembl_gff_url.md),
[`ensembl_genome_url()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ensembl_gff_url.md)

## Examples

``` r
ensembl_species()
#> # A tibble: 59 × 3
#>    species                                     assembly       collection        
#>    <chr>                                       <chr>          <chr>             
#>  1 plasmodium_berghei                          GCA900002375v2 ""                
#>  2 plasmodium_berghei_gca_900044335            PbK173         "protists_alveola…
#>  3 plasmodium_berghei_gca_900088445            PbNK65E        "protists_alveola…
#>  4 plasmodium_berghei_gca_900095585            PbSP11RLL      "protists_alveola…
#>  5 plasmodium_berghei_gca_900095635            PbSP11A        "protists_alveola…
#>  6 plasmodium_chabaudi                         PCHAS01        ""                
#>  7 plasmodium_chabaudi_adami_gca_900095565     PchDK          "protists_alveola…
#>  8 plasmodium_chabaudi_chabaudi_gca_900095605  PchCB          "protists_alveola…
#>  9 plasmodium_coatneyi_gca_001680005           ASM168000v1    "protists_alveola…
#> 10 plasmodium_cynomolgi_strain_b_gca_000321355 PcynB_1.0      "protists_alveola…
#> # ℹ 49 more rows

# the several assemblies of one species
ensembl_species(pattern = "ovale")
#> # A tibble: 5 × 3
#>   species                                  assembly collection                  
#>   <chr>                                    <chr>    <chr>                       
#> 1 plasmodium_ovale_curtisi_gca_900088555   Poc1     protists_alveolata1_collect…
#> 2 plasmodium_ovale_curtisi_gca_900088565   Poc2     protists_alveolata1_collect…
#> 3 plasmodium_ovale_gca_900090025           PowCR01  protists_alveolata1_collect…
#> 4 plasmodium_ovale_wallikeri_gca_900088485 Pow1     protists_alveolata1_collect…
#> 5 plasmodium_ovale_wallikeri_gca_900088545 Pow2     protists_alveolata1_collect…
```
