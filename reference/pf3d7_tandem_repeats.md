# Tandem repeats of the Pf3D7 reference genome

Every short tandem repeat found in the *P. falciparum* 3D7 reference, as
a
[`tandem_repeats()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/tandem_repeats.md)
table, so nobody has to run the finders again. Feed it to
[`tandem_repeats_to_avoid()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/tandem_repeats_to_avoid.md)
and subtract the result from a set of targets with
[`bed_subtract()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/bed_subtract.md).

## Usage

``` r
pf3d7_tandem_repeats()
```

## Value

The
[`tandem_repeats()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/tandem_repeats.md)
table for Pf3D7: about 272,000 rows.

## Details

Built once (April 2021) from `Pf3D7.fasta` (PlasmoDB, the assembly
[PF3D7_GENES](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF3D7_GENES.md)
is annotated on) by the union of two finders:

- Tandem Repeats Finder (Benson 1999),
  `trf Pf3D7.fasta 2 7 7 80 10 50 1000 -f -d -m -h` – match 2, mismatch
  7, indel 7, match / indel probabilities 80 / 10, minimum score 50,
  maximum period 1000 – converted to BED with
  `elucidator TandemRepeatFinderOutputToBed`;

- `elucidator findSimpleTandemRepeatLocations --maxRepeatUnitSize 10`,
  an exhaustive search for perfect repeats of units up to 10 bp, which
  `trf` is inconsistent about around short simple repeats.

`copies` is the finder's own copy number, which for `trf` reflects an
alignment with indels and so is not always `width / unit_size`. All 16
sequences are covered, the apicoplast and mitochondrion included.

The table is stored compactly under `inst/extdata` (about 1.1 MB) and
expanded on first use; later calls return the cached copy.
`data-raw/PF3D7_TANDEM_REPEATS.R` rebuilds it.

## References

Benson G (1999). Tandem repeats finder: a program to analyze DNA
sequences. *Nucleic Acids Research* 27(2):573-580.

## See also

[`tandem_repeats_to_avoid()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/tandem_repeats_to_avoid.md),
[`bed_subtract()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/bed_subtract.md)

## Examples

``` r
tr <- pf3d7_tandem_repeats()
table(pmin(tr$period, 5), useNA = "ifany")     # 1..4, and "5" for anything longer
#> 
#>      1      2      3      4      5 
#>  37476  56813  16337  15388 145922 
```
