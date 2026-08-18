# Write an interval table as a BED file

Three columns, tab separated, no header, `start` 0-based half-open –
what `bedtools` and `bcftools mpileup -R` expect. A fourth `name` column
is written when the table has one, since a BED that says what each
interval is survives being looked at later.

## Usage

``` r
write_bed(x, file, name = NULL, sort = TRUE, chrom = NULL)
```

## Arguments

- x:

  An interval table (`chrom`/`chr`, `start`, `end`), e.g. from
  [`bed_subtract()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/bed_subtract.md).

- file:

  Path to write.

- name:

  Column to use as the BED name field, or `NULL` for none. Defaults to
  `"name"` when the table has it.

- sort:

  Sort by chromosome and start (default `TRUE`), which is what the tools
  want.

- chrom:

  Column holding the chromosome name to write. Defaults to `"chrom"`
  when the table has it, else `"chr"`.

## Value

`file`, invisibly.

## Details

**`chrom` is preferred over `chr`.** Tables in this package carry both:
`chr` normalised for matching (`"7"`), and `chrom` as the source file
spells it (`"Pf3D7_07_v3"`). A BED is read by other tools against a real
reference, so it has to carry the name the FASTA and the BAMs use –
writing the normalised one produces a file that matches nothing,
silently.

## See also

[`bed_subtract()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/bed_subtract.md)

## Examples

``` r
iv <- data.frame(chr = c("Pf3D7_07_v3", "Pf3D7_07_v3"), start = c(403623, 403700),
                 end = c(403626, 403703), name = c("pfcrt-76", "pfcrt-102"))
write_bed(iv, file.path(tempdir(), "targets.bed"))
```
