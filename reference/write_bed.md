# Write an interval table as a BED file

Three columns, tab separated, no header, `start` 0-based half-open –
what `bedtools` and `bcftools mpileup -R` expect. A fourth `name` column
is written when the table has one, since a BED that says what each
interval is survives being looked at later.

## Usage

``` r
write_bed(
  x,
  file,
  name = NULL,
  sort = TRUE,
  chrom = NULL,
  name_is_coords = FALSE,
  make_names_unique = FALSE
)
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
  want. Numbering for uniqueness follows the written order.

- chrom:

  Column holding the chromosome name to write. Defaults to the first of
  a `<assembly>_chrom` column, `"chrom"`, and `"chr"` that the table has
  (see above).

- name_is_coords:

  Write `chrom-start-end` (`Pf3D7_07_v3-403221-403363`) as the name
  field instead of a column (default `FALSE`).

- make_names_unique:

  Make the name field unique within the file (default `FALSE`): a name
  that occurs more than once gets `_1`, `_2`, ... in file order,
  zero-padded to the width of that name's count (`_01` to `_12` for
  twelve). Combines with `name_is_coords`, so identical coordinates
  written twice get distinct names. With no name column and
  `name_is_coords = FALSE`, the coordinates are used as the names to
  make unique.

## Value

`file`, invisibly.

## Details

**The reference's own spelling is preferred.** A BED is read by other
tools against a real reference, so it has to carry the name the FASTA
and the BAMs use – writing a normalised one produces a file that matches
nothing, silently. Tables in this package carry two kinds of chromosome
column, and the default picks the one with the full name: a
`<assembly>_chrom` column first (`Pf3D7_chrom` in
[PF3D7_GENES](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF3D7_GENES.md)
and the other bundled datasets, where `chrom` is the short `"7"`), then
`chrom` (which
[`aa_intervals()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/aa_intervals.md),
[`tandem_repeats()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/tandem_repeats.md)
and friends fill with the source spelling, keeping `chr` for the
normalised `"7"`), then `chr`. Pieces from
[`bed_subtract()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/bed_subtract.md)
keep their `locs1` columns, so a gene table cut by a mask still writes
the right names.

## See also

[`bed_subtract()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/bed_subtract.md)

## Examples

``` r
iv <- data.frame(chr = c("Pf3D7_07_v3", "Pf3D7_07_v3"), start = c(403623, 403700),
                 end = c(403626, 403703), name = c("pfcrt-76", "pfcrt-102"))
write_bed(iv, file.path(tempdir(), "targets.bed"))
```
