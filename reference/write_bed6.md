# Write an interval table as a six-column BED, with optional metadata

The BED6 layout: `chrom`, `start`, `end`, `name`, `score`, `strand`, tab
separated, no header, `start` 0-based half-open. An optional seventh
column carries any other columns of the table as
`[field=value;field=value;]`, the form `elucidator` reads metadata in
(`meta_cols`; the name `meta` is kept for sample-metadata tables
throughout the package).

## Usage

``` r
write_bed6(
  x,
  file,
  name = NULL,
  score = NULL,
  strand = NULL,
  meta_cols = NULL,
  semicolon = ":",
  chrom = NULL,
  sort = TRUE,
  name_is_coords = FALSE,
  make_names_unique = FALSE
)
```

## Arguments

- x, file, chrom, sort:

  As for
  [`write_bed()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/write_bed.md).

- name, score, strand:

  Columns to write in those fields, or `NULL` for the defaults above.

- meta_cols:

  Columns to carry in the seventh column, as a character vector of
  names, or `TRUE` for every column not already written. `NULL` (the
  default) writes six columns.

- semicolon:

  What replaces a `;` inside a metadata field or value.

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

Defaults follow the table: `name` is its `name` column, else
`chrom-start-end`; `score` is its `score` column, else the interval's
width in bp; `strand` is its `strand` column, else `+`. The chromosome
column is chosen as
[`write_bed()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/write_bed.md)
does, preferring the reference's own spelling. Missing values in `name`,
`score` and `strand` are written as `.`.

Metadata fields and values are written as-is, whitespace included,
except that `;` – the separator – is replaced by `semicolon` (default
`:`) so a value can never split a field. Numbers are written in full,
never in scientific notation.

## See also

[`write_bed()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/write_bed.md)
for the three- or four-column form.

## Examples

``` r
targets <- bed_subtract(PF3D7_GENES[PF3D7_GENES$name == "pfcrt", ],
                        tandem_repeats_to_avoid(pf3d7_tandem_repeats()), pad = 10)
f <- file.path(tempdir(), "pfcrt_pieces.bed")
write_bed6(targets, f, meta_cols = c("gene_id", "piece"))
readLines(f)[1:2]
#> [1] "Pf3D7_07_v3\t403221\t403363\tpfcrt\t142\t+\t[gene_id=PF3D7_0709000;piece=1;]"
#> [2] "Pf3D7_07_v3\t403480\t403823\tpfcrt\t343\t+\t[gene_id=PF3D7_0709000;piece=2;]"

# every piece is called "pfcrt": number them, or name them by their coordinates
write_bed6(targets, f, make_names_unique = TRUE)
readLines(f)[1:2]
#> [1] "Pf3D7_07_v3\t403221\t403363\tpfcrt_01\t142\t+"
#> [2] "Pf3D7_07_v3\t403480\t403823\tpfcrt_02\t343\t+"
write_bed6(targets, f, name_is_coords = TRUE)
readLines(f)[1:2]
#> [1] "Pf3D7_07_v3\t403221\t403363\tPf3D7_07_v3-403221-403363\t142\t+"
#> [2] "Pf3D7_07_v3\t403480\t403823\tPf3D7_07_v3-403480-403823\t343\t+"
```
