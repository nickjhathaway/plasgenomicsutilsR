# IBD post-analysis results

A container for the per-SNP, pairwise-region, and selection-statistic
tables produced by the Python `plasgenomicsutils ibd` tools, plus the
shared cumulative-genome coordinate the genome-wide plots use. Pass file
paths or data frames; each argument is optional so you can hold only
what you plan to plot. The `plot_*()` functions read from an object of
this class.

Expected columns (superset; extras are kept):

- `per_snp_region`: `chr`, `pos`, `frac_pairs_ibd`, optionally `region`
  (output of `analyze_ibd_matrix` per-SNP / per-SNP-per-region).

- `pairwise_region`: `chr`, `pos`, `region_a`, `region_b`,
  `frac_pairs_ibd` (per-SNP pairwise-region output).

- `selection`: `chr`, `pos`, a metric column (`neg_log10_p`,
  `chi2_stat`, or `z_score`), optionally `region` and `significant`
  (output of `ibd_selection_statistic`).

## Methods

### Public methods

- [`IbdResults$new()`](#method-IbdResults-initialize)

- [`IbdResults$get_per_snp_region()`](#method-IbdResults-get_per_snp_region)

- [`IbdResults$get_pairwise_region()`](#method-IbdResults-get_pairwise_region)

- [`IbdResults$get_selection()`](#method-IbdResults-get_selection)

- [`IbdResults$get_thresholds()`](#method-IbdResults-get_thresholds)

- [`IbdResults$get_genes()`](#method-IbdResults-get_genes)

- [`IbdResults$chrom_layout()`](#method-IbdResults-chrom_layout)

- [`IbdResults$reference_id()`](#method-IbdResults-reference_id)

- [`IbdResults$print()`](#method-IbdResults-print)

- [`IbdResults$clone()`](#method-IbdResults-clone)

------------------------------------------------------------------------

### `IbdResults$new()`

Create an IbdResults object.

#### Usage

    IbdResults$new(
      per_snp_region = NULL,
      pairwise_region = NULL,
      selection = NULL,
      threshold = NULL,
      genes = NULL,
      reference = "pf3d7"
    )

#### Arguments

- `per_snp_region`:

  Per-SNP (optionally per-region) IBD table: path or data frame.

- `pairwise_region`:

  Per-SNP pairwise-region IBD table: path or data frame.

- `selection`:

  IBD selection-statistic table: path or data frame.

- `threshold`:

  Significance threshold(s): a scalar, a named vector or
  `region`/`threshold` data frame (per region), or a path to a
  one-number file.

- `genes`:

  Optional gene-annotation track (`name`, `chr`, `start`, `end`) drawn
  as reference lines on the Manhattan plots.

- `reference`:

  Reference id for chromosome lengths (default `"pf3d7"`).

#### Returns

An `IbdResults` object (invisibly self).

------------------------------------------------------------------------

### `IbdResults$get_per_snp_region()`

Per-SNP (optionally per-region) IBD table with `cum_pos`.

#### Usage

    IbdResults$get_per_snp_region()

------------------------------------------------------------------------

### `IbdResults$get_pairwise_region()`

Per-SNP pairwise-region IBD table with `cum_pos`.

#### Usage

    IbdResults$get_pairwise_region()

------------------------------------------------------------------------

### `IbdResults$get_selection()`

Selection-statistic table with `cum_pos`.

#### Usage

    IbdResults$get_selection()

------------------------------------------------------------------------

### `IbdResults$get_thresholds()`

Threshold tibble (`region`, `threshold`) or `NULL`.

#### Usage

    IbdResults$get_thresholds()

------------------------------------------------------------------------

### `IbdResults$get_genes()`

Gene-annotation track with `cum_mid`, or `NULL`.

#### Usage

    IbdResults$get_genes()

------------------------------------------------------------------------

### `IbdResults$chrom_layout()`

Chromosome layout tibble (offsets, bands, axis mid-points).

#### Usage

    IbdResults$chrom_layout()

------------------------------------------------------------------------

### `IbdResults$reference_id()`

The reference id used for chromosome lengths.

#### Usage

    IbdResults$reference_id()

------------------------------------------------------------------------

### `IbdResults$print()`

Compact summary of what the object holds.

#### Usage

    IbdResults$print(...)

------------------------------------------------------------------------

### `IbdResults$clone()`

The objects of this class are cloneable with this method.

#### Usage

    IbdResults$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
