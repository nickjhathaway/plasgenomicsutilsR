# Save a plot, preferring the cairo PDF device

A thin wrapper around
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)
that, for `.pdf` output, defaults to the cairo PDF device (better font
embedding) via
[`pdf_device()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pdf_device.md),
falling back to the standard `pdf` device where cairo is unavailable or
unreliable (e.g. Windows). Non-PDF outputs, and any explicit `device`,
pass straight through to `ggsave()`.

## Usage

``` r
save_plot(
  filename,
  plot = ggplot2::last_plot(),
  device = NULL,
  width = NULL,
  height = NULL,
  fit = TRUE,
  ...
)
```

## Arguments

- filename:

  Output path.

- plot:

  A single plot (defaults to the last plot displayed), **or a
  list/vector of plots** – e.g. the list returned by
  `plot_pairwise_ibd_for_genes(individual = TRUE)` – which is written as
  a **multi-page PDF**, one plot per page (`filename` must end in
  `.pdf`).

- device:

  Graphics device. `NULL` (default) auto-selects for `.pdf` output (see
  [`pdf_device()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pdf_device.md));
  pass [`grDevices::cairo_pdf`](https://rdrr.io/r/grDevices/cairo.html),
  `"pdf"`, or any device
  [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)
  accepts to force a choice.

- width, height:

  Output size in inches. **Give one and the other is computed** from the
  plot's contents; give neither and both are worked out; give both and
  they are used as-is. `NULL` for both falls back to the size the
  `plot_*()` function attached (see
  [`plot_dims()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_dims.md)),
  then to `ggsave()`'s default. For a multi-page list, one page size
  serves every page: the width every page needs, and the tallest height
  any page needs at that width.

- fit:

  Size the canvas to the drawing (default `TRUE`). A plot whose panel
  has a locked aspect ratio – anything using
  [`ggplot2::coord_fixed()`](https://ggplot2.tidyverse.org/reference/coord_fixed.html),
  such as
  [`plot_ibd_network()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_network.md)
  or the gene triangles – only fills a canvas of one particular shape;
  on any other shape the leftover appears as blank margin above and
  below (or beside) the drawing. Fitting measures the built plot's fixed
  furniture (titles, legends, axes, margins) and its panel ratio, then
  solves for the dimension you did not supply so there is no leftover.
  Set `FALSE` to use the requested / attached numbers verbatim.

- ...:

  Passed to
  [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)
  (e.g. `dpi`, `units`) for a single plot; ignored for a multi-page list
  (which honours only `width` / `height` / `device`).

## Value

`filename`, invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
# size is chosen automatically from the number of groups / chromosomes drawn
save_plot("ibd_manhattan.pdf", plot_ibd_sharing_manhattan(example_ibd_results()))
# a list of plots -> one multi-page PDF
save_plot("triangles.pdf",
          plot_pairwise_ibd_for_genes(example_ibd_results(), individual = TRUE))
} # }
```
