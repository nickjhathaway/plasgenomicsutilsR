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
  ...
)
```

## Arguments

- filename:

  Output path.

- plot:

  Plot to save (defaults to the last plot displayed).

- device:

  Graphics device. `NULL` (default) auto-selects for `.pdf` output (see
  [`pdf_device()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pdf_device.md));
  pass [`grDevices::cairo_pdf`](https://rdrr.io/r/grDevices/cairo.html),
  `"pdf"`, or any device
  [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)
  accepts to force a choice.

- width, height:

  Output size in inches. `NULL` (default) uses the size the `plot_*()`
  function attached to the plot (see
  [`plot_dims()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_dims.md)),
  falling back to `ggsave()`'s default if none is present.

- ...:

  Passed to
  [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)
  (e.g. `dpi`, `units`).

## Value

`filename`, invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
# size is chosen automatically from the number of regions / chromosomes drawn
save_plot("ibd_manhattan.pdf", plot_ibd_manhattan(example_ibd_results()))
} # }
```
