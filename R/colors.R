# Categorical colour palettes and per-metadata colour assignment, so the same
# level -> colour mapping can be shared across plots (UMAP points, admixture bars).
# Palettes mirror the colour-blind-friendly sets used in HaplotypeRainbows.

# Colour-blind-friendly base palettes (8, 12, 15 colours), mirroring the sets used in
# HaplotypeRainbows. Kept internal (dot-prefixed) so these generic names never collide with
# a user's globals; the exported `color_palette()` is the only entry point.
.color_palettes <- list(
  `8` = c("#2271B2", "#3DB7E9", "#F748A5", "#359B73",
          "#D55E00", "#E69F00", "#F0E442", "#000000"),
  `12` = c("#9F0162", "#009F81", "#FF5AAF", "#00FCCF", "#8400CD", "#008DF9",
           "#00C2F9", "#FFB2FD", "#A40122", "#E20134", "#FF6E3A", "#FFC33B"),
  `15` = c("#68023F", "#008169", "#EF0096", "#00DCB5", "#FFCFE2",
           "#003C86", "#9400E6", "#009FFA", "#FF71FD", "#7CFFFA",
           "#6A0213", "#008607", "#F60239", "#00E307", "#FFDC3D")
)

#' A colour-blind-friendly categorical palette
#'
#' Returns `n` colours from the package's colour-blind-friendly sets: the 8-, 12- or
#' 15-colour set that fits (interpolated beyond 15). Use it to colour metadata levels
#' consistently across plots. A single, namespaced entry point replaces the older
#' `colorPalette_08` / `colorPalette_12` / `colorPalette_15` objects.
#'
#' @param n Number of colours to return.
#' @return A character vector of `n` hex colours.
#' @examples
#' color_palette(5)
#' if (requireNamespace("scales", quietly = TRUE)) scales::show_col(color_palette(12))
#' @export
color_palette <- function(n) {
  if (n <= 0) return(character(0))
  if (n <= 8)  return(.color_palettes[["8"]][seq_len(n)])
  if (n <= 12) return(.color_palettes[["12"]][seq_len(n)])
  if (n <= 15) return(.color_palettes[["15"]][seq_len(n)])
  grDevices::colorRampPalette(.color_palettes[["15"]])(n)
}

# internal alias used throughout the plotting code
.pick_palette <- function(n) color_palette(n)

.levels_of <- function(x) {
  if (is.factor(x)) levels(droplevels(x)) else sort(unique(as.character(x[!is.na(x)])))
}

#' Assign colours to the levels of metadata columns
#'
#' Builds a named list (one entry per column) of `level -> hex colour` vectors, so a
#' single mapping can colour points in one plot and bars in another. Levels are taken
#' in factor order (or sorted); a colour-blind-friendly palette is chosen by level
#' count (interpolated past 15). `overrides` replaces specific colours.
#'
#' @param meta A data frame of metadata.
#' @param cols Columns to build colours for (default: all but a `sample` column).
#' @param overrides Optional named list `column -> (level -> colour)` to override
#'   individual assignments (e.g. `list(country = c(Ethiopia = "#FFDC3D"))`).
#' @return A named list, one `level -> colour` named vector per column.
#' @examples
#' meta_colors(data.frame(sample = 1:3, region = c("A", "B", "A")))
#' @export
meta_colors <- function(meta, cols = NULL, overrides = NULL) {
  if (is.null(cols)) cols <- setdiff(names(meta), "sample")
  auto <- stats::setNames(lapply(cols, function(cc) {
    lv <- .levels_of(meta[[cc]])
    stats::setNames(.pick_palette(length(lv)), lv)
  }), cols)
  if (!is.null(overrides)) {
    for (nm in names(overrides)) {
      if (nm %in% names(auto)) auto[[nm]][names(overrides[[nm]])] <- overrides[[nm]]
      else auto[[nm]] <- overrides[[nm]]
    }
  }
  auto
}
