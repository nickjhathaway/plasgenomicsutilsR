# The preferred PDF graphics device

Returns
[`grDevices::cairo_pdf()`](https://rdrr.io/r/grDevices/cairo.html) when
a cairo device can actually be opened and the platform is not Windows
(cairo embeds fonts more reliably, but its PDF output can be unreliable
on Windows), otherwise the string `"pdf"`. Availability is settled by
opening a throwaway device rather than by `capabilities("cairo")`, which
reports what R was built with and so can claim a cairo that fails to
load. Use it with
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html):
`ggsave(file, plot, device = pdf_device())`.

## Usage

``` r
pdf_device()
```

## Value

A device function
([grDevices::cairo_pdf](https://rdrr.io/r/grDevices/cairo.html)) or the
string `"pdf"`.

## Examples

``` r
pdf_device()
#> function (filename = if (onefile) "Rplots.pdf" else "Rplot%03d.pdf", 
#>     width = 7, height = 7, pointsize = 12, onefile = TRUE, family = "sans", 
#>     bg = "white", antialias = c("default", "none", "gray", "subpixel"), 
#>     fallback_resolution = 300, symbolfamily) 
#> {
#>     if (!checkIntFormat(filename)) 
#>         stop("invalid 'filename'")
#>     if (!capabilities("cairo")) 
#>         stop("cairo_pdf: Cairo-based devices are not available for this platform")
#>     antialiases <- eval(formals()$antialias)
#>     antialias <- match(match.arg(antialias, antialiases), antialiases)
#>     if (missing(symbolfamily)) 
#>         symbolfamily <- symbolfamilyDefault(family)
#>     invisible(.External(C_devCairo, filename, 6L, 72 * width, 
#>         72 * height, pointsize, bg, NA_integer_, antialias, onefile, 
#>         family, fallback_resolution, checkSymbolFont(symbolfamily)))
#> }
#> <bytecode: 0x564092f5f8a0>
#> <environment: namespace:grDevices>
```
