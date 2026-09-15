# EHH decay around one SNP

Extended haplotype homozygosity either side of a focal SNP, one curve
per allele: how far the haplotype carrying each allele stays identical
as you walk away from it. A sweep shows as one allele holding EHH near 1
far past the point where the other has decayed – the picture behind a
single point on an
[`run_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_ihs.md)
scan.

## Usage

``` r
plot_ehh(
  x,
  focal,
  group = NULL,
  span = 50000,
  min_haplotypes = 10,
  polarized = FALSE,
  limehh = 0.05,
  genes = NULL,
  gene_track = NULL,
  gene_label_angle = 0,
  colours = NULL,
  show_freq = TRUE,
  freq_position = c("topleft", "topright", "bottomleft", "bottomright"),
  reference = DEFAULT_REFERENCE,
  title = NULL,
  subtitle = NULL,
  add_ihs = NULL,
  ihs_args = list(),
  colors = NULL
)
```

## Arguments

- x:

  A
  [`parasite_haplotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/parasite_haplotypes.md)
  object, or a
  [PopStructure](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PopStructure.md)
  (haplotypes are then built with
  [`parasite_haplotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/parasite_haplotypes.md)
  defaults, which is worth doing yourself when the Fws or MAF cutoffs
  matter).

- focal:

  The SNP to measure from: a `chr:pos` id, a bare position, or a gene
  name from `genes`. A gene holding several SNPs resolves to the one
  with the most balanced alleles, reported in a message – name a
  `chr:pos` to pick a particular mutation.

- group:

  Optional metadata column; one panel per level. `NULL` (default) pools
  every haplotype.

- span:

  How far either side of the focal SNP to draw, in base pairs (default
  50 kb). One value is symmetric, two are the left and the right (named
  `left` / `right` if you like), as elsewhere in the package.

- min_haplotypes:

  Skip a group with fewer haplotypes than this (default 10). EHH from a
  handful of haplotypes is mostly noise.

- polarized:

  Treat the alleles as ancestral / derived (default `FALSE`, matching
  [`run_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_ihs.md)
  on unpolarized calls, where they are simply the two states).

- limehh:

  Stop each curve once EHH falls below this (rehh's `limehh`, default
  0.05).

- genes:

  Gene table for the track and for resolving `focal` (e.g.
  [PF3D7_GENES](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF3D7_GENES.md)).

- gene_track:

  Draw the gene track underneath (default `FALSE`). `genes` is usually
  supplied only to resolve `focal`, and an EHH window is wide enough
  that a full annotation would crowd a hundred names under it, so this
  is opt-in.

- gene_label_angle:

  Rotation for the gene names, in degrees.

- colours, colors:

  Named colours overriding the focal alleles' defaults. A biallelic
  marker's levels are `reference` and `alternate`; a multiallelic one's
  are `reference`, `alternate 1`, `alternate 2`, ..., numbered by
  descending frequency so `alternate 1` is the commonest alternate. The
  defaults are one shared colour-blind-safe palette across every EHH
  plot – `reference` always the same blue, `alternate` and `alternate 1`
  the same vermillion – so separate biallelic and multiallelic panels
  read together and the reference curve is the same colour in each. Name
  any subset to override, e.g. `colours = c("alternate 2" = "grey50")`.

- show_freq:

  Note each panel's haplotype count and allele frequencies inside it
  (default `TRUE`); `FALSE` leaves the panel clean.

- freq_position:

  Which corner that note sits in: `"topleft"` (default), `"topright"`,
  `"bottomleft"` or `"bottomright"`. The top corners are usually clear,
  since EHH is 1 at the focal SNP and both curves have flattened along
  the bottom by the window's edges.

- reference:

  Reference id, used when `focal` names a whole chromosome.

- title:

  Plot title: `NULL` (default) uses `"EHH around <snp>"`, a string sets
  a custom one, and `NA`/`FALSE` draws none. Set it here rather than
  adding `labs(title = )` to the result – with `gene_track = TRUE` the
  result is a patchwork, and `+ labs()` lands on the gene track at the
  bottom, putting the title in the middle of the figure.

- subtitle:

  Line under the title; `NULL` (default) draws none.

- add_ihs:

  Add the focal SNP's iHS to that corner note. A
  [`run_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_ihs.md)
  result is read for the focal SNP – the cheap path, and the one to
  prefer, since it reuses a scan you already have and so the number in
  the corner is the same one the genome-wide figures were drawn from.
  `TRUE` runs
  [`run_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_ihs.md)
  here instead, on `x`, with this plot's `group` and `polarized` and
  anything in `ihs_args`; that is a whole-genome scan per call, so it is
  slow and worth doing once into a variable rather than once per plot.
  `NULL` (default) or `FALSE` adds nothing. iHS is standardised against
  the whole genome and cannot be recovered from the window drawn here,
  which is why there is no third option. Read against the same caution
  [`run_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_ihs.md)
  carries: unpolarized, only the magnitude is shown, because the sign is
  major-versus-minor and not ancestral-versus-derived.

  A `contrast` column in the table asks the same focal marker more than
  one question – the reference against each alternate of a multiallelic
  codon, say. Each named set gets its own line in the corner, prefixed
  by its name, and a panel with no value for one of them simply omits
  that line rather than printing a blank. Without the column the table
  is read as a single unnamed contrast, as before.

- ihs_args:

  Extra arguments for the
  [`run_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_ihs.md)
  call made by `add_ihs = TRUE`, as a named list – `maxgap`,
  `maf_bands`, `min_maf` and the rest. `group` and `polarized` come from
  this plot so the two halves cannot disagree, and naming either here is
  an error. Ignored, with a warning, when `add_ihs` is a scan you
  computed yourself.

## Value

A patchwork of the curves over the gene track, or a plain ggplot without
one.

## Details

The alleles are the two states at the focal SNP itself, so this is the
mutant-versus- reference comparison without needing the SNPs annotated:
`reference` is the allele coded 0 and `alternate` the one coded 1.
`group` adds a panel per metadata group; without it every haplotype is
pooled, which is usually what you want first, since EHH knows nothing
about population structure and a group with few carriers gives a ragged
curve.

## Examples

``` r
ps <- example_pop_structure(umap = FALSE)
hap <- parasite_haplotypes(ps, maf = 0.05)
plot_ehh(hap, "pfcrt", genes = PF_EXAMPLE_DRUG_GENES, span = 30000)
#> `pfcrt` holds 11 SNPs; measuring from Pf3D7_07_v3:403624 (minor allele 0.5) -- name a `chr:pos` to pick another, or see ehh_candidates() for the shortlist
```
