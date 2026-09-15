# Integrated haplotype score (iHS)

Scans each group for recent positive directional selection,
standardising the integrated EHH ratio within allele-frequency bins so
scores are comparable along the genome. Large `abs(ihs)` marks a SNP
whose haplotype background is unusually long for its frequency.

## Usage

``` r
run_ihs(
  hap,
  group = NULL,
  meta = NULL,
  polarized = FALSE,
  freqbin = NULL,
  min_maf = 0.05,
  maf_bands = NULL,
  min_samples = 4,
  maxgap = NA,
  scalegap = NA,
  discard_at_border = NULL,
  threads = 1,
  contrast = c("ref", "pairwise", "none")
)
```

## Arguments

- hap:

  A
  [`parasite_haplotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/parasite_haplotypes.md)
  object.

- group:

  Metadata column naming the grouping, a vector aligned to the haplotype
  rows, or `NULL` to scan every sample as one population.

- meta:

  Metadata (defaults to the one carried by `hap`).

- polarized:

  Treat allele 1 as derived (needs a real ancestral state).

- freqbin:

  Width of the allele-frequency bins iHS is standardised within, or
  `NULL` (default) to pick one from `polarized`: **1** (a single bin)
  when unpolarized, `0.05` when polarized. The binning exists to control
  for *derived* allele frequency, and an unpolarized scan has no
  ancestral state – only `FREQ_MAJ`/`FREQ_MIN` – so major/minor is not
  derived/ancestral and binning by it controls nothing. rehh warns about
  this and about the resulting sparse bins above 0.5; a single bin
  silences both because it is the right answer, not because the warning
  was noise. Not cosmetic: on a real 249-sample cohort the two settings
  share only 25 of their top 50 \|iHS\| hits, and 4-12% of SNPs change
  sign (most in the 0.05-0.1 and \>0.3 MAF bands).

- min_maf:

  Minor-allele frequency floor applied at standardisation.

- maf_bands:

  Standardise within this many minor-allele-frequency bands, cut at
  quantiles of the observed frequencies so each band holds a similar
  number of markers. `NULL` (default) leaves the standardisation to rehh
  and `freqbin`. The spread of the log iHH ratio depends on how common
  the minor allele is – a rarer allele is carried by fewer haplotypes,
  so its integrals are noisier – and a single bin cannot remove that,
  which leaves rarer SNPs over-represented in the tail. Bands do remove
  it, without running into the empty bins that make rehh's own binning
  unusable on an unpolarized scan. Around 10 is reasonable; it replaces
  `freqbin` rather than combining with it, and it changes every score,
  so a result computed with it is not comparable to one computed
  without.

- min_samples:

  Skip groups smaller than this.

- maxgap:

  Largest gap between consecutive SNPs, in base pairs, that the EHH
  integration may cross; `NA` (the default, and rehh's) lets it cross
  any gap. This matters more than its default suggests. A region with no
  SNPs – a centromere, a masked hypervariable block – has nothing to
  break the haplotype, so EHH runs flat across it and the integral
  accumulates `EHH x gap length`. The SNPs flanking such a hole then
  score on the width of the hole rather than on their haplotypes, in
  either direction: the ratio is diluted towards zero when both alleles
  carry EHH into the gap, and inflated when only one does. Pick a value
  from the data's own spacing (several dozen times the median gap leaves
  ordinary density untouched while stopping at a real hole) rather than
  from a round number.

- scalegap:

  Gaps wider than this are counted as being exactly this wide, rather
  than stopping the integration outright; `NA` (default) does not
  rescale. A softer form of `maxgap` – it caps a hole's contribution
  instead of refusing to cross it.

- discard_at_border:

  Return `NA` instead of a truncated integral when the integration runs
  into the end of a chromosome or a gap wider than `maxgap`. `NULL`
  (default) ties it to `maxgap`: off when no `maxgap` is set (so the
  markers nearest the telomeres are still scored), on when one is. On
  sparse markers this can empty the scan – if EHH never decays before
  the data runs out, every marker is at a border – so a scan that comes
  back mostly `NA` says so.

- threads:

  Threads for rehh.

- contrast:

  For a multiallelic marker, which allele pairs to score. `"ref"`
  (default) contrasts each alternate against the reference – one `"0>k"`
  row per alternate, kept in the returned `contrast` column;
  `"pairwise"` scores every pair of alleles; `"none"` takes the single
  major-versus-minor value rehh reports and drops the multiallelic
  distinction. A biallelic marker has one contrast, so all three agree
  there.

## Value

A tibble with `group`, `chr`, `pos`, `snp_id`, `freq_minor`, `unihs`,
`ihs` and `neg_log10_p`.

`unihs` is the **un**standardised statistic,
`log(iHH_major / iHH_minor)` (ancestral over derived when
`polarized = TRUE`), so `exp(unihs)` is the integrated-EHH ratio itself.
`ihs` is that value z-scored within its frequency band, which is what
makes scores comparable along the genome but also throws the scale away
– the band's mean and sd are not recoverable from `ihs` alone, so keep
`unihs` if you ever want the ratio back. Note that `ihs = 0` does
**not** mean a ratio of 1: it means average for that group, and the
average is below 1 wherever minor-allele haplotypes are systematically
longer.

`ihs` and `neg_log10_p` are `NA` wherever the integral could not be
formed for one of the two alleles – see the note on missing scores
below.

## Details

Without an outgroup there is no ancestral state to polarise by, so
`polarized = FALSE` by default and the comparison is major versus minor
allele rather than ancestral versus derived. That is the standard
treatment for *P. falciparum* and it means the *sign* of `ihs` should
not be read as "selection on the derived allele" – use `abs(ihs)` and
`neg_log10_p`.

## Contrasts

`log(iHH_A / iHH_B)` is a ratio, so it needs exactly two terms and there
is no k-allele iHS. At a marker carrying three alleles the honest answer
is **k-1 contrasts**, each saying which two it compared, and the
`contrast` column names them as `"0>1"`, `"0>2"` and so on in
allele-index order (`$sites$alt` on the panel says which base each index
is).

`"ref"` compares each alternate against the reference only. That is what
keeps independent origins apart: at a codon where three changes arose
separately, pooling them into one "not reference" class merges exactly
the distinction the scan exists to draw.

What makes this cheap is that **per-allele iHH does not depend on which
other alleles are at the marker** – EHH for an allele class only ever
involves that class's haplotypes. So every contrast at a marker comes
out of one
[`rehh::calc_ehh()`](https://rdrr.io/pkg/rehh/man/calc_ehh.html) call,
and the numbers are identical to what an explicitly subset and recoded
panel would give (verified to nine decimal places). Only multiallelic
markers cost anything extra.

A biallelic marker has one contrast, so `"ref"` reproduces `"none"`
exactly, value for value. The standardisation is over every row at once,
so a marker contributing two rows has both judged against the same
genome-wide distribution.

Downstream, `contrast` is a grouping key:
[`ihs_windows()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ihs_windows.md),
[`ihs_genes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ihs_genes.md),
[`selection_peaks()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/selection_peaks.md)
and
[`plot_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ihs.md)
all split on it, because two contrasts at one position are two
measurements rather than two SNPs at one site.

## Why a SNP can appear for one group only, or score `NA`

The scan is per group, so a SNP is tested in a group only where it is
polymorphic there and clears `min_maf` there. A variant private to one
region therefore has one row, not one per region, and that is a
statement about the cohort rather than a fault.

A row can be present with `ihs` and `neg_log10_p` both `NA`. That means
the SNP passed the frequency filter but rehh could not integrate EHH for
at least one of its two alleles, so the log ratio is undefined. The
usual causes are too few haplotypes carrying the minor allele for the
decay to be estimated, and EHH that never falls below the cutoff before
the data runs out – a chromosome end, or a gap wider than `maxgap`. Both
get more common in small groups and at low minor-allele counts, which is
the same corner where a score that *is* returned deserves the least
trust. Treat `NA` as "not measurable here", not as "no selection".

## References

Voight, B. F., Kudaravalli, S., Wen, X. & Pritchard, J. K. (2006) A map
of recent positive selection in the human genome. *PLoS Biology* 4, e72.
[doi:10.1371/journal.pbio.0040072](https://doi.org/10.1371/journal.pbio.0040072)

Gautier, M., Klassmann, A. & Vitalis, R. (2017) rehh 2.0: a
reimplementation of the R package rehh to detect positive selection from
haplotype structure. *Molecular Ecology Resources* 17, 78-90.
[doi:10.1111/1755-0998.12634](https://doi.org/10.1111/1755-0998.12634)

## See also

[`ihs_windows()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ihs_windows.md),
[`ihs_genes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ihs_genes.md),
[`plot_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ihs.md),
[`run_rsb()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_rsb.md),
[`beta_score()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/beta_score.md)

## Examples

``` r
ps <- example_pop_structure(umap = FALSE)
hap <- parasite_haplotypes(ps, maf = 0.05)
run_ihs(hap, group = "country")
#> # A tibble: 365 × 8
#>    group    chr             pos snp_id    freq_minor   unihs     ihs neg_log10_p
#>    <fct>    <chr>         <dbl> <chr>          <dbl>   <dbl>   <dbl>       <dbl>
#>  1 Cambodia Pf3D7_02_v3  273786 Pf3D7_02…     0.467  -0.0681 -0.215       0.0811
#>  2 Cambodia Pf3D7_04_v3   92596 Pf3D7_04…     0.1    -0.492  -0.499       0.209 
#>  3 Cambodia Pf3D7_04_v3  544672 Pf3D7_04…     0.467  -0.0253 -0.186       0.0695
#>  4 Cambodia Pf3D7_04_v3  898667 Pf3D7_04…     0.0667 -0.542  -0.533       0.226 
#>  5 Cambodia Pf3D7_04_v3 1006054 Pf3D7_04…     0.133  -0.0181 -0.182       0.0676
#>  6 Cambodia Pf3D7_07_v3   80465 Pf3D7_07…     0.367   0.203  -0.0332      0.0117
#>  7 Cambodia Pf3D7_07_v3  375090 Pf3D7_07…     0.0667 -0.216  -0.314       0.123 
#>  8 Cambodia Pf3D7_07_v3  375470 Pf3D7_07…     0.0667  4.87    3.09        2.71  
#>  9 Cambodia Pf3D7_07_v3  375817 Pf3D7_07…     0.0667  4.87    3.09        2.71  
#> 10 Cambodia Pf3D7_07_v3  376527 Pf3D7_07…     0.333   0.740   0.327       0.128 
#> # ℹ 355 more rows
```
