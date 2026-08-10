# Diversity, linkage and selection

Three questions about a cohort, each answered from the same genotype
matrix and the same metadata grouping: how diverse is each population,
how far does linkage reach, and which loci are under selection.

The example fixture is deliberately tiny (60 samples, 49 SNPs), so the
numbers below are shaped by that, not by biology. Everything scales to a
real callset unchanged.

``` r

ps <- example_pop_structure(umap = FALSE)
```

## Diversity

[`pop_diversity()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diversity.md)
reports, per group: nucleotide diversity, expected heterozygosity,
Watterson’s theta, Tajima’s D, and the haplotype / multilocus-genotype
summaries.

``` r

pop_diversity(ps, group = "country", accessible = PF3D7_CORE_REGIONS)
#> # A tibble: 2 × 24
#>   group    chr   start   end unit   n_samples n_snps  n_sites seg_sites     he
#>   <fct>    <chr> <dbl> <dbl> <chr>      <int>  <int>    <dbl>     <int>  <dbl>
#> 1 Cambodia NA       NA    NA genome        30    719 20857252       158 0.0538
#> 2 Ghana    NA       NA    NA genome        30    719 20857252       546 0.205 
#> # ℹ 14 more variables: pi <dbl>, theta_w <dbl>, tajima_d <dbl>, tajima_p <dbl>,
#> #   n_taj_snps <int>, n_taj_samples <dbl>, n_hap_snps <int>,
#> #   n_hap_samples <int>, n_hap <int>, hap_div <dbl>, shannon_h <dbl>,
#> #   simpson_lambda <dbl>, evenness <dbl>, tajima_percentile <dbl>
```

### pi and He are not the same number

This is the one thing worth being careful about. **`pi` is per
accessible base pair** – the sum of per-site heterozygosity divided by
the number of callable sites. **`he`** is that same per-site quantity
averaged over the **SNPs**. They differ by a factor of several thousand
and only one of them is comparable between windows:

``` r

d <- pop_diversity(ps, group = "country", accessible = PF3D7_CORE_REGIONS)
d[c("group", "n_snps", "n_sites", "he", "pi")]
#> # A tibble: 2 × 5
#>   group    n_snps  n_sites     he         pi
#>   <fct>     <int>    <dbl>  <dbl>      <dbl>
#> 1 Cambodia    719 20857252 0.0538 0.00000185
#> 2 Ghana       719 20857252 0.205  0.00000707
```

Dividing by the SNP count instead of the base count is the usual
mistake, and it makes windows with different SNP density or missingness
silently incomparable. `accessible` sets the denominator; without it
every base of the unit is assumed callable, which is only right if your
VCF really was called across the whole span.

### Per gene, or in windows

Per gene is the convention for *P. falciparum* Tajima’s D – coding SNPs,
at least three of them, no MAF filter:

``` r

pop_diversity(ps, group = "country", by = "gene",
              genes = PF_EXAMPLE_DRUG_GENES)[1:4, c("group", "unit", "n_snps", "he")]
#> # A tibble: 4 × 4
#>   group    unit   n_snps      he
#>   <fct>    <chr>   <int>   <dbl>
#> 1 Cambodia pfcrt      27  0.0701
#> 2 Cambodia pfdhfr      0 NA     
#> 3 Cambodia pfmdr1      0 NA     
#> 4 Cambodia pfdhps      7  0.231
```

Windows give a genome-wide track. `step` slides the window rather than
tiling it, which is the usual way to scan Tajima’s D — a fixed grid can
split a signal across two windows and dilute it in both:

``` r

win <- pop_diversity(ps, group = "country", by = "window", window = 500000, step = 250000)
plot_diversity(win, metric = "tajima_d")
```

![](diversity-and-selection_files/figure-html/unnamed-chunk-5-1.png)

On a real callset `window = 5000, step = 2500` is a reasonable scanning
resolution. Overlapping windows are not independent, so count *peaks*,
not windows — \[selection_peaks()\] does that.

`tajima_p` tests each unit against the standard neutral model, and
`tajima_percentile` places its D within the group’s own distribution.
Prefer the percentile for shortlisting: see
[`?tajima_d_pvalue`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/tajima_d_pvalue.md)
for why the parametric test is conservative here.

### Haploid genotypes, mixed infections

The parasite is haploid, so a heterozygous call is a mixed infection at
that site rather than a diploid genotype. By default those calls are
dropped at that site (`het = "missing"`) and every allele count is a
count of samples; `het = "dosage"` splits the call between the two
alleles instead.

Tajima’s D needs one sample count for its variance term, but sites
differ in how many calls they have. Rather than discard every sample
with a gap, all sites are kept and the mean number of calls is used as
*n* – `n_taj_snps` and `n_taj_samples` report what went in. The
haplotype-based columns (`hap_div`, `shannon_h`, `simpson_lambda`) do
need complete data, and there something has to give: either the samples
with a gap or the sites. Sites are the plentiful axis, so the sites with
any missing call are dropped and every sample stays in the comparison
(`n_hap_snps` reports what survived). Dropping samples instead would be
hopeless at scale — at 0.05% missingness over 27,000 SNPs a sample
survives with probability ~1e-6.

Note that haplotype diversity is a *per-locus* statistic: genome-wide it
saturates at 1, because over thousands of SNPs every sample is its own
haplotype. Use `by = "gene"` or `by = "window"` for it to say anything.

## Linkage disequilibrium

Two views, and they live on different sides of the split. **r-squared
decay** with physical distance is the one genuinely quadratic statistic
here – every SNP against every other SNP within `max_dist` – so it runs
in the Python package, where it takes about 7 seconds on a 249-sample,
28k-SNP callset:

``` bash
plasgenomicsutils ld_decay --vcf cohort.bcf --meta meta.tsv --group-col region \
  --max-dist 50000 --bins 20 --output ld/decay.tsv.gz --half-decay-output ld/half.tsv
```

[`read_ld_decay()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/read_ld_decay.md)
brings it back, with the half-decay distance and the scan settings
attached so the thinning behind a curve travels with it:

``` r

ld <- read_ld_decay("ld/decay.tsv.gz", "ld/half.tsv")
attr(ld, "ld_half_decay")
plot_ld_decay(ld)
```

**[`ld_index()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ld_index.md)**
is the multilocus view and stays here, being linear and quick: `Ia` and
its standardised form `rbarD`, both 0 under free recombination and
rising as reproduction becomes clonal. Compare `rbarD` between datasets,
since `Ia` grows with the number of loci.

``` r

ld_index(ps, group = "country")
#> # A tibble: 2 × 5
#>   group    n_samples n_loci    ia  rbar_d
#>   <fct>        <int>  <int> <dbl>   <dbl>
#> 1 Cambodia        30     92 11.8  0.141  
#> 2 Ghana           30    425  2.74 0.00672
```

Both are **upward-biased in small samples**, so read them next to
`n_samples` – a group of ten will score above a group of a hundred drawn
from the same population.

## Selection

Two complementary scans. Directional sweeps (drug resistance) show up in
extended haplotype homozygosity; long-term balancing selection
(antigens) shows up in clustered allele frequencies.

### Haplotypes first

needs phased haplotypes with no missing calls, and a *P. falciparum* VCF
has neither.
[`parasite_haplotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/parasite_haplotypes.md)
is the bridge, and it reports every sample and SNP it removed – how the
haplotypes were made determines what the scan can honestly claim.

``` r

hap <- parasite_haplotypes(ps, maf = 0.05)
hap
#> <parasite_haplotypes> 60 haplotypes x 351 SNPs
#>   from            : 60 samples x 719 SNPs
#>   mixed calls     : 814 resolved by allele draw 
#>   SNPs dropped    : 221 missing, 147 MAF
#>   samples dropped : 0 missing
#>   imputed calls   : 178  seed: 42
```

With per-sample Fws to hand (from `plasgenomicsutils calculate_fws`),
gate to monoclonal infections first, since a polyclonal infection is a
mixture of haplotypes rather than one:

``` r

fws <- readr::read_tsv("fws.tsv")
hap <- parasite_haplotypes(ps, fws = fws, min_fws = 0.95, maf = 0.03)
```

Both the mixed-call resolution and the imputation are random draws at
the population allele frequency. `seed` makes a run reproducible – and a
signal worth reporting should survive repeating the whole thing under a
different seed.

### iHS, Rsb and XP-EHH

``` r

ihs <- run_ihs(hap, group = "country")
head(ihs, 3)
#> # A tibble: 3 × 7
#>   group    chr            pos snp_id             freq_minor    ihs neg_log10_p
#>   <fct>    <chr>        <dbl> <chr>                   <dbl>  <dbl>       <dbl>
#> 1 Cambodia Pf3D7_02_v3 273786 Pf3D7_02_v3:273786      0.467 -0.215      0.0811
#> 2 Cambodia Pf3D7_04_v3  92596 Pf3D7_04_v3:92596       0.1   -0.499      0.209 
#> 3 Cambodia Pf3D7_04_v3 544672 Pf3D7_04_v3:544672      0.467 -0.186      0.0695
```

**How to read it.** At each SNP, iHS contrasts how far haplotype
homozygosity extends from one allele against the other. An allele that
rose recently and fast still carries the long stretch it arose on,
because recombination has not had time to break it up.

Four things follow from how it is built:

- **It is a z-score.** Raw haplotype length depends heavily on allele
  frequency, so bins SNPs by frequency and standardises within each bin.
  `ihs` says how unusual a SNP is *among SNPs at the same frequency* —
  which is why a warning about bins holding fewer than ten markers
  matters: the standardisation there has nothing to stand on.
- **`neg_log10_p` is a rank, not a test.** It is the normal tail of that
  z-score, assuming the bulk of the genome is neutral. Do not
  Bonferroni-correct it; use the empirical tail, `abs(ihs) > 4` or the
  top 1%.
- **The sign means nothing here.** With no outgroup the contrast is
  major-versus-minor, not ancestral-versus-derived. Read `abs(ihs)`.
- **Read runs, not points.** A real sweep is a *cluster* of elevated
  SNPs spanning tens of kb, since the whole haplotype is carried along.
  An isolated tall SNP beside ordinary neighbours is usually noise —
  which is what the Manhattan plot is for.

And one thing it cannot do: a sweep that went to **completion** leaves
no minor allele to contrast against, so iHS is blind to it. That is what
the cross-population scans are for.

``` r

plot_ihs(ihs, genes = PF_EXAMPLE_DRUG_GENES)
```

![](diversity-and-selection_files/figure-html/unnamed-chunk-11-1.png)

`zoom` crops any of these plots to one interval – a chromosome, a
`"chr:start-end"` range, or a gene name from `genes` – and names every
gene in the window in a track underneath:

``` r

plot_ihs(ihs, genes = PF_EXAMPLE_DRUG_GENES, zoom = "7", zoom_pad = 0)
```

![](diversity-and-selection_files/figure-html/unnamed-chunk-12-1.png)

### EHH around one SNP

A scan point says a locus looks selected; the EHH decay says what the
haplotype is doing there.
[`plot_ehh()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ehh.md)
draws one curve per allele **at the focal SNP** – so it is the
mutant-versus-reference comparison without needing the SNPs annotated –
and a sweep shows as one allele holding EHH near 1 far past where the
other has decayed:

``` r

plot_ehh(hap, "pfcrt", genes = PF_EXAMPLE_DRUG_GENES, span = 30000)
#> `pfcrt` holds several SNPs; measuring from Pf3D7_07_v3:403624 (minor allele 0.5) -- name a `chr:pos` to pick another
```

![](diversity-and-selection_files/figure-html/unnamed-chunk-13-1.png)

A marker that only segregates in one part of the cohort is buried by
pooling, so restrict the haplotypes to where it varies.
[`subset_haplotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/subset_haplotypes.md)
matches metadata the way `PopStructure$subset()` does,
`column = values`, and takes several values:

``` r

gh <- subset_haplotypes(hap, country = "Ghana")
plot_ehh(gh, "pfcrt", genes = PF_EXAMPLE_DRUG_GENES, span = 30000)
#> `pfcrt` holds several SNPs; measuring from Pf3D7_07_v3:403336 (minor allele 0.333) -- name a `chr:pos` to pick another
```

![](diversity-and-selection_files/figure-html/unnamed-chunk-14-1.png)

The SNP panel is left exactly as
[`parasite_haplotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/parasite_haplotypes.md)
built it, so two subsets remain comparable to each other and to the
whole; a SNP that is monomorphic within the subset is dropped when the
curve is computed anyway. [`print()`](https://rdrr.io/r/base/print.html)
reports the restriction, since it changes what every scan off that
object means.

`focal` also takes a bare position or a gene name from `genes`; a gene
holding several SNPs resolves to the most balanced one and says which,
since EHH from a near-singleton is a flat line. `group` gives a panel
per metadata group – worth doing after the pooled view, because EHH
knows nothing about population structure and a group with few carriers
gives a ragged curve, which is why `min_haplotypes` drops the smallest
ones. `span` crops the window, one value or two (`c(left =, right =)`).

### The genotypes behind a signal

A scan says a locus looks selected; the genotypes say what the
haplotypes actually are.
[`plot_region_haplotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_region_haplotypes.md)
draws one row per sample and one column per SNP over a window, clusters
the samples and puts the dendrogram beside them, with the gene track
underneath:

``` r

plot_region_haplotypes(ps, "pfcrt", pad = 20000, genes = PF_EXAMPLE_DRUG_GENES)
```

![](diversity-and-selection_files/figure-html/unnamed-chunk-15-1.png)

`split` blocks the rows by a metadata column and clusters *within* each
block, so a haplotype shared across a group reads as a solid band
instead of being scattered by one global ordering (this is what
`ComplexHeatmap`’s `row_split` does). `mark_snps` puts a line on the
SNPs worth pointing at:

``` r

plot_region_haplotypes(ps, "pfcrt", pad = 20000, split = "country",
                       genes = PF_EXAMPLE_DRUG_GENES, annotations = "country")
```

![](diversity-and-selection_files/figure-html/unnamed-chunk-16-1.png)

`annotations` adds coloured strips down the right, one column per
metadata variable, each with its own legend and taking its colours from
the object’s shared maps – so a level keeps the colour it has in the
UMAP or the admixture bars:

``` r

plot_region_haplotypes(ps, "pfdhps", pad = 20000, split = "country",
                       annotations = "country", genes = PF_EXAMPLE_DRUG_GENES)
```

![](diversity-and-selection_files/figure-html/unnamed-chunk-17-1.png)

`spacing` decides what the horizontal axis means, and the two answers
show different things. `"even"` (the default) gives every SNP the same
width, which is how the haplotype structure is easiest to read but says
nothing about distance. `"genomic"` keeps every mark the same width and
moves it to its real coordinate, so the gaps between SNPs are what you
see:

Either way a SNP is only ever drawn over the genes it really falls in.
Under `"even"` the axis counts SNPs, so a gene’s box is exactly the
columns it holds — its width is how many SNPs are in it, *not* how long
the gene is — and a gene with no genotyped SNP in the window has no
width at all and is left off the track with a message. Under `"genomic"`
the boxes are true extents and the marks are clipped at gene edges
instead.

``` r

plot_region_haplotypes(ps, "pfcrt", pad = 20000, split = "country",
                       spacing = "genomic", genes = PF_EXAMPLE_DRUG_GENES)
```

![](diversity-and-selection_files/figure-html/unnamed-chunk-18-1.png)

One thing to be sure of before reading the colours: `allele` says
whether a dosage of 2 means two alternate alleles or two reference ones.
The plot asks the object, which knows when it came from
\[load_genotypes()\]; an object built from a bare matrix cannot say, and
the two codings are indistinguishable afterwards, so it assumes alt
dosage and tells you it did.

Rsb and XP-EHH compare two populations rather than two alleles, so they
see exactly those completed sweeps:

``` r

head(run_rsb(hap, group = "country"), 3)
#> # A tibble: 3 × 8
#>   pair              pop1     pop2  chr            pos snp_id   value neg_log10_p
#>   <chr>             <chr>    <chr> <chr>        <dbl> <chr>    <dbl>       <dbl>
#> 1 Cambodia vs Ghana Cambodia Ghana Pf3D7_04_v3  92596 Pf3D7_0… -2.14        1.49
#> 2 Cambodia vs Ghana Cambodia Ghana Pf3D7_04_v3 401090 Pf3D7_0… -2.16        1.51
#> 3 Cambodia vs Ghana Cambodia Ghana Pf3D7_04_v3 544672 Pf3D7_0… -2.04        1.38
```

`pairs = list(c("a", "b"))` restricts the comparison; the default is
every pair.

[`ihs_genes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ihs_genes.md)
reduces any of these to one row per gene – the peak SNP inside it:

``` r

ihs_genes(ihs, genes = PF_EXAMPLE_DRUG_GENES, within = 10000, min_snps = 1)
#> # A tibble: 6 × 10
#>   group  gene  chr    start    end n_snps max_neg_log10_p max_abs_value peak_pos
#>   <fct>  <chr> <chr>  <dbl>  <dbl>  <int>           <dbl>         <dbl>    <dbl>
#> 1 Cambo… pfcrt 7     4.03e5 4.06e5     11           1.54          2.18    411005
#> 2 Cambo… pfdh… 8     5.48e5 5.51e5      9           0.448         0.922   546840
#> 3 Cambo… pfke… 13    1.72e6 1.73e6      1           0.132         0.335  1727144
#> 4 Ghana  pfke… 13    1.72e6 1.73e6     11           3.06          3.33   1734259
#> 5 Ghana  pfcrt 7     4.03e5 4.06e5     61           2.01          2.59    414656
#> 6 Ghana  pfdh… 8     5.48e5 5.51e5     18           0.454         0.932   548163
#> # ℹ 1 more variable: peak_in_gene <lgl>
```

### Beta

[`beta_score()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/beta_score.md)
looks for neighbourhoods where allele frequencies pile up around an
intermediate-frequency core – the footprint of long-term balancing
selection, and in this parasite usually an antigen rather than a drug
target.

``` r

b <- beta_score(ps, group = "country", window = 300000, min_window_snps = 1)
plot_beta(b)
```

![](diversity-and-selection_files/figure-html/unnamed-chunk-21-1.png)

Use `window = 1000` (the default) on a real, dense callset; the fixture
needs a far wider window just to find any neighbours.
[`beta_genes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/beta_genes.md)
summarises per gene, giving the per-gene table these scans are usually
reported as – mean beta beside the iHS peak and Tajima’s D from
`pop_diversity(by = "gene")`.

## Coverage QC

The depth tables come from the Python package
(`plasgenomicsutils coverage_depth_stats` and
`coverage_dropout_regions`, which read the BAMs); this package reads and
plots them.

``` r

cov <- read_coverage("coverage_by_sample.tsv.gz")
coverage_qc(cov, threshold = 10, min_mean = 5, min_breadth = 80)   # one row per sample
plot_coverage_summary(cov)     # mean vs breadth, failures labelled
plot_coverage_by_chrom(cov)    # sample x chromosome, relative to each sample's own mean

plot_coverage_dropout(read_coverage("coverage_windows.tsv.gz"),
                      genes = PF_EXAMPLE_DRUG_GENES)
```

Breadth matters more than mean depth: selective whole-genome
amplification can give a respectable average while leaving much of the
genome at zero, and only the `pct_ge_10x`-style column shows that. On
one real cohort the sample that failed QC had a mean of 122x — and a
median of 33x with a quarter of the core genome under 10x.

Two things to keep straight when generating those tables. The depth
engines do not agree by definition: `mosdepth` counts **fragments** (an
overlapping mate pair once) while `pysam`/`samtools depth` count
**reads**, which puts mosdepth a couple of percent lower everywhere.
Fragment depth is the better measure of independent evidence; either is
fine as long as a cohort uses one, and the `engine` column records
which. And
[`plot_coverage_dropout()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_coverage_dropout.md)
earns its place on the cross-sample question — a window empty in
*everyone* is not missing data, it silently reads as invariant.
