# Load genotypes from a VCF, optionally LD-pruned

Converts a VCF to GDS (only when needed) and returns the genotype
matrix, LD-pruned by default. The backend is SeqArray, whose SeqVarGDS
stores allele **indices** and so can hold a record with any number of
alleles.

## Usage

``` r
load_genotypes(
  vcf,
  gds = NULL,
  prune = NULL,
  ld_threshold = 0.2,
  slide_max_bp = 20000,
  slide_max_n = 200,
  autosome_only = FALSE,
  maf = NaN,
  missing_rate = NaN,
  seed = 42,
  vcf_dir = NULL,
  allele = c("alt", "ref"),
  variants = c("biallelic_snvs", "all"),
  encoding = c("dosage", "allele_index", "allele_set"),
  star = c("missing", "allele"),
  refresh = c("stale", "never", "always")
)
```

## Arguments

- vcf:

  Path to a (bgzipped) VCF, or a **BCF** – converted to VCF text first
  with `bcftools`, reusing any VCF already sitting next to it rather
  than making another copy.

- gds:

  Optional GDS path; derived from `vcf` if `NULL`.

- prune:

  LD-prune. Defaults to `TRUE` for `encoding = "dosage"` and `FALSE`
  otherwise, because pruning measures correlation between dosages and an
  allele index is not a count of anything. `FALSE` returns every record
  `variants` admits (see *Which records reach the panel*), unpruned –
  use this for the genotype matrix fed to
  [`pop_diff()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff.md)
  /
  [`pop_diff_table()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff_table.md),
  since LD-pruning removes the very SNPs that carry the differentiation
  signal.

- ld_threshold, slide_max_bp, slide_max_n, autosome_only:

  Passed to
  [`SNPRelate::snpgdsLDpruning()`](https://rdrr.io/pkg/SNPRelate/man/snpgdsLDpruning.html)
  (defaults 0.2 / 20000 / 200 / `FALSE`); ignored when `prune = FALSE`.

- maf, missing_rate:

  Optional MAF / per-SNP missing-rate cutoffs for pruning.

- seed:

  Random seed for the pruning.

- vcf_dir:

  Where to put the VCF converted from a BCF (default: alongside the
  BCF). Point it somewhere scratch to keep converted copies out of the
  data directory.

- allele:

  Which allele the returned dosage counts, `"alt"` (default) or `"ref"`.
  Only reported allele frequencies (and the arbitrary sign of a PCA
  axis) depend on this – every diversity, differentiation, LD and
  selection statistic here is symmetric in `p` and `1 - p`. Meaningless
  for `encoding = "allele_index"`, which names an allele rather than
  counting copies of one, and refused there.

- variants:

  Which records to read. `"biallelic_snvs"` (default) keeps biallelic
  SNVs only; `"all"` keeps every record, multiallelic sites and indels
  included.

- encoding:

  How `genotype` is coded. `"dosage"` (default) counts alternate copies,
  0/1/2 – the classical matrix, and what PCA, admixture, diversity and
  differentiation read. `"allele_index"` names **which** allele each
  sample carries, 0 for the reference and 1..k for the alternates, which
  can carry a multiallelic site but sets a mixed call to missing.
  `"allele_set"` names the **set** of alleles a sample carries – a
  singleton set for a pure call, a two-element set for a mix – so it
  keeps both a mixed call and its multiallelic identity. It is a strict
  superset of the other two: a
  [PopStructure](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PopStructure.md)
  derives an `allele_index` (or, through it, a `dosage`) view from it on
  demand (`$genotype(needs = )`), so one stored panel serves both the
  haplotype plot and every analysis. See *Which encoding to ask for*.

- star:

  How to treat the `*` spanning-deletion allele: `"missing"` (default)
  records it as `has_spanning_del` and leaves it out of the allele
  count, so a call naming it becomes missing; `"allele"` keeps it as a
  literal allele.

- refresh:

  What to do with a derived file older than what it was built from – the
  text VCF beside a BCF, and the GDS beside either. `"stale"` (default)
  rebuilds it, `"always"` rebuilds regardless, `"never"` reuses it and
  warns. The default is a rebuild because the alternative is silent:
  update the BCF and every result below comes from the old records, with
  the GDS looking current because it is newer than the stale VCF it was
  built from. A GDS left over from the older SNPRelate backend is a
  different format entirely and is rebuilt whatever `refresh` says.

## Value

A list with `genotype` (matrix; sample row names and `chr:pos0` column
names – 0-based, like every other position in the package), `sample.id`,
`snp.id`, and the facts the matrix itself cannot carry: `allele` (which
allele a dosage counts), `encoding`, `pruned`, `positions` and
`variants`.

For `encoding = "allele_set"` the codes in `genotype` are a 0-based
index into per-column state lists, so two more fields come back:
`state_levels` (named by column: the state names a code resolves to,
e.g. `reference` / `alternate 1` / `alternate 1 + alternate 2`) and
`state_sets` (the same states as integer allele vectors, which is what
lets an allele index be derived from the set).
[PopStructure](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PopStructure.md)
stores both.
[PopStructure](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PopStructure.md)
keeps them, so anything that names a call, refuses a non-dosage panel or
warns about pruning can ask instead of assuming.

Also `sites`: one row per genotype column, in the same order, with
`site_key`, `chr`, `pos` (0-based), `ref`, `alt`, `n_alt`, `n_alt_real`
and `has_spanning_del`.

`alt` is a list column holding the record's **own** ALT list, in its own
order and `*` included, because that is what an allele index names –
`alt[[i]][k]` is the allele a genotype of `k` refers to, and dropping
anything from it would point the index at the wrong base. `n_alt` counts
that list for the same reason.

Whether a record is *multiallelic* is a different question, since `*` is
a missingness annotation rather than an allele: `n_alt_real` excludes it
and is the one that answers, with `has_spanning_del` recording the `*`
separately. So `A > *,T` is an ordinary biallelic SNV with `n_alt = 2`
and `n_alt_real = 1`.

**Ask this table before reading a per-allele number off a dosage
matrix**: a column with `n_alt_real > 1` has had its alternates
collapsed onto one number and cannot answer one.

## Details

Pruned or not depends on the question. Pruning is right for PCA, UMAP
and admixture, where correlated SNPs would let one locus dominate the
structure. It is wrong wherever the correlation between neighbouring
SNPs *is* the signal – differentiation
([`pop_diff()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff.md))
and haplotypes
([`plot_region_haplotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_region_haplotypes.md))
– because it keeps one SNP out of each correlated run and drops the
rest. Holding both is cheap: the GDS is reused, so a second call with
`prune = FALSE` only re-reads it.

## Which records reach the panel

`prune = FALSE` means unpruned, not every record: the default
`variants = "biallelic_snvs"` keeps only records whose single ALT is a
substitution, so the panel is usually smaller than the VCF's record
count. Three kinds are skipped, and the message says **how many of
each**:

- sites with **no ALT allele** (`ALT="."`) – the reference positions an
  all-sites caller emits. Usually the biggest share by far, and the
  easiest to miss, since nothing about them says "variant":
  `bcftools view --exclude-types indels` leaves every one in place,
  because they are not indels.

- **indels** and other non-SNV records.

- sites with **more than one ALT** – which `variants = "all"` keeps.

The tally is exact and costs nothing: SeqArray has already read every
record by the time the selection happens, so the classification that
chose the panel is the one reported. A record carrying a `*` is its own
class and is skipped: part of the cohort has no base there to compare,
so `A > *,T` is not a clean SNP site however well T behaves. The
companion Python package's `spanning_del_filter` is the way to keep such
a site – it recodes the deleted calls as missing and drops the allele,
after which the record really is an ordinary SNP.

Sites that are **invariant across the loaded samples are kept**, as long
as the record lists an ALT. A site every sample calls `1/1`, or every
sample calls `0/0`, comes through – which is why the functions needing
variable sites
([`parasite_haplotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/parasite_haplotypes.md),
[`run_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_ihs.md))
apply their own `maf` cutoff instead of trusting the panel.

## Keeping every variant

`variants = "all"` keeps every record – multiallelic sites, indels and
`ALT="."` positions included. Pair it with `encoding = "allele_index"`
to keep the alternates apart; with the default `encoding = "dosage"` it
warns, because a dosage cannot say which alternate a non-reference call
carries:

    encoding = "dosage"                  encoding = "allele_index"
      pos  alleles  s1 s2 s3               pos  alleles  s1 s2 s3
      400  C/T,G     0  2  2               400  C/T,G     0  1  2
                        ^  ^                                  ^  ^
                        1/1 and 2/2                           T and G, told apart
                        land on one number

One GDS serves every combination: SeqArray reads the whole callset and
the selection happens in R, so asking for a different `variants` or
`encoding` re-reads the same file rather than rebuilding it.

## Which encoding to ask for

A **dosage** says how many copies of one allele a sample has. That is
two numbers, and a site with three alleles needs three, so a dosage
cannot say *which* alternate a non-reference call carries – at a codon
where D384A, D384G and D384Y arose separately, every carrier lands on
the same number. Ask for it when the panel is biallelic and the consumer
is PCA, UMAP, admixture, diversity or differentiation.

An **allele index** names the allele: 0 for the reference, 1..k for the
alternates in the record's own ALT order (`$sites$alt` says which is
which). It carries a multiallelic site faithfully, and the dosage
statistics refuse it rather than mangling it – there is no faithful
conversion, so contrast one allele at a time instead.

The parasite is haploid, so a heterozygous call is a **mixed infection**
rather than a diploid genotype, and one index cannot name two clones.
Those calls are missing under `allele_index`, which is the same reading
`pop_diversity(het = "missing")` takes.

## See also

[PopStructure](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PopStructure.md),
[`pop_structure()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_structure.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# LD-pruned for PCA / UMAP / admixture
geno <- load_genotypes("clean.vcf.gz", gds = "clean.gds")

# and the full panel for anything where SNP correlation IS the signal
full <- load_genotypes("clean.vcf.gz", gds = "clean.gds", prune = FALSE)
} # }
```
