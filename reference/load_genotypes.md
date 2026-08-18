# Load genotypes from a VCF, optionally LD-pruned

Converts a VCF to GDS (only when needed) and returns the genotype matrix
(samples x SNPs, coded 0/1/2, `NA` for missing) via SNPRelate, LD-pruned
by default.

## Usage

``` r
load_genotypes(
  vcf,
  gds = NULL,
  prune = TRUE,
  ld_threshold = 0.2,
  slide_max_bp = 20000,
  slide_max_n = 200,
  autosome_only = FALSE,
  maf = NaN,
  missing_rate = NaN,
  seed = 42,
  vcf_dir = NULL,
  allele = c("alt", "ref"),
  variants = c("biallelic_snvs", "all")
)
```

## Arguments

- vcf:

  Path to a (bgzipped) VCF, or a **BCF** – SNPRelate reads VCF text
  only, so a BCF is converted first with `bcftools`, reusing any VCF
  already sitting next to it rather than making another copy.

- gds:

  Optional GDS path; derived from `vcf` if `NULL`.

- prune:

  LD-prune (default `TRUE`). `FALSE` returns every record `variants`
  admits (see *Which records reach the panel*), unpruned – use this for
  the genotype matrix fed to
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

  Which allele the returned dosage counts. SNPRelate counts the
  **reference** allele; the default `"alt"` flips that so the matrix
  means what the rest of the package says it means. Only reported allele
  frequencies (and the arbitrary sign of a PCA axis) depend on this –
  every diversity, differentiation, LD and selection statistic here is
  symmetric in `p` and `1 - p`.

- variants:

  Which records to read. `"biallelic_snvs"` (default) keeps biallelic
  SNVs only; `"all"` keeps every record, multiallelic sites and indels
  included, at the cost of a dosage that cannot say which ALT it counts.
  See the two sections below – nothing in this package handles `"all"`,
  and it warns.

## Value

A list with `genotype` (matrix; sample row names and `chr:pos0` column
names – 0-based, like every other position in the package), `sample.id`,
`snp.id`, and the facts the matrix itself cannot carry: `allele` (which
allele the dosages count), `pruned`, `positions` and `variants`.
[PopStructure](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PopStructure.md)
keeps them, so anything that names a call or warns about pruning can ask
instead of assuming.

## Details

Which you want depends on the question. Pruning is right for PCA, UMAP
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

`prune = FALSE` means unpruned, not every record: the panel is usually
smaller than the VCF's record count whatever `prune` is, because the
default `variants = "biallelic_snvs"` has SNPRelate read the file with
`method = "biallelic.only"`. Three kinds of record are skipped,
silently:

- sites with **no ALT allele** (`ALT="."`) – the reference positions an
  all-sites caller emits. Usually the biggest share by far, and the
  easiest to miss, since nothing about them says "variant":
  `bcftools view --exclude-types indels` leaves every one of them in
  place, because they are not indels.

- **indels** and other non-SNV records.

- sites with **more than one ALT**.

So a VCF of 28,927 records carrying 9,158 `ALT="."` positions and 48
multiallelic sites loads as 19,721 SNPs. When a panel comes out short,
count what is actually there rather than the total, and drop the no-ALT
records upstream if you would rather the two numbers agree:

    bcftools view -H -m2 -M2 -v snps file.vcf.gz | wc -l   # what will load
    bcftools view -e 'ALT="."' -Ob -o out.bcf in.bcf       # or --min-ac 1

Sites that are **invariant across the loaded samples are kept**, as long
as the VCF lists an ALT allele there: `"biallelic.only"` asks how many
alleles the record declares, not whether these samples differ. A site
every sample calls `1/1`, or every sample calls `0/0`, comes through –
which is why the functions needing variable sites
([`parasite_haplotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/parasite_haplotypes.md),
[`run_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_ihs.md))
apply their own `maf` cutoff instead of trusting the panel. SNPRelate's
own help calls this "excluding monomorphic variants", which is easy to
read as the stronger promise.

## Keeping every variant

`variants = "all"` reads the VCF with `method = "copy.num.of.ref"`
instead, which keeps every record – multiallelic sites, indels and
`ALT="."` positions included – and stores the **copy number of the
reference allele**. Nothing in this package reads that correctly, so it
warns; the switch is here to hand the matrix to something that does.

What breaks is the coding, not the reading. A dosage says how many
reference copies a sample has and nothing about *which* alternate allele
makes up the rest, so at `C -> T,G` a sample called `1/1` and a sample
called `2/2` are both 0 reference copies and land on the same number
despite carrying different alleles:

    variants = "biallelic_snvs"        variants = "all"
      pos  allele  s1 s2 s3              pos  allele  s1 s2 s3
      100  A/G      2  0  1              100  A/G      2  0  1
      500  G/A      0  0  0              200  T/.      2  2  2   <- no ALT
                                         300  AT/A     2  0  2   <- indel
                                         400  C/T,G    2  0  0   <- 1/1 and 2/2 both 0
                                         500  G/A      0  0  0

The allele strings themselves survive in the GDS
(`SNPRelate::snpgdsSNPList()$allele` gives `"C/T,G"`), so which alleles
exist is recoverable even though the dosage cannot express them. The
returned list records the choice as `variants`, and the GDS is tagged
with it, so the two panels never get confused for one another through a
reused `.gds`.

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
