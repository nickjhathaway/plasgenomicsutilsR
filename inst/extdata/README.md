# Example IBD dataset

A small, **public** fixture for exercising the IBD plots, loaded by
`example_ibd_results()`.

**Provenance.** Derived from publicly available *P. falciparum* whole-genome
samples across five African countries (Tanzania, Kenya, DRC, Ethiopia, Sudan;
205 samples). Identity-by-descent was computed with the `plasgenomicsutils` IBD
tools and downsampled to ~1,357 SNPs for a compact fixture. Values are real (not
synthetic). Regions are the sample countries.

Files (all group-level; no per-sample genotype data). Here the grouping is the sample
country, so the `group` columns hold country names. `pos` is the variant's 0-based
position, matching the package convention (`?plasgenomicsutilsR-coordinates`):

- `ibd_per_snp_group.tsv.gz` — per-SNP fraction of pairs IBD, per group
  (`chr`, `pos`, `group`, `frac_pairs_ibd`)
- `ibd_pairwise_group.tsv.gz` — per-SNP IBD between group pairs
  (`chr`, `pos`, `group_a`, `group_b`, `frac_pairs_ibd`)
- `ibd_selection_per_group.tsv.gz` (+ `..._threshold.tsv`) — per-group IBD
  selection statistic
- `ibd_selection_global.tsv.gz` (+ `..._threshold.txt`) — genome-wide selection
  statistic
- `sample_regions.tsv` — the public sample → country map

## Example population-structure dataset

`pop_structure_ghana_cambodia.rds`, loaded by `example_pop_structure()`, is a small
**public** genotype matrix for the PCA / UMAP / sNMF-admixture tools.

**Provenance.** 60 publicly available *P. falciparum* samples (30 Ghana, 30 Cambodia).
Genotypes are real alt-allele dosages (`0`/`1`/`2`, `NA` where uncalled); the two countries
are strongly differentiated, so they separate on PC1 and sNMF picks K = 2.

**Two SNP panels**, because the two kinds of plot want opposite things:

- `genotype` — 49 biallelic SNPs thinned genome-wide, extracted from the Pf7 example BCF that
  also ships with the companion Python package
  (`plasgenomicsutils/tests/data/ghana_cambodia.pf7.tiny.bcf`). This is what PCA / UMAP /
  admixture read, and the panel those demos were tuned on.
- `full` — the same 49 plus **670** biallelic SNPs (MAF ≥ 2%) covering **pfcrt**, **pfdhps**
  and **pfkelch13** ±30 kb, pulled from the same public Ghana/Cambodia Pf7 callset. 719 SNPs
  in total. Locus-scale plots — `plot_ibd_locus()`, `zoom =`, `plot_region_haplotypes()`,
  `plot_ehh()` — need SNPs *close together* to show anything, and a genome-wide thinning at
  one SNP per few kb leaves 0–2 in a 50 kb window. Registered as the object's `"full"` panel,
  so those plots pick it up on their own while the thinned set stays the primary.

Ids are `chr:pos` with **1-based** VCF positions (SNPRelate's convention, which
`load_genotypes()` follows). Multiallelic sites split by `bcftools norm -m-` can leave two
biallelic records at one position; only the commonest ALT is kept, since a `chr:pos` id can
name just one SNP.

Contents: a list with `genotype` (60 × 49), `full` (60 × 719), `meta` (`sample`, `country`),
and `allele` / `pruned` recording what the dosages count and that the primary panel is
thinned.

`pop_structure_africa.rds` (`example_pop_structure("africa")`) is a richer, multi-region
**public** dataset for the combined UMAP + admixture figure (`plot_structure_figure()`):
258 published East/Central-African samples (`country`: DRC, Kenya, Tanzania, Uganda;
finer `site`: Kenya_East/West, Tanzania_East/West, three historical Uganda groups, DRC; a
macro `region`). Sample names and metadata come from a published supplemental table of
other-study samples; genotypes are extracted from that study's public combined BCF. To
make the regional structure clear in a compact fixture, the 2,000 SNPs kept are the ones
that **most differentiate the sites** — the top per-pair Jost's D markers
(`top_differentiating_snps()`), selected round-robin across all site pairs from ~20k
unpruned biallelic SNPs (an intentional, documented ascertainment for illustration).
Contents: a list of `genotype` (258 × 2000 integer matrix) and `meta` (`sample`,
`country`, `site`, `region`).

`pf3d7_drug_gene_cds.gff` is every feature (gene, mRNA, exon, CDS and UTR rows) of six
*P. falciparum* 3D7 drug-resistance genes, copied verbatim out of `PlasmoDB-55_Pfalciparum3D7.gff` (VEuPathDB): *pfdhfr*
`PF3D7_0417200`, *pfmdr1* `PF3D7_0523000`, *pfaat1* `PF3D7_0629500`, *pfcrt* `PF3D7_0709000`,
*pfdhps* `PF3D7_0810800` and *pfkelch13* `PF3D7_1343700`. Real coordinates, so codon intervals
computed from it are the published ones (*pfcrt* 76 at `Pf3D7_07_v3:403,624-403,626`).

The six were chosen to cover what the conversion has to get right: *pfcrt* has 13 coding
exons, four of whose codons straddle an intron; *pfkelch13* and *pfaat1* are on the minus
strand; and every gene carries markers quoted by residue in the literature. Used by
`read_gff_cds()` / `aa_intervals()` / `snp_aa_positions()` examples and the
*Amino acids and genomic coordinates* article, and by `read_gff_features()` (exons and
introns, for the *Avoiding tandem repeats* article).

`pf3d7_drug_gene_regions.fasta.gz` exists so the **tests** can check `snp_aa_positions()`'s
`ref_codon` / `ref_aa` against real published residues without touching the network. It is not
the documented way to supply sequence — the docs point `fasta =` at the released genome by URL,
which is what an analysis should use, since a fixture only ever covers the genes someone thought
to put in it. Each chromosome record keeps its **full length and true coordinates**, with real
bases across the span of the gene's CDS and `N` everywhere else, which costs almost nothing once
gzipped (35 KB for 23 Mb of nominal genome). Positions outside those spans return `NA` rather
than a base. Built from `Pf3D7.fasta` (PlasmoDB 3D7, 2015-06-18 assembly), matching the GFF
above.

The reference residues it yields are the published ones — *pfcrt* 76 `AAA`/K, *pfkelch13* 580
`TGT`/C, *pfdhfr* 108 `AGC`/S, *pfmdr1* 86 `AAT`/N. *pfdhps* 437 reads `GGT`/G, because 3D7 itself
carries the 437G allele the A437G marker is named for — the reference is not the wild type, here
or in general, so the test asserts `G` deliberately rather than treating it as a failure.

## Pf3D7 tandem repeats

`pf3d7_tandem_repeats.rds`, loaded by `pf3d7_tandem_repeats()`, is every short tandem repeat
found in the *P. falciparum* 3D7 reference (`Pf3D7.fasta`, PlasmoDB), so nobody has to run the
finders again. About 272,000 records over all 16 sequences, apicoplast and mitochondrion
included.

**Provenance.** The `combined.bed` of the HEOME redesign (April 2021): the union of Tandem
Repeats Finder, `trf Pf3D7.fasta 2 7 7 80 10 50 1000 -f -d -m -h` (match 2, mismatch 7, indel
7, match/indel probabilities 80/10, minimum score 50, maximum period 1000), converted with
`elucidator TandemRepeatFinderOutputToBed`, and `elucidator findSimpleTandemRepeatLocations
--maxRepeatUnitSize 10`, an exhaustive search for perfect repeats of units up to 10 bp that
`trf` is inconsistent about. Sorted and de-duplicated; the same span can appear twice with
different units where the two finders aligned it differently.

**Storage.** Only the first four BED fields carry information, and the name is
`chrom-start-end__UNIT_xCOPIES`, so the file holds `chrom` and `repeat_unit` as factors,
`start`/`end` as integers and `copies`, xz-compressed, and the name is rebuilt on load.
`copies` is stored only for the 4% of records (`trf` alignments with indels) where it is not
`width / unit_size` to the finder's six significant digits; the rest are `NA` in the file
and recomputed. 18 MB of BED becomes about 1.1 MB, and `data-raw/PF3D7_TANDEM_REPEATS.R`
checks the rebuilt names are byte-identical to the source. Coordinates are 0-based half-open
as in the BED.
