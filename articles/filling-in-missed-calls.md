# Filling in missed calls

``` r

library(plasgenomicsutilsR)
cds <- read_gff_cds(system.file("extdata", "pf3d7_drug_gene_cds.gff",
                                package = "plasgenomicsutilsR"))
ps <- example_pop_structure(umap = FALSE)
```

Filtering removes SNPs, and some of what it removes is a residue you
meant to keep. Going the other way finds them: take the codons you care
about, subtract what the panel already holds, and what is left is what
was never called.

[`bed_subtract()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/bed_subtract.md)
cuts at base resolution, which is the part that matters. A codon with
one of its three bases in the panel comes back as the other two – a
residue cannot be read from one base, so the partly covered codon is
exactly the one worth re-calling.

``` r

want <- aa_intervals(data.frame(transcript_id = c("pfcrt", "pfcrt", "pfkelch13"),
                                aa_position = c(72, 76, 580)), cds)

# a panel that happens to hold only the middle base of the K76T codon
have <- paste0("Pf3D7_07_v3:", want$start[want$aa_position == 76] + 1)

gaps <- bed_subtract(want, have)
gaps[, c("name", "aa_position", "start", "end", "piece", "width")]
#> # A tibble: 4 × 6
#>   name                  aa_position   start     end piece width
#>   <chr>                       <int>   <dbl>   <dbl> <int> <dbl>
#> 1 PF3D7_0709000.1-AA72           72  403611  403614     1     3
#> 2 PF3D7_0709000.1-AA76           76  403623  403624     1     1
#> 3 PF3D7_0709000.1-AA76           76  403625  403626     2     1
#> 4 PF3D7_1343700.1-AA580         580 1725257 1725260     1     3
```

`locs2` can be named by whatever holds the panel – `chr:pos` ids as
above, a genotype matrix, a
[`load_genotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/load_genotypes.md)
list, or a `PopStructure`:

``` r

sum(bed_subtract(want, ps)$width)   # bases of those three codons the panel lacks
#> [1] 9
```

[`write_bed()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/write_bed.md)
writes the result for `bcftools mpileup -R`:

``` r

write_bed(gaps, "fill_in_snps.bed")
```

``` bash
bcftools mpileup -f Pf3D7.fasta -R fill_in_snps.bed --bam-list bams.txt   | bcftools call --ploidy 2 -m -Ob -o fill_in_calls.bcf
```

It writes `chrom` in preference to `chr`. Tables here carry both – `chr`
normalised for matching (`"7"`) and `chrom` as the source spells it
(`"Pf3D7_07_v3"`) – and a BED is read by other tools against a real
reference, so the normalised name produces a file that matches nothing,
silently.

One caveat on intron-spanning codons:
[`aa_intervals()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/aa_intervals.md)
reports `[start, end)`, which for a codon split by an intron reaches
across the intron. The exact bases are in `codon_positions`, and
`spans_intron` flags the rows where the two differ.
