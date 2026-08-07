# Pf3D7 gene coordinates

Every protein-coding gene in the *Plasmodium falciparum* 3D7 reference,
for use as the `genes` track of an
[IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
object (gene reference lines on Manhattan / tug-of-war plots, the
pairwise-IBD triangles, and
[`pos_selection_genes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pos_selection_genes.md)).

## Usage

``` r
PF3D7_GENES
```

## Format

A data frame with 5318 rows and columns:

- Pf3D7_chrom:

  sequence id, e.g. `Pf3D7_01_v3`

- start, end:

  CDS-span bounds, 0-based half-open `[start, end)` – the translated
  extent (min CDS start to max CDS end across isoforms), excluding UTRs

- chrom:

  short chromosome (`"1"`..`"14"`, or `"API"` / `"MIT"`)

- gene_id:

  PlasmoDB gene id, e.g. `PF3D7_0709000`

- name:

  friendly display name (see above)

## Source

VEuPathDB / PlasmoDB Pf3D7 GFF, release version 2020-09-01. The GFF is
1-based inclusive; `start` is shifted down by one when the dataset is
built so the shipped coordinates follow the package convention (see
[plasgenomicsutilsR-coordinates](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plasgenomicsutilsR-coordinates.md)).

## Details

Coordinates are the **transcript (mRNA) span** – the minimum transcript
start to the maximum transcript end across a gene's isoforms – rather
than the wider `gene` feature, so they exclude the untranslated flanks
recorded on the gene record.

`name` is a friendly `"pf"`-prefixed lower-casing of the GFF `Name`
attribute (e.g. `CRT` -\> `pfcrt`, `MSP1` -\> `pfmsp1`); genes with no
`Name` fall back to their `gene_id` (so `name` equals `gene_id` for
those). The two folate genes, whose GFF Names are long compounds, are
surgically set to `pfdhfr` (`PF3D7_0417200`) and `pfdhps`
(`PF3D7_0810800`). Family names (`pfvar`, `pfrif`, ...) repeat across
the family's members, so `name` is not unique; `gene_id` is.

## See also

[PF_EXAMPLE_DRUG_GENES](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF_EXAMPLE_DRUG_GENES.md),
[`pos_selection_genes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pos_selection_genes.md)
