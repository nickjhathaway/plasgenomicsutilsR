# Documentation for the bundled Pf3D7 gene-coordinate datasets (built by
# data-raw/PF3D7_GENES.R). Lazy-loaded; see that script for provenance.

#' Pf3D7 gene coordinates
#'
#' Every protein-coding gene in the *Plasmodium falciparum* 3D7 reference, for use as
#' the `genes` track of an [IbdResults] object (gene reference lines on Manhattan /
#' tug-of-war plots, the pairwise-IBD triangles, and [pos_selection_genes()]).
#'
#' Coordinates are the **transcript (mRNA) span** -- the minimum transcript start to the
#' maximum transcript end across a gene's isoforms -- rather than the wider `gene`
#' feature, so they exclude the untranslated flanks recorded on the gene record.
#'
#' `name` is a friendly `"pf"`-prefixed lower-casing of the GFF `Name` attribute
#' (e.g. `CRT` -> `pfcrt`, `MSP1` -> `pfmsp1`); genes with no `Name` fall back to their
#' `gene_id` (so `name` equals `gene_id` for those). The two folate genes, whose GFF
#' Names are long compounds, are surgically set to `pfdhfr` (`PF3D7_0417200`) and
#' `pfdhps` (`PF3D7_0810800`). Family names (`pfvar`, `pfrif`, ...) repeat across the
#' family's members, so `name` is not unique; `gene_id` is.
#'
#' @format A data frame with `r nrow(plasgenomicsutilsR::PF3D7_GENES)` rows and columns:
#' \describe{
#'   \item{Pf3D7_chrom}{sequence id, e.g. `Pf3D7_01_v3`}
#'   \item{start, end}{CDS-span bounds, 0-based half-open `[start, end)` -- the translated
#'     extent (min CDS start to max CDS end across isoforms), excluding UTRs}
#'   \item{chrom}{short chromosome (`"1"`..`"14"`, or `"API"` / `"MIT"`)}
#'   \item{gene_id}{PlasmoDB gene id, e.g. `PF3D7_0709000`}
#'   \item{name}{friendly display name (see above)}
#' }
#' @source VEuPathDB / PlasmoDB Pf3D7 GFF, release version 2020-09-01. The GFF is 1-based
#'   inclusive; `start` is shifted down by one when the dataset is built so the shipped
#'   coordinates follow the package convention (see [plasgenomicsutilsR-coordinates]).
#' @seealso [PF_EXAMPLE_DRUG_GENES], [pos_selection_genes()]
"PF3D7_GENES"

#' Pf3D7 drug-resistance / selection example genes
#'
#' A small, curated subset of [PF3D7_GENES]: well-known *P. falciparum*
#' drug-resistance and selection loci (`pfcrt`, `pfdhfr`, `pfmdr1`, `pfdhps`,
#' `pfkelch13`, `pfaat1`, `pfgch1`, `pfpx1`). Used as the default `genes` track of
#' [example_ibd_results()] and throughout the documentation.
#'
#' @format A data frame with 8 rows and the same columns as [PF3D7_GENES].
#' @source VEuPathDB / PlasmoDB Pf3D7 GFF, release version 2020-09-01.
#' @seealso [PF3D7_GENES]
"PF_EXAMPLE_DRUG_GENES"

#' Pf3D7 core genome regions
#'
#' The **core** (non-subtelomeric, non-hypervariable) intervals of the *P. falciparum* 3D7
#' core chromosomes. A locus is "core" when it overlaps one of these intervals and
#' "subtelomeric / internally hypervariable" otherwise -- intersect a gene or SNP table
#' against this with [bed_intersect()].
#'
#' @format A data frame with columns `Pf3D7_chrom`, `start`, `end`, `chrom`; coordinates are
#'   0-based half-open (see [plasgenomicsutilsR-coordinates]).
#' @source MalariaGEN `regions-20130225` core/non-core boundaries (with a few core
#'   extensions bringing `hrp2` / `hrp3` into core), as shipped by the companion Python
#'   package (`plasgenomicsutils` `builtin:pf3d7_core_regions`).
#' @seealso [PF3D7_PARALOG_GENES], [bed_intersect()]
"PF3D7_CORE_REGIONS"

#' Pf3D7 paralogous / hypervariable gene families
#'
#' Members of the *P. falciparum* multi-gene / hypervariable families (var, rifin, stevor,
#' surfin, ...) whose short reads mismap; commonly masked in population-genetic analyses.
#' Intersect against this with [bed_intersect()] to flag genes that fall in these families.
#'
#' @format A data frame with columns `Pf3D7_chrom`, `start`, `end`, `chrom`, `gene_id`,
#'   `description`; coordinates are 0-based half-open (see
#'   [plasgenomicsutilsR-coordinates]).
#' @source The companion Python package (`plasgenomicsutils` `builtin:pf3d7_paralog_genes`).
#' @seealso [PF3D7_CORE_REGIONS], [bed_intersect()]
"PF3D7_PARALOG_GENES"
