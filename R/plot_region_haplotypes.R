# Genotypes over one region, sample by SNP, with the samples clustered.

.GENO_LEVELS <- c("reference", "mixed", "alternate")
.HAP_ROWS_IN <- 9        # inches the rows get when they are not individually labelled
.HAP_ANN_IN <- 0.22      # inches per annotation strip
# Every call gets a colour of its own, from the package's colour-blind-friendly set, and only
# missing data is grey: a near-white fill for one of the calls is hard to tell from both the
# panel and the grey of a missing call, which is the one distinction that must stay obvious.
# The genotype legend sits above the annotation legends, which then follow the order the caller
# listed them in. Without an explicit order ggplot sorts guides by an internal hash of their
# labels, so the stack rearranged itself between datasets and two figures stopped being
# comparable.
.HAP_LEGEND_CALL <- 1L

.GENO_FILL <- c(reference = "#2271B2", mixed = "#359B73", alternate = "#D55E00")


# A marker's calls as *allele sets*, read from the callset rather than from a dosage matrix.
# A dosage says how many copies of one allele a sample carries, which at a site with three
# alleles cannot say which of them it is; the set can. States are named from the alleles they
# contain, and only the states that occur are returned -- the full enumeration of a
# triallelic site is seven, and most of them are usually empty.
.read_genotype_sets <- function(vcf, names = c("index", "base")) {
  names <- match.arg(names)
  if (!nzchar(Sys.which("bcftools")))
    stop("reading allele sets needs bcftools on PATH", call. = FALSE)
  if (!file.exists(vcf)) stop("no such file: ", vcf, call. = FALSE)
  samples <- system2("bcftools", c("query", "-l", shQuote(vcf)), stdout = TRUE, stderr = FALSE)
  # REF comes back in the same pass: fetching it per record with `-r` would need the callset
  # to be indexed, which a hand-built VCF need not be.
  lines <- system2("bcftools", c("query", "-f",
                                 shQuote("%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n"),
                                 shQuote(vcf)), stdout = TRUE, stderr = FALSE)
  if (!length(lines)) stop("no records in ", basename(vcf), call. = FALSE)
  parts <- strsplit(lines, "\t", fixed = TRUE)

  codes <- matrix(NA_integer_, nrow = length(samples), ncol = length(parts),
                  dimnames = list(samples, NULL))
  levs <- vector("list", length(parts))
  ids <- character(length(parts))
  for (j in seq_along(parts)) {
    chrom <- parts[[j]][1]
    pos <- as.integer(parts[[j]][2])
    ref <- parts[[j]][3]
    alts <- strsplit(parts[[j]][4], ",", fixed = TRUE)[[1]]
    n_alleles <- 1L + length(alts)
    gt <- parts[[j]][-(1:4)]
    sets <- lapply(strsplit(gt, "[/|]"), function(v) {
      v <- suppressWarnings(as.integer(v[v != "."]))
      if (!length(v)) NA_integer_ else sort(unique(v))
    })
    key <- vapply(sets, function(v) if (length(v) == 1 && is.na(v[1])) NA_character_
                  else paste(v, collapse = ","), character(1))
    seen <- sort(unique(key[!is.na(key)]))
    nm <- vapply(seen, function(k) .allele_set_name(
      as.integer(strsplit(k, ",", fixed = TRUE)[[1]]), n_alleles,
      alleles = if (identical(names, "base")) c(ref, alts) else NULL), character(1))
    codes[, j] <- match(key, seen) - 1L
    levs[[j]] <- unname(nm)
    ids[j] <- paste0(chrom, ":", pos - 1L)      # 0-based, as every id in this package is
  }
  colnames(codes) <- ids
  names(levs) <- ids
  list(codes = codes, levels = levs)
}

# Colours for the extra call states, chosen to stay apart from the three the plot already
# uses. Taking the next entries off the shared palette is not enough: it hands back an orange
# for `alternate 1` that sits next to the existing `alternate`, and a triallelic column then
# reads as an ordinary one. Pick greedily on worst-case CIEDE2000 across all three
# dichromacies instead, which is the same measure `colour_blind_distance()` reports.
.distinct_fills <- function(used, n) {
  if (n <= 0) return(character(0))
  cand <- setdiff(color_palette(min(12L, max(8L, n + 5L))), unname(used))
  picked <- character(0)
  for (i in seq_len(n)) {
    if (!length(cand)) break
    score <- vapply(cand, function(cc)
      min(colour_blind_distance(c(unname(used), picked), cc), na.rm = TRUE), numeric(1))
    best <- cand[which.max(score)]
    picked <- c(picked, best)
    cand <- setdiff(cand, best)
  }
  # a palette that ran out is better short than recycled into a duplicate
  picked
}

# Attach the extra fills, saying so when there are not enough. Two distinct allele states
# sharing a colour is exactly the confusion the greedy pick above exists to avoid, and it
# fails silently -- the plot looks perfectly fine. So the states past the end of the palette
# get no fill (ggplot draws them grey and lists them in the legend) and the caller is told.
.assign_extra_fills <- function(fills, extra) {
  if (!length(extra)) return(fills)
  got <- .distinct_fills(fills, length(extra))
  if (length(got) < length(extra)) {
    warning(sprintf(paste("%d call state(s) need a fill and the colour-blind-safe palette",
                          "offers %d that stay clear of the ones already in use, so %s",
                          "%s no colour of %s own. Narrow the window, or pass `colours =`",
                          "with a scale you have checked."),
                    length(extra), length(got),
                    paste(utils::head(extra[-seq_along(got)], 3), collapse = ", "),
                    if (length(extra) - length(got) == 1L) "has" else "have",
                    if (length(extra) - length(got) == 1L) "its" else "their"),
            call. = FALSE)
  }
  c(fills, stats::setNames(got, extra[seq_along(got)]))
}

# Row-clustering distance that does not put an ordered scale on a nominal column.
#
# An `additional_genotypes` marker's codes are an arbitrary sorted index into that marker's
# own states: `alternate 2` is not twice as far from `reference` as `alternate 1` is, and for
# alleles of independent origin no ordering exists at all. Euclidean distance on those codes
# also lets one triallelic marker carry up to 4 units of distance where a biallelic SNP
# carries 2, so it outweighs several SNPs in the Ward ordering that decides row order.
#
# Dosage columns keep the ordinary squared difference, which is meaningful there. Nominal
# columns contribute a mismatch indicator scaled to the same range, so one disagreement at a
# nominal marker weighs the same as one homozygous difference at a SNP.
.geno_dist <- function(G, nominal_ids) {
  nominal <- colnames(G) %in% nominal_ids
  if (!any(nominal)) return(stats::dist(G))
  num <- G[, !nominal, drop = FALSE]
  nom <- G[, nominal, drop = FALSE]
  d2 <- if (ncol(num)) as.matrix(stats::dist(num))^2 else
    matrix(0, nrow(G), nrow(G), dimnames = list(rownames(G), rownames(G)))
  for (j in seq_len(ncol(nom))) {
    v <- nom[, j]
    mism <- outer(v, v, function(a, b) as.numeric(a != b))
    mism[is.na(mism)] <- NA_real_
    # 2 to match a homozygous difference on the 0/1/2 dosage scale
    d2 <- d2 + ifelse(is.na(mism), 0, mism) * 4
  }
  stats::as.dist(sqrt(d2))
}

# What to call the set of alleles a sample carries at one marker. A biallelic marker keeps
# the wording the plot has always used, so adding a multiallelic marker beside biallelic ones
# does not rename the calls they were already showing.
.allele_set_name <- function(a, n_alleles, alleles = NULL) {
  # With `alleles` the state is named by the bases themselves -- what a carrier contrast
  # wants, since `carrier = "C"` names something the reader can check against the callset.
  # Without them it is the positional wording a legend wants, and a biallelic marker keeps
  # exactly the words the plot has always used.
  if (!is.null(alleles)) {
    idx <- a + 1L
    idx[idx < 1L | idx > length(alleles)] <- NA_integer_
    return(paste(alleles[idx], collapse = " + "))
  }
  one <- function(i) if (i == 0L) "reference"
                     else if (n_alleles <= 2) "alternate" else paste("alternate", i)
  if (length(a) == 1L) return(one(a))
  if (n_alleles <= 2 && identical(a, 0:1)) return("mixed")
  paste(vapply(a, one, character(1)), collapse = " + ")
}

# What a dosage means depends on which allele it counts, and the two codings are
# indistinguishable from the matrix alone -- 2 is homozygous alternate under alt dosage and
# homozygous reference under ref dosage. Getting it backwards silently mislabels the whole
# plot, so the object is asked rather than assumed.
.geno_calls <- function(v, allele, snp_id = NULL, state_levels = NULL) {
  out <- rep(NA_character_, length(v))

  # A marker read as allele sets carries its own states, and its `value` is a *nominal*
  # index into them -- not a dosage. Resolve those rows first and mark them, so the dosage
  # arithmetic below never sees them: `3L - v` on a 4-state column yields a 0 subscript,
  # which R drops silently and shortens the result, and on a 5-state column it yields a
  # negative one, which errors outright.
  is_state <- rep(FALSE, length(v))
  for (id in names(state_levels)) {
    hit <- which(snp_id == id)
    if (length(hit)) {
      lv <- state_levels[[id]]
      k <- v[hit] + 1L
      k[!is.na(k) & (k < 1L | k > length(lv))] <- NA_integer_   # unknown state, not a drop
      out[hit] <- lv[k]
      is_state[hit] <- TRUE
    }
  }

  d <- which(!is_state & !is.na(v))
  if (length(d)) {
    idx <- if (identical(allele, "ref")) 3L - v[d] else v[d] + 1L
    idx[idx < 1L | idx > length(.GENO_LEVELS)] <- NA_integer_
    out[d] <- .GENO_LEVELS[idx]
  }

  extra <- setdiff(unlist(state_levels, use.names = FALSE), .GENO_LEVELS)
  factor(out, levels = c(.GENO_LEVELS, extra))
}

.resolve_allele <- function(x, allele) {
  if (!is.null(allele)) return(match.arg(allele, c("alt", "ref")))
  recorded <- if (is.function(x$allele)) x$allele() else NULL
  if (!is.null(recorded)) return(recorded)
  message("this object does not record which allele its dosages count; assuming `alt` -- ",
          "pass `allele = \"ref\"` if it counts reference alleles instead")
  "alt"
}

# Samples in clustered order within each split block. Clustering per block rather than over
# everything is the point of splitting: the blocks are fixed by the annotation and the
# ordering inside each one is still learned from the genotypes, which is what
# ComplexHeatmap's row_split does.
.cluster_within <- function(G, blocks, cluster = TRUE, nominal_ids = character(0)) {
  ord <- lapply(levels(blocks), function(b) {
    ids <- rownames(G)[blocks == b]
    if (!cluster || length(ids) < 3) return(list(ids = ids, hc = NULL))
    sub <- G[ids, , drop = FALSE]
    # a block whose samples are identical (or all-missing) has nothing to cluster on
    d <- try(.geno_dist(sub, nominal_ids), silent = TRUE)
    if (inherits(d, "try-error") || !any(is.finite(d)) || max(d, na.rm = TRUE) == 0)
      return(list(ids = ids, hc = NULL))
    d[!is.finite(d)] <- max(d, na.rm = TRUE)
    hc <- stats::hclust(d, method = "ward.D2")
    list(ids = ids[hc$order], hc = hc)
  })
  names(ord) <- levels(blocks)
  ord
}

# Dendrogram segments for every block, in the heatmap's row coordinates: leaves sit on the
# row they label, so the two panels line up facet by facet.
.block_dendro <- function(ord, rows) {
  segs <- lapply(names(ord), function(b) {
    hc <- ord[[b]]$hc
    if (is.null(hc)) return(NULL)
    s <- .dendro_segments(hc)
    # .dendro_segments() puts leaf k of the *clustered order* at x = k; the heatmap numbers
    # its rows the same way inside the block, so the leaf index IS the row index
    base <- min(rows$.row[rows$.split == b]) - 1
    cbind(data.frame(y = s$x + base, yend = s$xend + base, x = s$y, xend = s$yend),
          .split_cols(rows, rep(b, length(s$x))))
  })
  segs <- do.call(rbind, segs[!vapply(segs, is.null, logical(1))])
  if (is.null(segs) || !nrow(segs)) return(NULL)
  rownames(segs) <- NULL
  segs
}

# ---- one block per level, or per combination of levels ---------------------------------
# The rows carry `.split`, one factor whose levels are the blocks in drawing order, and --
# when `split` named more than one column -- `.split1`, `.split2`, ... holding each column's
# own factor. The panels facet on those, so ggplot draws one strip per column ("Northern" |
# "REF") instead of a single pasted label, while `.split` stays the one key every other
# piece of bookkeeping (row numbering, dendrogram anchors) is written against.
.split_vars <- function(rows) {
  nested <- grep("^\\.split[0-9]+$", names(rows), value = TRUE)
  if (length(nested)) nested else ".split"
}

# The split columns for a set of blocks, by the block's `.split` label. Every panel that
# facets alongside the heatmap builds its data with this, so they all share the same factors.
.split_cols <- function(rows, blocks) {
  vars <- unique(c(".split", .split_vars(rows)))
  out <- rows[match(as.character(blocks), as.character(rows$.split)), vars, drop = FALSE]
  rownames(out) <- NULL
  out
}

.split_facet <- function(rows) {
  facets <- lapply(.split_vars(rows), function(v) call("[[", quote(.data), v))
  ggplot2::facet_grid(rows = do.call(ggplot2::vars, facets), scales = "free_y",
                      space = "free_y")
}

# Which SNPs to mark: `chr:pos` ids, bare positions, a gene name resolved in `genes`, or an
# interval table (every SNP inside each interval -- an aa_intervals() codon table, say, which
# saves the caller doing coordinate arithmetic on it).
.resolve_marks <- function(marks, loci, genes, reference) {
  if (is.null(marks) || !length(marks)) return(numeric(0))
  if (is.data.frame(marks)) {
    if (!all(c("start", "end") %in% names(marks)))
      stop("a data frame `mark_snps` needs start and end columns", call. = FALSE)
    ch <- if ("chr" %in% names(marks)) normalise_chr(marks$chr)
          else if ("chrom" %in% names(marks)) normalise_chr(marks$chrom) else NA_character_
    out <- unlist(lapply(seq_len(nrow(marks)), function(i) {
      keep <- loci$pos >= marks$start[i] & loci$pos < marks$end[i]
      if (!is.na(ch[i])) keep <- keep & loci$chr == ch[i]
      loci$pos[keep]
    }), use.names = FALSE)
    if (!length(out))
      message("no genotyped SNP inside any `mark_snps` interval, so nothing is marked")
    return(out)
  }
  if (is.numeric(marks)) return(marks)
  out <- lapply(marks, function(m) {
    hit <- which(loci$snp_id == m | loci$id_norm == m)
    if (length(hit)) return(loci$pos[hit])
    iv <- try(.resolve_region(m, genes, reference), silent = TRUE)
    if (inherits(iv, "try-error"))
      stop("`mark_snps = \"", m, "\"` is not a SNP in the window, a position, or a gene ",
           "name in `genes`", call. = FALSE)
    hit <- loci$pos[loci$pos >= iv$start & loci$pos < iv$end]
    # a gene can resolve fine and still hold no genotyped SNP, which would otherwise just
    # draw nothing and look like the argument was ignored
    if (!length(hit)) message("no genotyped SNP inside ", m, ", so nothing is marked there")
    hit
  })
  unlist(out, use.names = FALSE)
}

#' Genotypes over one region, clustered by sample
#'
#' A genotype heatmap for a single interval: one row per sample, one column per SNP, with the
#' samples clustered and a dendrogram beside them, and a gene track underneath. Splitting the
#' rows by a metadata column fixes the blocks and clusters *within* each one, so a haplotype
#' shared across a group shows up as a solid band rather than being scattered by a
#' genome-wide ordering.
#'
#' The genotypes should be the **full** panel, not an LD-pruned one. Pruning keeps one SNP out
#' of each correlated run, and a correlated run is what a shared haplotype is -- so a pruned
#' panel shows fewer SNPs, chosen to be as uncorrelated as possible, and the blocks come out
#' thinner than they are. Build the object with `load_genotypes(prune = FALSE)` for this plot
#' and keep the pruned one for PCA / UMAP / admixture, where pruning is what you want; the
#' plot says so when the object records that it was pruned.
#'
#' `spacing` decides what the horizontal axis means, and the three answers show different
#' things. `"even"` gives every SNP the same width, which is how the haplotype structure is
#' easiest to read but says nothing about distance. `"genomic"` puts each SNP at its real
#' coordinate, so a dense cluster of SNPs looks dense -- correct about position, but sparse
#' stretches become wide empty bands, and SNPs closer together than one mark width merge into
#' a single block that cannot be told from one wide SNP.
#'
#' `"gapped"` is the middle ground: every SNP keeps a full readable column as under `"even"`,
#' and an empty stretch of genome buys blank columns -- one per `gap_unit`, up to `gap_max`
#' for any single gap. A desert then reads as a visible gap instead of vanishing, without a
#' long one taking over the panel. `gap_unit` defaults to a fiftieth of the window, so
#' ordinary spacing between SNPs costs nothing and only a stretch noticeably emptier than the
#' rest opens up. The cap is what keeps it a compression rather than a coordinate: distances
#' come out ordered and roughly proportional, not to scale, so read it for "there is a lot of
#' nothing here", not for how much.
#'
#' It also gives the gene track somewhere to draw. Under `"even"` a gene with no genotyped SNP
#' has no columns and is dropped from the track with a message; under `"gapped"` the blank
#' columns are real space, so a gene sitting in a desert still appears, in about the right
#' place.
#'
#' Either way a SNP is only ever drawn over the genes it really falls in. Under `"even"` the
#' axis counts SNPs, so a gene's box is exactly the columns it holds: its width says how many
#' SNPs are in it, **not** how long it is, and a gene with no genotyped SNP in the window has no
#' width at all and is left off the track with a message. (Interpolating genomic bounds onto a
#' SNP-index axis instead puts each gene edge at a fractional column, and since a tile occupies a
#' whole column, a SNP just outside a gene ends up drawn over it.) Under `"genomic"` the boxes
#' are true extents and it is the fixed-width SNP marks that are clipped, so a mark near a gene
#' edge stops there rather than reaching into its neighbour.
#'
#' @param x A [PopStructure] object (its genotypes supply the calls, its metadata the split).
#' @param region The interval to draw: a gene name from `genes`, a range
#'   (`"13:1,720,000-1,730,000"`), a whole chromosome, or a one-row data frame with
#'   chr/start/end.
#' @param split Optional metadata column(s) whose levels block the rows. Samples are
#'   clustered inside each block, and the blocks keep the column's level order -- a factor's
#'   levels are honoured, so a `region` ordered geographically stays geographic. More than
#'   one column nests the blocks in the order given: `split = c("region", "PIN_variant")`
#'   makes one block per region, each divided by variant, with a strip per column and only
#'   the combinations that hold a sample drawn. Samples missing any of the columns are
#'   dropped with a message.
#' @param genotypes Optional alternative calls to draw: a genotype matrix (samples x SNPs
#'   with `chr:pos` column names), a [load_genotypes()] list, or another [PopStructure]. `NULL`
#'   (default) uses `x`'s own matrix. Metadata, grouping and the active sample set always come
#'   from `x`, so one pruned object can supply the annotations while the full panel supplies
#'   the calls -- which is the combination this plot wants.
#' @param annotations Optional metadata columns to draw as coloured strips down the right,
#'   one column each, sharing the object's colour maps (see [meta_colors()]) so a level keeps
#'   the colour it has in the other plots. Each gets its own legend.
#' @param annotation_colours,annotation_colors Per-annotation palettes, as a named list
#'   (`list(region = c(north = "#1B9E77", south = "grey60"))`). Named colours override those
#'   levels and leave the rest of that annotation's shared map alone, so recolouring one
#'   level does not mean respelling all of them; an unnamed vector is taken in level order.
#'   Only this plot changes -- `PopStructure$set_colors()` is still the way to move a colour
#'   everywhere at once, which is what keeps a level looking the same across figures.
#' @param border Outline every call (default `TRUE`), which is what makes single SNPs
#'   readable as cells rather than a wash of colour.
#' @param border_colour,border_color Colour of that outline.
#' @param allele Which allele the dosages count, `"alt"` or `"ref"`. `NULL` (default) asks
#'   the object; if it does not record it, alt is assumed and a message says so. The two are
#'   indistinguishable from the matrix, and getting it backwards mislabels every call.
#' @param samples Optional sample ids to keep.
#' @param spacing `"even"` (default) gives every SNP equal width; `"genomic"` places each at
#'   its real coordinate; `"gapped"` keeps the equal widths and inserts blank columns for
#'   empty stretches of genome, so distance is visible without the SNPs shrinking.
#' @param gap_unit,gap_max Under `"gapped"` spacing, base pairs of empty genome per blank
#'   column and the most blank columns any one gap may claim (default `10`). `gap_unit`
#'   defaults to `NULL`, a fiftieth of the window; lower it to exaggerate distance, raise it
#'   to play it down. `gap_max` stops a single desert from squeezing every SNP into the
#'   margin.
#' @param cluster Cluster the samples (default `TRUE`). `FALSE` keeps them in the order they
#'   arrive, which is worth doing when the metadata order is the point.
#' @param dendrogram Draw the dendrogram beside the rows (needs `cluster`).
#' @param dend_width Width of the dendrogram panel, as a fraction of the heatmap's.
#' @param mark_snps Optional SNPs to mark with a vertical line: `chr:pos` ids, bare positions,
#'   a gene name from `genes`, or an interval table with `start`/`end` (and `chr`) -- every SNP
#'   inside each interval is marked, so an [aa_intervals()] codon table can be passed straight
#'   in without converting its coordinates.
#' @param mark_colour Colour for those lines.
#' @param genes Gene table for the track and for resolving names (e.g. [PF3D7_GENES]).
#' @param gene_track Draw the gene track under the heatmap (default `TRUE` when `genes` is
#'   given). Under `"even"` spacing the boxes are mapped onto the SNP columns they cover, so
#'   the track still says which columns sit in which gene.
#' @param gene_label_angle Rotation for the gene names, in degrees; `45` or `90` keeps long
#'   systematic ids from colliding.
#' @param pad,min_span Context around `region`, as in [plot_ibd_locus()]: one value pads both
#'   sides, two the left and the right (named `left` / `right` if you like).
#' @param max_snps Refuse to draw more than this many columns (default 2000). A window with
#'   thousands of SNPs is unreadable as tiles; narrow it rather than have it silently thinned.
#' @param snp_width Width of each mark under `"genomic"` spacing, in base pairs. `NULL`
#'   (default) uses 0.5% of the window, wide enough to see and narrow enough to leave the
#'   gaps between SNPs visible.
#' @param additional_genotypes Optional VCF/BCF path(s) holding markers the genotype matrix
#'   cannot express. A dosage counts copies of one allele, so a site with three alleles
#'   collapses every non-reference call to the same number -- which is why such sites are
#'   usually dropped from a callset in the first place. Markers given here are read as the
#'   **set of alleles** each sample carries, so `reference`, `alternate 1`, `alternate 2` and
#'   the mixed states between them stay distinct. Needs `bcftools` on `PATH`; samples are
#'   matched by name, and a position already in the genotypes is replaced with the richer
#'   allele-set form read here (a dosage column cannot keep two alternates apart). Only the states
#'   that actually occur are added to the legend: the full enumeration of a triallelic site
#'   is seven, and most of them are ordinarily empty.
#' @param colours,colors Named fill colours for `reference` / `mixed` / `alternate`, and for
#'   any state an `additional_genotypes` marker contributes.
#' @param na_colour,na_color Fill for missing calls.
#' @param show_sample_names Label the rows. `NULL` (default) labels them when there are at
#'   most 40 samples.
#' @param reference Reference id, used when `region` names a whole chromosome.
#' @return A patchwork of the dendrogram, heatmap and gene track, or a plain ggplot when
#'   neither of those two is drawn.
#' @examples
#' ps <- example_pop_structure(umap = FALSE)
#' plot_region_haplotypes(ps, "pfcrt", split = "country", pad = 20000,
#'                        genes = PF_EXAMPLE_DRUG_GENES)
#' @export
plot_region_haplotypes <- function(x, region, split = NULL, annotations = NULL,
                                   genotypes = NULL, samples = NULL,
                                   spacing = c("even", "genomic", "gapped"),
                                   cluster = TRUE, dendrogram = TRUE, dend_width = 0.15,
                                   border = TRUE, border_colour = "grey45",
                                   allele = NULL,
                                   mark_snps = NULL, mark_colour = "#B2182B",
                                   genes = NULL, gene_track = NULL, gene_label_angle = 0,
                                   pad = 0, min_span = 0, max_snps = 2000, snp_width = NULL,
                                   gap_unit = NULL, gap_max = 10,
                                   colours = NULL, na_colour = "grey85",
                                   show_sample_names = NULL,
                                   reference = DEFAULT_REFERENCE,
                                   annotation_colours = NULL,
                                   additional_genotypes = NULL,
                                   border_color = NULL, colors = NULL, na_color = NULL,
                                   annotation_colors = NULL) {
  annotation_colours <- .alias_arg("annotation_colours", "annotation_colors")
  border_colour <- .alias_arg("border_colour", "border_color")
  colours <- .alias_arg("colours", "colors")
  na_colour <- .alias_arg("na_colour", "na_color")
  .need_package("ggplot2", "plot_region_haplotypes()")
  spacing <- match.arg(spacing)
  if (is.null(gene_track)) gene_track <- !is.null(genes)

  gt <- .haplotype_genotypes(x, genotypes)
  G <- gt$G
  # Markers the dosage matrix cannot express -- a site with three alleles collapses every
  # non-reference call to one number -- read from their own callset as allele sets and
  # spliced in here, before the window is chosen, so they are laid out like any other column.
  state_levels <- list()
  for (v in additional_genotypes) {
    add <- .read_genotype_sets(v)
    have <- unique(normalise_chr(sub(":.*", "", colnames(G))))
    lookup <- stats::setNames(unique(sub(":.*", "", colnames(G))), have)
    as_have <- lookup[normalise_chr(sub(":.*", "", colnames(add$codes)))]
    ids <- ifelse(is.na(as_have), colnames(add$codes),
                  paste0(as_have, ":", sub(".*:", "", colnames(add$codes))))
    colnames(add$codes) <- ids
    names(add$levels) <- ids
    miss <- setdiff(rownames(G), rownames(add$codes))
    if (length(miss))
      stop(sprintf("%d of the genotypes' samples are not in %s (e.g. %s)", length(miss),
                   basename(v), paste(utils::head(miss, 3), collapse = ", ")), call. = FALSE)
    clash <- intersect(colnames(G), ids)
    if (length(clash)) {
      # The panel may already carry the biallelic position of a codon (a combined callset
      # that keeps the codon sites inline does). The allele-set form read here is the richer
      # one -- it keeps alternate 1 / alternate 2 distinct, which a dosage column cannot --
      # so replace the existing column rather than refusing.
      message(sprintf("%s: %d position(s) already in the panel replaced with their allele-set form (%s)",
                      basename(v), length(clash), paste(utils::head(clash, 5), collapse = ", ")))
      G <- G[, setdiff(colnames(G), clash), drop = FALSE]
    }
    # A left join on the genotypes being plotted: the marker is usually called on the whole
    # cohort while the figure shows a subset, so extras are expected -- but say how many, or
    # a name mismatch that silently drops half the callset looks like a clean merge.
    spare <- setdiff(rownames(add$codes), rownames(G))
    if (length(spare))
      message(sprintf("%s has %d sample(s) not in the genotypes; they are left out",
                      basename(v), length(spare)))
    G <- cbind(G, add$codes[rownames(G), , drop = FALSE])
    state_levels <- c(state_levels, add$levels)
  }
  if (!is.null(samples)) {
    keep <- rownames(G) %in% samples
    if (!any(keep)) stop("none of `samples` are in the genotypes", call. = FALSE)
    G <- G[keep, , drop = FALSE]
  }

  # ---- the window and the SNPs in it ---------------------------------------
  loci <- .parse_snp_ids(colnames(G))
  loci$snp_id <- colnames(G)
  loci$chr <- normalise_chr(loci$chr)
  loci$id_norm <- paste0(loci$chr, ":", loci$pos)
  iv <- .pad_region(.resolve_region(region, genes, reference), pad,
                    .chrom_layout(reference), min_span)
  sel <- loci[loci$chr == iv$chr & loci$pos >= iv$start & loci$pos <= iv$end, , drop = FALSE]
  if (!nrow(sel))
    stop("no genotyped SNPs in ", iv$chr, ":", format(iv$start, scientific = FALSE), "-",
         format(iv$end, scientific = FALSE), "; widen the window with `pad`", call. = FALSE)
  if (nrow(sel) > max_snps)
    stop(nrow(sel), " SNPs in this window is more than `max_snps` (", max_snps, "); narrow ",
         "the window, or raise the limit if you really want that many columns",
         call. = FALSE)
  sel <- sel[order(sel$pos), , drop = FALSE]
  G <- G[, sel$snp_id, drop = FALSE]

  # ---- rows: blocked by the metadata, clustered inside each block ----------
  blocks <- .row_blocks(x, split, rownames(G))
  if (!is.null(blocks$dropped))
    message("dropped ", blocks$dropped, " sample(s) with no ",
            paste(split, collapse = " / "))
  G <- G[names(blocks$f), , drop = FALSE]
  # the allele-set columns are nominal, so the distance must not scale them
  ord <- .cluster_within(G, blocks$f, cluster, nominal_ids = base::names(state_levels))
  row_ids <- unlist(lapply(ord, `[[`, "ids"), use.names = FALSE)
  rows <- data.frame(
    sample = row_ids,
    .split = factor(as.character(blocks$f[row_ids]), levels = levels(blocks$f)),
    stringsAsFactors = FALSE)
  # nested splits keep each column's own factor beside the combined key (see .split_vars)
  if (!is.null(blocks$parts))
    rows <- cbind(rows, blocks$parts[match(row_ids, rownames(blocks$parts)), , drop = FALSE])
  # rows are numbered top to bottom within the whole plot; each facet then shows its own slice
  rows$.row <- seq_len(nrow(rows))

  # ---- long form, with the x axis the chosen spacing means -----------------
  long <- data.frame(
    sample = rep(row_ids, times = nrow(sel)),
    snp_id = rep(sel$snp_id, each = length(row_ids)),
    pos = rep(sel$pos, each = length(row_ids)),
    col = rep(seq_len(nrow(sel)), each = length(row_ids)),
    value = as.vector(G[row_ids, , drop = FALSE]),
    stringsAsFactors = FALSE)
  long$.row <- rows$.row[match(long$sample, rows$sample)]
  long <- cbind(long, .split_cols(rows, rows$.split[match(long$sample, rows$sample)]))
  long$call <- .geno_calls(long$value, allele %||% gt$allele %||% .resolve_allele(x, NULL),
                           long$snp_id, state_levels)
  # LD pruning keeps one SNP out of each correlated run, which is exactly what a shared
  # haplotype is made of, so a pruned panel understates the very structure this plot is for.
  if (isTRUE(gt$pruned))
    message("these genotypes are LD-pruned, so this window shows only the SNPs that survived ",
            "pruning and the haplotype blocks will look thinner than they are; rebuild with ",
            "`load_genotypes(..., prune = FALSE)` for a haplotype view")

  tiles <- .snp_tile_x(sel, spacing, snp_width, gap_unit = gap_unit, gap_max = gap_max)
  gene_boxes <- .gene_boxes_in_x(genes, iv, sel, tiles, spacing)
  # a mark must never straddle into a gene the SNP is not in (see .clip_tiles_to_genes)
  if (spacing == "genomic") tiles <- .clip_tiles_to_genes(tiles, sel, gene_boxes)
  long$xmin <- tiles$xmin[long$col]
  long$xmax <- tiles$xmax[long$col]
  # under genomic spacing the marks no longer reach the window edges, but the axis should
  xlim <- if (spacing == "genomic") c(min(iv$start, min(tiles$xmin)),
                                      max(iv$end, max(tiles$xmax)))
          else c(min(tiles$xmin), max(tiles$xmax))
  # a gene box interpolated past the outermost SNP would be clipped away; let the axis hold it
  if (spacing == "gapped" && !is.null(gene_boxes) && nrow(gene_boxes))
    xlim <- c(min(xlim[1], gene_boxes$.gene_xmin), max(xlim[2], gene_boxes$.gene_xmax))

  # ---- the heatmap ---------------------------------------------------------
  fills <- .GENO_FILL
  extra <- intersect(setdiff(levels(long$call), names(fills)), as.character(long$call))
  fills <- .assign_extra_fills(fills, extra)
  if (!is.null(colours)) fills[names(colours)] <- unname(colours)
  labels <- if (is.null(show_sample_names)) nrow(rows) <= 40 else isTRUE(show_sample_names)

  p <- ggplot2::ggplot(long) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                   ymin = .data$.row - 0.5, ymax = .data$.row + 0.5, fill = .data$call),
      colour = if (border) border_colour else NA,
      linewidth = if (border) 0.06 else 0) +
    ggplot2::scale_fill_manual(values = fills, na.value = na_colour,
                               # the three base calls always, plus only those extra states a
                               # marker actually produced -- a triallelic site has seven and
                               # is normally missing most of them
                               limits = unique(c(.GENO_LEVELS,
                                                 intersect(levels(long$call),
                                                           as.character(long$call)))),
                               drop = FALSE, name = "call",
                               guide = ggplot2::guide_legend(order = .HAP_LEGEND_CALL)) +
    ggplot2::scale_y_reverse(
      breaks = if (labels) rows$.row, labels = if (labels) rows$sample,
      expand = ggplot2::expansion(0)) +
    .haplotype_x_scale(sel, spacing, iv, xlim, x_of = tiles$x) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   panel.spacing.y = grid::unit(2, "pt"),
                   axis.ticks.y = if (labels) ggplot2::element_line() else
                     ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(),
                   strip.background = ggplot2::element_rect(fill = "grey95", colour = NA),
                   strip.text.y.right = ggplot2::element_text(angle = 0))

  marks <- .resolve_marks(mark_snps, sel, genes, reference)
  if (length(marks)) {
    mx <- .marks_to_x(marks, sel, tiles, spacing)
    if (length(mx))
      p <- p + ggplot2::geom_vline(xintercept = mx, colour = mark_colour, linewidth = 0.35)
  }
  faceted <- nlevels(rows$.split) > 1
  if (faceted) p <- p + .split_facet(rows)

  # ---- dendrogram beside the rows, gene track underneath -------------------
  dend <- if (dendrogram && cluster) .dendro_panel(.block_dendro(ord, rows), rows,
                                                   faceted) else NULL
  ann <- .hap_annotation_panel(x, annotations, rows, faceted, border, border_colour,
                               annotation_colours)
  # only the rightmost panel names the blocks, or every block is labelled twice
  if (!is.null(ann) && faceted)
    p <- p + ggplot2::theme(strip.text.y = ggplot2::element_blank(),
                            strip.text.y.right = ggplot2::element_blank(),
                            strip.background = ggplot2::element_blank())
  track <- if (gene_track) {
    boxes <- gene_boxes
    .gene_track_panel(boxes, xlim, angle = gene_label_angle, width_in = .ZOOM_WIDTH_IN,
                      min_width = spacing == "genomic")
  } else NULL

  out <- .assemble_haplotype_panels(p, dend, ann, track, dend_width)
  # A labelled row has to be tall enough for its text; an unlabelled one only has to be
  # visible, so hundreds of samples compress into a readable page instead of a 30-inch one.
  row_in <- if (labels) 0.16 else max(0.015, min(0.09, .HAP_ROWS_IN / nrow(rows)))
  attr(out, "plasgenomics_dims") <- c(
    width = .ZOOM_WIDTH_IN + (if (labels) 1 else 0) +
      # the strips sit beside the annotation panel, one per split column
      (if (is.null(ann)) 0 else attr(ann, "ann_in") + 0.6 * length(.split_vars(rows))),
    height = max(3, row_in * nrow(rows) + 0.3 * max(1L, nlevels(rows$.split)) + 1) +
      if (is.null(track)) 0 else attr(track, "track_in"))
  out
}

# Coloured strips down the right, one column per annotation. Each needs its own fill scale --
# a level of "region" and a level of "year" are unrelated and must not share a palette -- which
# is what ggnewscale is for: without it a single panel can only carry one fill mapping.
# Colours come from the object's own maps, so a level keeps the colour it has in the UMAP or
# the admixture bars.
.hap_annotation_panel <- function(x, cols, rows, faceted, border, border_colour,
                                  ann_colours = NULL) {
  if (is.null(cols) || !length(cols)) return(NULL)
  .need_package("ggnewscale", "annotation strips in plot_region_haplotypes()")
  meta <- x$get_meta()
  if (is.null(meta)) stop("`annotations` needs metadata on the object", call. = FALSE)
  missing <- setdiff(cols, names(meta))
  if (length(missing))
    stop("not a metadata column: ", paste(missing, collapse = ", "), call. = FALSE)

  key <- if ("sample" %in% names(meta)) meta$sample else rownames(meta)
  maps <- x$get_colors()
  if (is.null(maps)) maps <- list()
  need <- setdiff(cols, names(maps))
  if (length(need)) maps[need] <- meta_colors(meta, cols = need)[need]

  p <- ggplot2::ggplot() +
    ggplot2::scale_x_continuous(breaks = seq_along(cols), labels = cols,
                                expand = ggplot2::expansion(0), position = "top") +
    ggplot2::scale_y_reverse(expand = ggplot2::expansion(0))
  if (!is.null(ann_colours)) {
    if (!is.list(ann_colours))
      stop("`annotation_colours` is a named list, one entry per annotation, e.g. ",
           "list(", cols[1], " = c(level = \"#1B9E77\"))", call. = FALSE)
    stray <- setdiff(names(ann_colours), cols)
    if (length(stray))
      warning("`annotation_colours` names something that is not an annotation: ",
              paste(stray, collapse = ", "), ". Annotations here: ",
              paste(cols, collapse = ", "), call. = FALSE)
  }

  for (k in seq_along(cols)) {
    cc <- cols[k]
    # A factor annotation keeps every level it was built with, so a level whose samples
    # are all outside this plot would otherwise be drawn as a legend key with nothing
    # behind it -- and blank, wherever the shared colour map has no entry for it.
    d <- cbind(
      data.frame(.row = rows$.row, x = k,
                 value = droplevels(.as_group_factor(meta[[cc]][match(rows$sample, key)]))),
      .split_cols(rows, rows$.split))
    lev <- levels(d$value)
    maps[[cc]] <- .fill_palette_gaps(maps[[cc]], lev, meta, cc)
    maps[[cc]] <- .merge_palette(maps[[cc]], ann_colours[[cc]], lev, cc)
    p <- p +
      ggplot2::geom_rect(
        data = d,
        ggplot2::aes(xmin = .data$x - 0.5, xmax = .data$x + 0.5,
                     ymin = .data$.row - 0.5, ymax = .data$.row + 0.5, fill = .data$value),
        colour = if (border) border_colour else NA,
        linewidth = if (border) 0.06 else 0) +
      # numbered by the position in `annotations`, so the legends read in the order asked for
      ggplot2::scale_fill_manual(values = maps[[cc]], na.value = "grey85", drop = FALSE,
                                 name = cc,
                                 guide = ggplot2::guide_legend(order = .HAP_LEGEND_CALL + k)) +
      ggnewscale::new_scale_fill()
  }
  p <- p +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(), panel.spacing.y = grid::unit(2, "pt"),
      axis.title = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text.x.top = ggplot2::element_text(angle = 90, hjust = 0, vjust = 0.5, size = 7),
      plot.margin = ggplot2::margin(l = 2, r = 2, t = 0, b = 0),
      strip.background = ggplot2::element_rect(fill = "grey95", colour = NA),
      strip.text.y.right = ggplot2::element_text(angle = 0))
  if (faceted) p <- p + .split_facet(rows)
  attr(p, "ann_in") <- .HAP_ANN_IN * length(cols)
  p
}

# The calls to draw. By default the object's own matrix, but a haplotype view usually wants
# the unpruned panel while PCA / admixture want the pruned one, and holding both is cheaper
# than rebuilding: `genotypes` takes a matrix, a load_genotypes() list, or another PopStructure,
# and metadata still comes from `x`. Whatever it is, it also carries whichever of the allele
# coding and the pruning flag it knows, since those are what the plot cannot infer.
.haplotype_genotypes <- function(x, genotypes) {
  src <- list(allele = NULL, pruned = NULL)
  if (is.null(genotypes)) {
    # the full panel when the object has one: pruning drops the correlated SNPs that make a
    # shared haplotype a solid band, which is the whole point of this plot
    panel <- if ("full" %in% x$panels()) "full" else NULL
    G <- x$genotype(prefer = "full")
    src$allele <- x$allele(panel)
    src$pruned <- x$pruned(panel)
  } else if (inherits(genotypes, "PopStructure")) {
    G <- genotypes$genotype()
    src$allele <- genotypes$allele(); src$pruned <- genotypes$pruned()
  } else if (is.list(genotypes) && !is.null(genotypes$genotype)) {
    G <- genotypes$genotype
    if (is.null(rownames(G)) && !is.null(genotypes$sample.id))
      rownames(G) <- genotypes$sample.id
    src$allele <- genotypes$allele; src$pruned <- genotypes$pruned
  } else {
    G <- as.matrix(genotypes)
  }
  if (is.null(G) || !nrow(G) || !ncol(G))
    stop("no genotypes to draw", call. = FALSE)
  if (is.null(rownames(G)) || is.null(colnames(G)))
    stop("`genotypes` needs sample row names and `chr:pos` column names", call. = FALSE)

  # keep the object's active samples, so restrict() / subset() still decide who is drawn
  if (!is.null(genotypes)) {
    want <- x$get_samples()
    keep <- intersect(want, rownames(G))
    if (!length(keep))
      stop("`genotypes` has none of this object's samples", call. = FALSE)
    if (length(keep) < length(want))
      message("`genotypes` is missing ", length(want) - length(keep), " of the object's ",
              length(want), " samples; drawing the ", length(keep), " it has")
    G <- G[keep, , drop = FALSE]
  }
  list(G = G, allele = src$allele, pruned = src$pruned)
}

# The metadata column(s) that block the rows, as a factor over the samples being drawn.
# Several columns nest: the first is the outer block and each later one divides it, with
# every column's own level order kept (so an ordered factor stays ordered) and only the
# combinations that actually hold a sample surviving. `parts` then carries each column's
# factor, one row per kept sample, so the panels can facet on all of them.
.row_blocks <- function(x, split, ids) {
  if (is.null(split) || !length(split)) {
    f <- factor(rep("all", length(ids)), levels = "all")
    names(f) <- ids
    return(list(f = f, dropped = NULL, parts = NULL))
  }
  meta <- x$get_meta()
  missing <- if (is.null(meta)) split else setdiff(split, names(meta))
  if (length(missing))
    stop("`split = \"", paste(missing, collapse = "\", \""), "\"` is not a metadata column",
         call. = FALSE)
  key <- if ("sample" %in% names(meta)) meta$sample else rownames(meta)
  parts <- lapply(split, function(cc) .as_group_factor(meta[[cc]][match(ids, key)]))
  keep <- Reduce(`&`, lapply(parts, function(v) !is.na(v)))
  parts <- lapply(parts, function(v) droplevels(v[keep]))
  f <- if (length(parts) == 1) parts[[1]] else
    droplevels(interaction(parts, sep = " / ", lex.order = TRUE))
  names(f) <- ids[keep]
  out <- list(f = f, dropped = if (all(keep)) NULL else sum(!keep), parts = NULL)
  if (length(parts) > 1) {
    out$parts <- as.data.frame(stats::setNames(parts, paste0(".split", seq_along(parts))),
                               stringsAsFactors = FALSE)
    rownames(out$parts) <- ids[keep]
  }
  out
}

# Tile edges per SNP. Even spacing gives each column the same width. So does genomic spacing,
# but centred on the SNP's real coordinate: equal marks at true positions are what let you see
# how far apart the SNPs are. Stretching each tile to meet its neighbours would fill the gaps
# back in and hide exactly that -- and would make an isolated SNP a wide block purely because
# nothing was called near it.
.snp_tile_x <- function(sel, spacing, snp_width = NULL, gap_unit = NULL, gap_max = 10) {
  n <- nrow(sel)
  if (spacing == "even" || n == 1)
    return(list(xmin = seq_len(n) - 0.5, xmax = seq_len(n) + 0.5))
  if (spacing == "gapped") {
    x <- .gapped_x(sel$pos, gap_unit, gap_max)
    return(list(xmin = x - 0.5, xmax = x + 0.5, x = x))
  }
  pos <- sel$pos
  span <- max(diff(range(pos)), 1)
  w <- if (!is.null(snp_width)) snp_width else max(span * 0.005, 1)
  list(xmin = pos - w / 2, xmax = pos + w / 2)
}

# Column index for each SNP under "gapped" spacing: one column per SNP, plus blank columns
# where the genome between two SNPs is empty.
#
# Even spacing draws a 40 kb desert and a 40 bp gap identically; genomic spacing draws the
# desert honestly and the SNPs in it too small to read. This keeps every SNP one readable
# column wide and spends `1` blank column per `gap_unit` of empty sequence, so a gap is
# visible and roughly proportional without a long one taking over the panel -- `gap_max`
# caps what any single gap can claim.
.gapped_x <- function(pos, gap_unit = NULL, gap_max = 10) {
  n <- length(pos)
  if (n <= 1) return(seq_len(n))
  gaps <- diff(pos)
  # A fixed number of base pairs per blank column cannot suit both a 40 kb gene window and a
  # whole chromosome, so the default is a fiftieth of the window: ordinary spacing between
  # SNPs stays under it and costs nothing, while a stretch of genome noticeably emptier than
  # the rest buys columns in proportion. Give a number to fix it instead.
  if (is.null(gap_unit)) gap_unit <- max(diff(range(pos)) / 50, 1)
  blanks <- pmin(floor(gaps / gap_unit), gap_max)
  cumsum(c(1, 1 + blanks))
}

# Genomic position -> plot x, for anything that is not a SNP (a gene edge, a marked site).
# Piecewise linear between the SNP anchors, so a coordinate between two SNPs lands between
# their columns, and the blank columns of a gap are real space something can be drawn in.
.gapped_pos_to_x <- function(v, pos, x) {
  if (!length(v)) return(numeric(0))
  if (length(pos) == 1) return(rep(x[1], length(v)))
  stats::approx(pos, x, xout = v, rule = 2, ties = "ordered")$y
}

# Under even spacing the axis counts SNPs, so it is labelled with the coordinates of a few of
# them rather than pretending the numbers are positions.
.haplotype_x_scale <- function(sel, spacing, iv, xlim, x_of = NULL) {
  if (spacing == "genomic") {
    unit <- if (diff(xlim) >= 1e4) 1e3 else 1
    return(ggplot2::scale_x_continuous(
      name = paste0("chromosome ", iv$chr, " position (", if (unit == 1) "bp" else "kb", ")"),
      labels = function(v) format(round(v / unit, if (unit == 1) 0 else 1), big.mark = ",",
                                  trim = TRUE, scientific = FALSE),
      expand = ggplot2::expansion(0)))
  }
  at <- unique(round(seq(1, nrow(sel), length.out = min(nrow(sel), 6))))
  # the axis counts columns either way, so it is labelled with the coordinates of a few SNPs
  # rather than pretending the numbers are positions
  breaks <- if (spacing == "gapped") x_of[at] else at
  name <- if (spacing == "gapped")
    paste0("chromosome ", iv$chr, ": ", nrow(sel), " SNPs, gaps compressed")
  else paste0("chromosome ", iv$chr, ": ", nrow(sel), " SNPs, evenly spaced")
  ggplot2::scale_x_continuous(
    name = name, breaks = breaks,
    labels = format(round(sel$pos[at] / 1000, 1), big.mark = ",", trim = TRUE),
    expand = ggplot2::expansion(0))
}

# Positions -> the plot's x units (a column index under even spacing).
.marks_to_x <- function(marks, sel, tiles, spacing) {
  if (spacing == "genomic") return(marks)
  # gapped spacing has blank columns between the SNPs, so a marked position with no SNP of
  # its own still has somewhere to sit: interpolate it rather than dropping it
  if (spacing == "gapped") return(.gapped_pos_to_x(marks, sel$pos, tiles$x))
  i <- match(marks, sel$pos)
  # a marked position with no genotyped SNP has no column to sit on
  (seq_len(nrow(sel)))[i[!is.na(i)]]
}

# Gene boxes in the plot's x units. Under even spacing a gene spans the columns it covers,
# interpolated so a gene between two SNPs still lands between their columns.
.gene_boxes_in_x <- function(genes, iv, sel, tiles, spacing) {
  if (is.null(genes)) return(NULL)
  g <- .gene_track(genes)
  g <- g[g$chr == iv$chr & as.numeric(g$end) >= iv$start & as.numeric(g$start) <= iv$end, ,
         drop = FALSE]
  if (!nrow(g)) return(NULL)
  if (spacing == "genomic") {
    g$.gene_xmin <- as.numeric(g$start)
    g$.gene_xmax <- as.numeric(g$end)
    return(g)
  }
  if (spacing == "gapped") return(.gapped_gene_boxes(g, sel, tiles$x))

  # Under even spacing the axis counts SNPs, so a gene's extent on it is exactly the columns it
  # holds -- nothing else is well defined. Interpolating its genomic bounds onto the axis instead
  # puts each edge at a fractional column, and since a tile occupies a whole column, a SNP just
  # outside a gene ends up drawn over that gene's box. Snap to the columns.
  #
  # A gene with no SNP in the window therefore has no extent at all, and any box drawn for it
  # would necessarily sit under some other gene's SNPs; those are dropped and reported.
  in_gene <- lapply(seq_len(nrow(g)), function(j)
    which(sel$pos >= as.numeric(g$start[j]) & sel$pos < as.numeric(g$end[j])))
  empty <- lengths(in_gene) == 0
  if (any(empty))
    message(sum(empty), " gene(s) in the window hold no genotyped SNP and are left off the ",
            "track under `spacing = \"even\"`, where they have no width: ",
            paste(utils::head(g$name[empty], 5), collapse = ", "),
            if (sum(empty) > 5) ", ..." else "")
  g <- g[!empty, , drop = FALSE]
  in_gene <- in_gene[!empty]
  if (!nrow(g)) return(NULL)
  g$.gene_xmin <- vapply(in_gene, function(i) min(i) - 0.5, numeric(1))
  g$.gene_xmax <- vapply(in_gene, function(i) max(i) + 0.5, numeric(1))
  g
}

# Gene boxes under gapped spacing: interpolated onto the axis, then held out of any SNP
# column the gene does not contain.
#
# The blank columns are real room, so a gene is not restricted to the columns it holds the way
# it is under even spacing -- one sitting in a desert with no genotyped SNP still gets drawn,
# in about the right place, which is most of the point of the mode. But interpolation alone
# would let a gene lying in a *small* gap reach under a neighbouring SNP's tile, and a tile
# drawn over a gene it is not in is the one thing the track must never say. So each box is
# clipped to the space between the flanking columns that are outside it, and a gene with no
# room left is dropped and reported, as under even spacing.
.gapped_gene_boxes <- function(g, sel, x) {
  lo <- .gapped_pos_to_x(as.numeric(g$start), sel$pos, x)
  hi <- .gapped_pos_to_x(as.numeric(g$end), sel$pos, x)
  for (j in seq_len(nrow(g))) {
    gs <- as.numeric(g$start[j]); ge <- as.numeric(g$end[j])
    inside <- which(sel$pos >= gs & sel$pos < ge)
    # a gene always covers the columns it holds, whatever interpolation says
    if (length(inside)) {
      lo[j] <- min(lo[j], min(x[inside]) - 0.5)
      hi[j] <- max(hi[j], max(x[inside]) + 0.5)
    }
    left <- x[sel$pos < gs]
    right <- x[sel$pos >= ge]
    if (length(left)) lo[j] <- max(lo[j], max(left) + 0.5)
    if (length(right)) hi[j] <- min(hi[j], min(right) - 0.5)
  }
  g$.gene_xmin <- lo
  g$.gene_xmax <- hi
  empty <- !(hi > lo)
  if (any(empty))
    message(sum(empty), " gene(s) in the window have no room on the axis under `spacing = ",
            "\"gapped\"` -- no genotyped SNP and no gap to sit in: ",
            paste(utils::head(g$name[empty], 5), collapse = ", "),
            if (sum(empty) > 5) ", ..." else "")
  g <- g[!empty, , drop = FALSE]
  if (!nrow(g)) NULL else g
}

# Keep a genomic-spacing mark from spilling into a gene it is not in: each SNP's fixed-width
# mark is clipped to the gene that holds it, or to the gap between the flanking genes when it
# holds none. Only marks within half their own width of a boundary move at all.
.clip_tiles_to_genes <- function(tiles, sel, boxes, min_frac = 0.25) {
  if (is.null(boxes) || !nrow(boxes)) return(tiles)
  gs <- as.numeric(boxes$.gene_xmin); ge <- as.numeric(boxes$.gene_xmax)
  w <- tiles$xmax - tiles$xmin
  for (k in seq_along(tiles$xmin)) {
    p <- sel$pos[k]
    inside <- which(p >= gs & p < ge)
    if (length(inside)) {
      lo <- max(gs[inside]); hi <- min(ge[inside])
    } else {
      before <- ge[ge <= p]; after <- gs[gs > p]
      lo <- if (length(before)) max(before) else -Inf
      hi <- if (length(after)) min(after) else Inf
    }
    # never shrink a mark to invisibility: keep a quarter of its width, centred on the SNP
    floor_w <- w[k] * min_frac
    if (hi - lo < floor_w) next
    tiles$xmin[k] <- min(max(tiles$xmin[k], lo), hi - floor_w)
    tiles$xmax[k] <- max(min(tiles$xmax[k], hi), lo + floor_w)
  }
  tiles
}

# The dendrogram panel: rows on the y axis so it lines up with the heatmap, distance growing
# leftwards away from it.
.dendro_panel <- function(segs, rows, faceted) {
  if (is.null(segs)) return(NULL)
  # Anchor each facet to the same y range the heatmap uses (its rows, plus the half row the
  # tiles extend past the first and last). Without them the dendrogram facet takes the range
  # of whatever segments it holds, and the leaves drift off the rows they label; a shared
  # global range is just as wrong, since `space = "free_y"` then gives every block the same
  # height regardless of how many samples it has.
  lim <- lapply(split(rows$.row, rows$.split), function(r) c(min(r) - 0.5, max(r) + 0.5))
  anchor <- cbind(data.frame(y = unlist(lim, use.names = FALSE), x = 0),
                  .split_cols(rows, rep(names(lim), each = 2)))

  p <- ggplot2::ggplot(segs) +
    ggplot2::geom_segment(ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y,
                                       yend = .data$yend), linewidth = 0.25,
                          colour = "grey25") +
    ggplot2::geom_blank(data = anchor, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::scale_x_reverse(expand = ggplot2::expansion(mult = c(0.04, 0))) +
    ggplot2::scale_y_reverse(expand = ggplot2::expansion(0)) +
    ggplot2::theme_void() +
    # theme_void() leaves the facet strips, and the blocks are already named once on the far
    # right -- naming them here as well repeats every label and pushes the dendrogram away
    # from the genotypes it is describing
    ggplot2::theme(plot.margin = ggplot2::margin(l = 2, r = 0, t = 0, b = 0),
                   panel.spacing.y = grid::unit(2, "pt"),
                   strip.text = ggplot2::element_blank(),
                   strip.background = ggplot2::element_blank())
  if (faceted) p <- p + .split_facet(rows)
  p
}

# dendrogram | heatmap | annotations, with the gene track under the heatmap. The track and the
# annotation strips keep the absolute sizes they measured for themselves, so their labels
# cannot be clipped however tall or wide the page ends up; the heatmap takes what is left.
#
# patchwork's `design` letters are positional -- A is the FIRST plot passed, B the second --
# so the letters have to be generated in the order the plots go in, not chosen to be mnemonic.
.assemble_haplotype_panels <- function(p, dend, ann, track, dend_width) {
  if (is.null(dend) && is.null(ann) && is.null(track)) return(p)
  .need_package("patchwork", "the panels around plot_region_haplotypes()")

  plots <- list(); widths <- list(); role <- character(0)
  add <- function(pl, w, what) {
    plots[[length(plots) + 1L]] <<- pl
    widths[[length(widths) + 1L]] <<- w
    role <<- c(role, what)
  }
  if (!is.null(dend)) add(dend, grid::unit(dend_width, "null"), "dend")
  add(p, grid::unit(1, "null"), "heat")                       # the one elastic column
  if (!is.null(ann)) add(ann, grid::unit(attr(ann, "ann_in"), "in"), "ann")

  design <- paste(LETTERS[seq_along(plots)], collapse = "")
  heights <- grid::unit(1, "null")
  if (!is.null(track)) {
    # the track spans only the heatmap column: it shares that x axis and nothing else does
    tl <- LETTERS[length(plots) + 1L]
    plots[[length(plots) + 1L]] <- track
    design <- paste(design, paste(ifelse(role == "heat", tl, "#"), collapse = ""), sep = "\n")
    heights <- grid::unit.c(heights, grid::unit(attr(track, "track_in"), "in"))
  }
  out <- patchwork::wrap_plots(plots, design = design) +
    patchwork::plot_layout(widths = do.call(grid::unit.c, widths), heights = heights,
                           guides = "collect")
  attr(out, "plasgenomics_track_in") <- if (is.null(track)) 0 else attr(track, "track_in")
  out
}

# Lay a caller's palette over the shared one: named colours replace those levels and leave
# the rest, so recolouring one level does not mean respelling all of them. An unnamed vector
# is positional, in level order, as elsewhere in the package.
.merge_palette <- function(base, override, levels, what) {
  if (is.null(override) || !length(override)) return(base)
  if (is.null(names(override))) {
    if (length(override) < length(levels))
      stop("`annotation_colours$", what, "` has ", length(override), " colour(s) for ",
           length(levels), " level(s) (", paste(levels, collapse = ", "),
           "); name them, or give one per level", call. = FALSE)
    override <- stats::setNames(override[seq_along(levels)], levels)
  }
  unknown <- setdiff(names(override), levels)
  if (length(unknown))
    warning("`annotation_colours$", what, "` names level(s) that are not in the data: ",
            paste(unknown, collapse = ", "), ". Levels here: ",
            paste(levels, collapse = ", "), call. = FALSE)
  base[names(override)] <- override
  base
}

# A shared colour map can predate a level -- metadata added later, a level renamed -- and
# `scale_fill_manual()` draws an uncoloured key rather than complaining. Fill any gap from
# the automatic palette for that column, so every key drawn has a colour behind it.
.fill_palette_gaps <- function(pal, levels, meta, col) {
  gaps <- setdiff(levels, names(pal))
  if (!length(gaps)) return(pal)
  auto <- tryCatch(meta_colors(meta, cols = col)[[col]], error = function(e) NULL)
  for (g in gaps) pal[g] <- if (!is.null(auto) && g %in% names(auto)) auto[[g]] else "grey85"
  pal
}
