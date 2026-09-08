# A genomic interval -> the residues it covers, per transcript. The companion of
# snp_aa_positions() for regions rather than points: which amino acids does a target window,
# an amplicon, or a piece of gene left after masking actually read?

#' The amino acids a genomic range covers
#'
#' For each interval, every transcript whose CDS it overlaps and the range of codons of
#' that transcript the interval covers -- the residues a target window, an amplicon, or a
#' piece of gene left after [bed_subtract()] will actually read. Covered means the codon's
#' bases lie inside the interval, so a window that starts 100 bp upstream of a gene reports
#' from codon 1, and one that starts inside the gene reports from the first codon it holds
#' whole; nothing is inferred from the interval's first and last positions.
#'
#' A codon is covered when all three of its bases are inside the interval. `partial = TRUE`
#' also counts a codon the interval holds only one or two bases of -- the codon it starts or
#' ends in -- and `partial_start` / `partial_end` say when that happened. Because an interval
#' is contiguous and a transcript's coding bases are in order along it, the covered codons of
#' one transcript are always a contiguous run, `aa_start` to `aa_end`.
#'
#' Codon numbers are **1-based**, counting the initiator methionine as 1, as everywhere in
#' the package ([plasgenomicsutilsR-coordinates]). A CDS that includes its stop codon reports
#' it as the last codon (425 for *pfcrt*, whose protein is 424 residues), as
#' [snp_aa_positions()] does; in `aa_seq` it is `*`. Genomic coordinates in `ranges` are 0-based
#' half-open, like every interval table.
#'
#' @param ranges An interval table (`chr` or `chrom`, `start`, `end`), such as
#'   [bed_subtract()] or [read_gff_features()] return, or a gene table.
#' @param gff Path or URL to a GFF3, or the table [read_gff_cds()] returns -- parse once and
#'   reuse it.
#' @param keep `"all"` (the default) keeps every input row, filling the result columns with
#'   `NA` for an interval that covers no codon; `"hits"` drops those rows.
#' @param partial Count a codon the interval covers only partly (default `FALSE`).
#' @param fasta Optional sequence, as for [snp_aa_positions()]: a path or URL to the genome
#'   (gzipped is fine), or a named vector of sequences. When given, `aa_seq` holds the
#'   reference residues of the covered codons, `aa_start` to `aa_end`, in transcript order;
#'   a stop is `*` and an unreadable codon `X`. A GFF that ends in a `##FASTA` section
#'   supplies its own sequence.
#' @param genes A gene table with `name` and `gene_id`, used to fill `gene_name`. Defaults to
#'   [PF3D7_GENES]; `NULL` leaves `gene_name` `NA`.
#' @param chrom,start,end Column names in `ranges` (defaults `"chr"`, `"start"`, `"end"`;
#'   `"chrom"` is accepted as a chromosome alias).
#' @return A tibble with one row per interval and covered transcript: every column of
#'   `ranges`, then `transcript_id`, `gene_id`, `gene_name`, `aa_start`, `aa_end`, `n_aa`,
#'   `partial_start`, `partial_end`, `strand`, `coding` (whether any codon was covered), and
#'   `aa_seq` when sequence is available. Rows follow the input order; an interval covering
#'   several transcripts has several rows.
#' @seealso [snp_aa_positions()] for single positions, [aa_intervals()] for the other
#'   direction, [bed_subtract()] and [read_gff_features()] for where the intervals come from.
#' @examples
#' gff <- system.file("extdata", "pf3d7_drug_gene_cds.gff", package = "plasgenomicsutilsR")
#' cds <- read_gff_cds(gff)
#'
#' # pfcrt with its tandem repeats removed: which residues does each piece still read?
#' crt <- PF3D7_GENES[PF3D7_GENES$name == "pfcrt", ]
#' pieces <- bed_subtract(crt, tandem_repeats_to_avoid(pf3d7_tandem_repeats()), pad = 10)
#' genomic_range_aa_positions(pieces, cds)[, c("piece", "start", "end", "aa_start", "aa_end",
#'                                             "n_aa")]
#'
#' # the K76T codon itself, and a window holding only its middle base
#' k76 <- data.frame(chr = "Pf3D7_07_v3", start = c(403623, 403624), end = c(403626, 403625))
#' genomic_range_aa_positions(k76, cds)[, c("start", "end", "aa_start", "aa_end", "coding")]
#' genomic_range_aa_positions(k76, cds, partial = TRUE)[, c("start", "end", "aa_start",
#'                                                          "partial_start", "partial_end")]
#'
#' # with sequence, the residues themselves
#' fa <- system.file("extdata", "pf3d7_drug_gene_regions.fasta.gz", package = "plasgenomicsutilsR")
#' genomic_range_aa_positions(k76[1, ], cds, fasta = fa)$aa_seq
#' @export
genomic_range_aa_positions <- function(ranges, gff, keep = c("all", "hits"), partial = FALSE,
                                       fasta = NULL, genes = PF3D7_GENES,
                                       chrom = "chr", start = "start", end = "end") {
  keep <- match.arg(keep)
  df <- as.data.frame(ranges, stringsAsFactors = FALSE)
  if (!nrow(df)) return(tibble::as_tibble(df))
  iv <- .as_interval_table(df, chrom, start, end, what = "ranges")
  if (anyNA(iv$start) || anyNA(iv$end) || any(iv$end < iv$start))
    stop("`ranges` has a missing or inverted start/end", call. = FALSE)
  chr <- normalise_chr(iv$chr)
  s1 <- iv$start + 1; e1 <- iv$end                 # the GFF is 1-based inclusive

  parsed <- if (is.data.frame(gff)) gff else read_gff_cds(gff)
  seqs <- .resolve_sequence(fasta, parsed)
  cds <- as.data.frame(parsed, stringsAsFactors = FALSE)
  cds$.chr <- normalise_chr(cds$chrom)

  name_of <- NULL
  if (!is.null(genes)) {
    gt <- as.data.frame(genes)
    if (all(c("name", "gene_id") %in% names(gt)))
      name_of <- stats::setNames(as.character(gt$name), as.character(gt$gene_id))
  }

  hits <- list()
  for (ch in intersect(unique(chr), unique(cds$.chr))) {
    ex_chr <- cds[cds$.chr == ch, , drop = FALSE]
    ri <- which(chr == ch)
    for (tx in unique(ex_chr$transcript_id)) {
      o <- .tx_offsets(ex_chr[ex_chr$transcript_id == tx, , drop = FALSE])
      ex <- o$ex
      cand <- ri[s1[ri] <= max(ex$end) & e1[ri] >= min(ex$start)]
      for (k in cand) {
        # every coding base of this transcript inside the interval, as a codon number
        aa_all <- integer(0)
        for (j in seq_len(nrow(ex))) {
          lo <- max(s1[k], ex$start[j]); hi <- min(e1[k], ex$end[j])
          if (lo > hi) next
          p <- seq.int(lo, hi)
          off <- if (ex$strand[j] == "-") ex$end[j] - p else p - ex$start[j]
          t_eff <- o$cum_start[j] + off - o$skip      # 1-based from the first whole codon
          t_eff <- t_eff[t_eff >= 1]
          aa_all <- c(aa_all, ceiling(t_eff / 3))
        }
        if (!length(aa_all)) next
        tab <- table(aa_all)
        aa_idx <- as.integer(names(tab))
        full <- as.integer(tab) >= 3L
        covered <- if (isTRUE(partial)) aa_idx else aa_idx[full]
        if (!length(covered)) next
        a0 <- min(covered); a1 <- max(covered)
        hit <- data.frame(
          .row = k, transcript_id = tx, gene_id = ex$gene_id[1],
          aa_start = as.integer(a0), aa_end = as.integer(a1), n_aa = length(covered),
          partial_start = !full[aa_idx == a0], partial_end = !full[aa_idx == a1],
          strand = ex$strand[1], stringsAsFactors = FALSE)
        if (!is.null(seqs)) {
          aa <- seq.int(a0, a1)
          t1 <- (aa - 1) * 3 + 1 + o$skip
          rc <- .ref_codons(rep(ex$chrom[1], length(aa)), .tx_to_genomic(o, t1),
                            .tx_to_genomic(o, t1 + 1), .tx_to_genomic(o, t1 + 2),
                            rep(ex$strand[1], length(aa)), seqs)
          res <- rc$aa
          res[is.na(res)] <- "X"
          hit$aa_seq <- paste(res, collapse = "")
        }
        hits[[length(hits) + 1L]] <- hit
      }
    }
  }
  empty <- data.frame(.row = integer(0), transcript_id = character(0), gene_id = character(0),
                      aa_start = integer(0), aa_end = integer(0), n_aa = integer(0),
                      partial_start = logical(0), partial_end = logical(0),
                      strand = character(0), stringsAsFactors = FALSE)
  if (!is.null(seqs)) empty$aa_seq <- character(0)
  ann <- if (length(hits)) do.call(rbind, hits) else empty

  if (keep == "hits") {
    if (!nrow(ann)) return(tibble::as_tibble(df[0, , drop = FALSE]))
  } else {
    miss <- setdiff(seq_len(nrow(df)), ann$.row)
    if (length(miss)) {
      filler <- empty[rep(1L, 0), , drop = FALSE]
      filler <- data.frame(.row = miss, transcript_id = NA_character_, gene_id = NA_character_,
                           aa_start = NA_integer_, aa_end = NA_integer_, n_aa = NA_integer_,
                           partial_start = NA, partial_end = NA, strand = NA_character_,
                           stringsAsFactors = FALSE)
      if (!is.null(seqs)) filler$aa_seq <- NA_character_
      ann <- rbind(ann, filler)
    }
  }
  ann <- ann[order(ann$.row, ann$transcript_id), , drop = FALSE]
  out <- cbind(df[ann$.row, , drop = FALSE], ann[, -1, drop = FALSE])
  out$gene_name <- if (is.null(name_of)) NA_character_ else unname(name_of[out$gene_id])
  out$coding <- !is.na(out$transcript_id)
  lead <- c("transcript_id", "gene_id", "gene_name", "aa_start", "aa_end", "n_aa",
            "partial_start", "partial_end", "strand", "coding")
  out <- out[, c(setdiff(names(out), c(lead, "aa_seq")), lead,
                 if (!is.null(seqs)) "aa_seq"), drop = FALSE]
  rownames(out) <- NULL
  tibble::as_tibble(out)
}
