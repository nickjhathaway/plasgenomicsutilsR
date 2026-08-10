# Amino-acid positions -> the genomic interval of their codon.

# GFF attribute value, length-preserving (NA when the key is absent).
.gff_attr <- function(a, key) {
  has <- grepl(paste0("(^|;)", key, "="), a)
  val <- sub(paste0(".*(^|;)", key, "=([^;]*).*"), "\\2", a)
  ifelse(has, val, NA_character_)
}

#' Read the CDS features of a GFF
#'
#' The coding exons only, with the transcript each belongs to -- what
#' [aa_intervals()] needs to walk a protein back onto the genome. Reading the GFF is the slow
#' part, so parse once and reuse the result across calls.
#'
#' @param gff Path to a GFF3 file.
#' @return A data frame of `transcript_id`, `gene_id`, `chrom`, `start`, `end` (1-based
#'   inclusive, as the GFF gives them), `strand` and `phase`, one row per CDS exon.
#' @seealso [aa_intervals()]
#' @export
read_gff_cds <- function(gff) {
  if (length(gff) != 1 || !is.character(gff) || !file.exists(gff))
    stop("`gff` must be the path to a GFF file", call. = FALSE)
  g <- utils::read.delim(gff, header = FALSE, comment.char = "#", quote = "",
                         stringsAsFactors = FALSE,
                         col.names = c("seqid", "source", "type", "start", "end",
                                       "score", "strand", "phase", "attr"))
  cds <- g[g$type == "CDS", , drop = FALSE]
  if (!nrow(cds)) stop("no CDS features in ", gff, call. = FALSE)
  out <- data.frame(
    transcript_id = .gff_attr(cds$attr, "Parent"),
    gene_id = .gff_attr(cds$attr, "gene_id"),
    chrom = cds$seqid,
    start = as.integer(cds$start),
    end = as.integer(cds$end),
    strand = cds$strand,
    # GFF phase: bases to drop from the start of this exon to reach its first whole codon
    phase = suppressWarnings(as.integer(cds$phase)),
    stringsAsFactors = FALSE)
  out$phase[is.na(out$phase)] <- 0L
  if (anyNA(out$transcript_id))
    stop("some CDS rows have no `Parent` attribute naming their transcript", call. = FALSE)
  out[order(out$transcript_id, out$start), , drop = FALSE]
}

# Resolve whatever the caller put in `transcript_id` to real transcript ids: a transcript id,
# a gene id (its transcripts), or a gene symbol from `genes`.
.resolve_transcripts <- function(ids, cds, genes) {
  tx <- unique(cds$transcript_id)
  by_gene <- split(cds$transcript_id, cds$gene_id)
  sym <- NULL
  if (!is.null(genes)) {
    g <- as.data.frame(genes)
    if (all(c("name", "gene_id") %in% names(g)))
      sym <- stats::setNames(as.character(g$gene_id), tolower(as.character(g$name)))
  }
  lapply(ids, function(id) {
    if (id %in% tx) return(id)
    gene <- if (id %in% names(by_gene)) id
            else if (!is.null(sym) && tolower(id) %in% names(sym)) sym[[tolower(id)]]
            else NA_character_
    if (is.na(gene) || !gene %in% names(by_gene)) return(character(0))
    unique(by_gene[[gene]])
  })
}

# Genomic positions of the three bases of codon `k`, in transcript order.
.codon_positions <- function(ex, k) {
  # transcript order: along the strand, so a minus-strand transcript starts at its highest
  # coordinate. Phase on that first exon says how many bases precede the first whole codon.
  ex <- ex[order(ex$start, decreasing = ex$strand[1] == "-"), , drop = FALSE]
  len <- ex$end - ex$start + 1L
  skip <- ex$phase[1]
  ends <- cumsum(len)
  starts <- ends - len + 1L
  t_at <- skip + (3L * (k - 1L) + 1L):(3L * k)     # 1-based along the concatenated CDS
  if (any(t_at > sum(len))) return(NULL)           # past the end of the protein
  vapply(t_at, function(t) {
    j <- which(t >= starts & t <= ends)[1]
    off <- t - starts[j]
    if (ex$strand[j] == "-") ex$end[j] - off else ex$start[j] + off
  }, numeric(1))
}

#' Genomic interval of an amino-acid position
#'
#' Turns "codon 76 of *pfcrt*" into a genomic interval, by walking the protein back through
#' the transcript's coding exons. Amino-acid positions are how resistance markers are named
#' and reported, while every plot and interval tool here works in genomic coordinates, and the
#' conversion is not something you can do by eye: it depends on the exon structure, the strand
#' and the CDS phase.
#'
#' The result is shaped like the package's other interval tables -- `chr`, `start`, `end`,
#' `name` -- so it drops straight into `genes =`, `mark_snps =`, [annotate_snps()] or
#' [bed_intersect()]. `name` is `<transcript_id>-AA<aa_position>`.
#'
#' A codon can straddle an intron, in which case its three bases are not contiguous: `start`
#' and `end` then span the intron as well, and `spans_intron` flags it so the width is not
#' mistaken for three bases. `codon_positions` always lists the three base positions
#' themselves.
#'
#' @param positions A data frame with `transcript_id` and `aa_position`. `transcript_id` may
#'   be a transcript id (`"PF3D7_0709000.1"`), a gene id (`"PF3D7_0709000"` -- every
#'   transcript of it is returned), or a gene symbol when `genes` is given (`"pfcrt"`).
#' @param gff A GFF path, or the result of [read_gff_cds()] (parse once, reuse).
#' @param genes Optional gene table with `name` and `gene_id` columns, so `transcript_id` can
#'   be a symbol; defaults to [PF3D7_GENES].
#' @param one_based_output Return 1-based inclusive coordinates instead of the package's
#'   0-based half-open convention. `FALSE` (default) keeps `start` 0-based so the output can
#'   be used as an interval table directly.
#' @return A data frame with `chr` (normalised), `chrom` (as the GFF spells it), `start`,
#'   `end`, `name`, `transcript_id`, `gene_id`, `aa_position`, `strand`, `codon_positions`
#'   (a comma-separated list of the three base positions, always 1-based as coordinates are
#'   usually quoted), `n_exons` and `spans_intron`. Positions past the end of a protein, and
#'   ids that match nothing, are dropped with a warning naming them.
#' @seealso [read_gff_cds()], [annotate_snps()]
#' @examples
#' \dontrun{
#' cds <- read_gff_cds("Pf3D7.gff")
#' aa <- data.frame(transcript_id = c("pfcrt", "pfkelch13"), aa_position = c(76, 580))
#' aa_intervals(aa, cds)
#' }
#' @export
aa_intervals <- function(positions, gff, genes = PF3D7_GENES, one_based_output = FALSE) {
  p <- as.data.frame(positions)
  need <- c("transcript_id", "aa_position")
  if (!all(need %in% names(p)))
    stop("`positions` needs `transcript_id` and `aa_position` columns", call. = FALSE)
  p$transcript_id <- as.character(p$transcript_id)
  p$aa_position <- suppressWarnings(as.integer(p$aa_position))
  if (anyNA(p$aa_position) || any(p$aa_position < 1, na.rm = TRUE))
    stop("`aa_position` must be a whole number 1 or greater", call. = FALSE)

  cds <- if (is.data.frame(gff)) gff else read_gff_cds(gff)
  if (!all(c("transcript_id", "chrom", "start", "end", "strand", "phase") %in% names(cds)))
    stop("`gff` must be a GFF path or a read_gff_cds() result", call. = FALSE)
  by_tx <- split(cds, cds$transcript_id)

  resolved <- .resolve_transcripts(unique(p$transcript_id), cds, genes)
  names(resolved) <- unique(p$transcript_id)
  unknown <- names(resolved)[!lengths(resolved)]
  if (length(unknown))
    warning("no transcript found for: ", paste(unknown, collapse = ", "), call. = FALSE)

  rows <- list(); past_end <- character(0)
  for (i in seq_len(nrow(p))) {
    for (tx in resolved[[p$transcript_id[i]]]) {
      ex <- by_tx[[tx]]
      pos <- .codon_positions(ex, p$aa_position[i])
      if (is.null(pos)) {
        past_end <- c(past_end, paste0(tx, " AA", p$aa_position[i]))
        next
      }
      rows[[length(rows) + 1L]] <- data.frame(
        chr = normalise_chr(ex$chrom[1]), chrom = ex$chrom[1],
        start = min(pos), end = max(pos),
        name = paste0(tx, "-AA", p$aa_position[i]),
        transcript_id = tx, gene_id = ex$gene_id[1],
        aa_position = p$aa_position[i], strand = ex$strand[1],
        codon_positions = paste(sort(pos), collapse = ","),
        n_exons = length(unique(vapply(pos, function(v)
          which(v >= ex$start & v <= ex$end)[1], integer(1)))),
        stringsAsFactors = FALSE)
    }
  }
  if (length(past_end))
    warning("past the end of the protein: ", paste(past_end, collapse = ", "), call. = FALSE)
  if (!length(rows))
    stop("no amino-acid position could be placed on the genome", call. = FALSE)

  out <- do.call(rbind, rows)
  out$spans_intron <- out$n_exons > 1L
  # 0-based half-open like every other interval here, unless asked otherwise
  if (!one_based_output) out$start <- out$start - 1L
  rownames(out) <- NULL
  out
}
