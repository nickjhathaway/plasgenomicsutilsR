# Reading a GFF once, and pulling a gene's exons, introns, CDS or span out of it as
# 0-based half-open intervals -- the pieces a target design starts from before the tandem
# repeats come out (gene -> exons -> bed_subtract(..., pad = )).

# Parse a GFF3 (path or URL, plain or gzipped) into its records, and any ##FASTA section.
.read_gff <- function(gff) {
  if (length(gff) != 1 || !is.character(gff))
    stop("`gff` must be one path or URL to a GFF file", call. = FALSE)
  is_url <- grepl("^(https?|ftp)://", gff)
  if (!is_url && !file.exists(gff))
    stop("no such file: ", gff, call. = FALSE)
  cols <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attr")
  if (is_url) {
    # read.delim() gunzips a local path on its own, but a remote stream has to be decompressed
    # here -- and read through readLines(), since read.table() cannot push back on the
    # binary-mode connection that gzcon() gives.
    con <- if (grepl("\\.gz$", gff)) gzcon(url(gff, open = "rb")) else url(gff)
    on.exit(try(close(con), silent = TRUE), add = TRUE)
    txt <- readLines(con, warn = FALSE)
  } else {
    txt <- readLines(gff, warn = FALSE)
  }
  # A GFF3 may end with a `##FASTA` directive and the sequences themselves. Those lines are not
  # tab-delimited records, so they have to come off before parsing -- and they are worth
  # keeping, since they are the reference bases snp_aa_positions() needs.
  fa_at <- which(grepl("^##FASTA", txt))
  seqs <- NULL
  if (length(fa_at)) {
    seqs <- .parse_fasta(txt[seq(fa_at[1] + 1L, length(txt))])
    txt <- txt[seq_len(fa_at[1] - 1L)]
  }
  g <- utils::read.delim(text = txt, header = FALSE, comment.char = "#", quote = "",
                         stringsAsFactors = FALSE, col.names = cols)
  list(records = g, sequence = seqs)
}

# Ensembl prefixes ids with their type (`transcript:PF3D7_0709000.1`, `gene:...`)
.gff_id <- function(v) sub("^[A-Za-z_]+:", "", v)

# What the caller named -> the gene ids to keep, and any transcripts named directly.
.resolve_gene_ids <- function(ids, gene_ids, gff_name, tx2gene, genes) {
  sym <- NULL
  if (!is.null(genes)) {
    g <- as.data.frame(genes)
    if (all(c("name", "gene_id") %in% names(g)))
      sym <- split(as.character(g$gene_id), tolower(as.character(g$name)))
  }
  has_name <- !is.na(gff_name)
  gsym <- split(gene_ids[has_name], tolower(gff_name[has_name]))
  want <- character(); tx <- character(); tx_gene <- character(); miss <- character()
  for (id in unique(as.character(ids))) {
    lid <- tolower(id)
    if (id %in% gene_ids) {
      want <- c(want, id)
    } else if (!is.null(sym) && lid %in% names(sym) && any(sym[[lid]] %in% gene_ids)) {
      want <- c(want, intersect(sym[[lid]], gene_ids))
    } else if (lid %in% names(gsym)) {
      want <- c(want, gsym[[lid]])
    } else if (id %in% names(tx2gene)) {
      tx <- c(tx, id); tx_gene <- c(tx_gene, unname(tx2gene[[id]]))
    } else {
      miss <- c(miss, id)
    }
  }
  if (length(miss))
    stop("not a gene id, gene `Name` or transcript id in the GFF",
         if (!is.null(sym)) ", nor a `name` in `genes`" else "", ": ",
         paste(miss, collapse = ", "), call. = FALSE)
  list(genes = unique(c(want, tx_gene)), tx = unique(tx),
       tx_only = setdiff(unique(tx_gene), want))
}

#' Exons, introns, CDS or spans of genes, from a GFF
#'
#' The intervals a target design starts from. Most of a gene's troublesome tandem repeats
#' sit in its introns, so a design will often want the exons alone; the rest are in the
#' exons, so the flow is gene, then exons, then [bed_subtract()] with the repeats from
#' [tandem_repeats_to_avoid()]. This reads a GFF3 and returns the pieces you ask for as
#' interval rows that drop straight into [bed_subtract()], [bed_intersect()] and
#' [write_bed()].
#'
#' `feature` is what to return:
#' * `"exon"` -- the exons (the GFF's `exon` features; a GFF with none falls back to its
#'   `CDS` features, with a message);
#' * `"intron"` -- the gaps between consecutive exons;
#' * `"cds"` -- the coding parts of the exons only;
#' * `"gene"` -- each gene's whole span, one row per gene.
#'
#' `per = "gene"` (the default) gives one set of intervals per gene, merging the exons of
#' every isoform, so an intron is a stretch no isoform transcribes; `per = "transcript"`
#' keeps each transcript's own exons and introns. `index` numbers them in transcript
#' orientation, so exon 1 of a minus-strand gene is the one with the highest coordinates.
#'
#' @section Coordinates: **0-based half-open**, like every interval table in the package
#'   (see [plasgenomicsutilsR-coordinates]) -- converted from the GFF's 1-based inclusive
#'   coordinates here, at the boundary. [read_gff_cds()] is the one reader that keeps the
#'   GFF's numbering, because [aa_intervals()] walks codons on it; do not mix the two.
#'
#' @param gff Path or URL to a GFF3 file, plain or gzipped, as for [read_gff_cds()].
#' @param feature One of `"exon"`, `"intron"`, `"cds"`, `"gene"`.
#' @param ids Which genes: gene ids (`"PF3D7_0709000"`), the GFF's `Name` attribute
#'   (`"CRT"`, case-insensitive), a `name` from `genes` (`"pfcrt"`), or transcript ids, which
#'   restrict that gene to the transcript named. `NULL` (the default) returns every gene.
#' @param per `"gene"` merges isoforms; `"transcript"` keeps them apart.
#' @param genes A gene table with `name` and `gene_id` columns, so `ids` can use the
#'   package's friendly names. Defaults to [PF3D7_GENES]; `NULL` to use the GFF's `Name`
#'   only.
#' @return A tibble with one row per interval: `gene_id`, `name` (the `genes` name, else
#'   the GFF `Name`, else `NA`), `transcript_id` when `per = "transcript"`, `feature`,
#'   `index`, `chrom` (as the GFF spells it), `chr` (normalised), `start`, `end`, `width`,
#'   `strand`. Genes in GFF order, intervals by position.
#' @seealso [tandem_repeats_to_avoid()] and [bed_subtract()] for what comes next,
#'   [read_gff_cds()] for the codon-walking reader.
#' @examples
#' gff <- system.file("extdata", "pf3d7_drug_gene_cds.gff", package = "plasgenomicsutilsR")
#' read_gff_features(gff, "exon", ids = "pfcrt")
#' read_gff_features(gff, "intron", ids = "pfcrt")[, c("index", "start", "end", "width")]
#' read_gff_features(gff, "gene")
#'
#' # the exons of two genes, minus the tandem repeats and 10 bp of clearance around each
#' exons <- read_gff_features(gff, "exon", ids = c("pfcrt", "pfkelch13"))
#' targets <- bed_subtract(exons, tandem_repeats_to_avoid(pf3d7_tandem_repeats()), pad = 10)
#' targets[, c("name", "index", "start", "end", "piece", "width")]
#' @export
read_gff_features <- function(gff, feature = c("exon", "intron", "cds", "gene"), ids = NULL,
                              per = c("gene", "transcript"), genes = PF3D7_GENES) {
  feature <- match.arg(feature)
  per <- match.arg(per)
  g <- .read_gff(gff)$records
  id <- .gff_id(.gff_attr(g$attr, "ID"))
  parent <- .gff_id(.gff_attr(g$attr, "Parent"))
  is_gene <- grepl("gene$", g$type)            # gene, protein_coding_gene, ncRNA_gene, pseudogene
  if (!any(is_gene)) stop("no gene features in ", gff, call. = FALSE)
  gene_ids <- id[is_gene]
  gff_name <- stats::setNames(.gff_attr(g$attr[is_gene], "Name"), gene_ids)
  tx_rows <- which(!is_gene & !is.na(parent) & parent %in% gene_ids)
  tx2gene <- stats::setNames(parent[tx_rows], id[tx_rows])

  if (feature == "gene") {
    rows <- which(is_gene)
  } else {
    type <- if (feature == "cds") "CDS" else "exon"
    rows <- which(g$type == type)
    if (!length(rows) && type == "exon") {
      rows <- which(g$type == "CDS")
      if (length(rows)) message("no exon features in the GFF; using its CDS features as the exons")
    }
    if (!length(rows)) stop("no ", type, " features in ", gff, call. = FALSE)
  }

  if (feature == "gene") {
    gid <- id[rows]
    tx <- rep(NA_character_, length(rows))
  } else {
    # VEuPathDB writes a `gene_id` on every exon; Ensembl only a `Parent=transcript:...`, so
    # walk transcript -> gene, and failing that strip the `.N` isoform suffix
    tx <- parent[rows]
    gid <- .gff_attr(g$attr[rows], "gene_id")
    i <- is.na(gid) & tx %in% gene_ids; gid[i] <- tx[i]
    i <- is.na(gid) & tx %in% names(tx2gene); gid[i] <- unname(tx2gene[tx[i]])
    i <- is.na(gid); gid[i] <- sub("\\.[0-9]+$", "", tx[i])
  }
  feat <- data.frame(gene_id = gid, transcript_id = tx, chrom = g$seqid[rows],
                     start = as.numeric(g$start[rows]) - 1, end = as.numeric(g$end[rows]),
                     strand = g$strand[rows], stringsAsFactors = FALSE)

  if (!is.null(ids)) {
    want <- .resolve_gene_ids(ids, gene_ids, gff_name, tx2gene, genes)
    keep <- feat$gene_id %in% want$genes &
      (is.na(feat$transcript_id) | !feat$gene_id %in% want$tx_only |
         feat$transcript_id %in% want$tx)
    feat <- feat[keep, , drop = FALSE]
  }
  if (!nrow(feat)) stop("no ", feature, " intervals for the genes asked for", call. = FALSE)

  by_tx <- per == "transcript" && feature != "gene"
  key <- if (by_tx) paste(feat$gene_id, feat$transcript_id) else feat$gene_id
  parts <- lapply(split(seq_len(nrow(feat)), factor(key, levels = unique(key))), function(i) {
    f <- feat[i, , drop = FALSE]
    if (by_tx) {
      o <- order(f$start, f$end); s <- f$start[o]; e <- f$end[o]
    } else {
      m <- .merge_spans(f$start, f$end); s <- m$s; e <- m$e
    }
    if (feature == "intron") {
      if (length(s) < 2L) return(NULL)
      s2 <- e[-length(e)]; e <- s[-1]; s <- s2
    }
    idx <- if (identical(f$strand[1], "-")) rev(seq_along(s)) else seq_along(s)
    data.frame(gene_id = f$gene_id[1], transcript_id = f$transcript_id[1], feature = feature,
               index = idx, chrom = f$chrom[1], start = s, end = e, strand = f$strand[1],
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, parts)
  if (is.null(out) || !nrow(out))
    stop("no ", feature, " intervals: ", if (feature == "intron") "every gene asked for has a single exon"
         else "nothing matched", call. = FALSE)

  nm <- unname(gff_name[out$gene_id])
  if (!is.null(genes)) {
    gt <- as.data.frame(genes)
    if (all(c("name", "gene_id") %in% names(gt))) {
      hit <- match(out$gene_id, as.character(gt$gene_id))
      nm[!is.na(hit)] <- as.character(gt$name)[hit[!is.na(hit)]]
    }
  }
  out$name <- nm
  out$chr <- normalise_chr(out$chrom)
  out$width <- out$end - out$start
  out <- out[order(match(out$gene_id, gene_ids), out$transcript_id, out$start), , drop = FALSE]
  cols <- c("gene_id", "name", if (by_tx) "transcript_id", "feature", "index", "chrom", "chr",
            "start", "end", "width", "strand")
  rownames(out) <- NULL
  tibble::as_tibble(out[, cols, drop = FALSE])
}
