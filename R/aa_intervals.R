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
#' @param gff Path to a GFF3 file, or a URL -- a plain or gzipped file is read straight from
#'   the web, so a released annotation can be used without keeping a copy.
#' @return A tibble of `transcript_id`, `gene_id`, `chrom`, `start`, `end` (1-based
#'   inclusive, as the GFF gives them), `strand` and `phase`, one row per CDS exon.
#' @seealso [aa_intervals()]
#' @examples
#' \dontrun{
#' # VEuPathDB / PlasmoDB, the source the bundled gene datasets were built from
#' cds <- read_gff_cds(paste0("https://plasmodb.org/common/downloads/Current_Release/",
#'                            "Pfalciparum3D7/gff/data/PlasmoDB-68_Pfalciparum3D7.gff"))
#' # Ensembl Protists works too, despite naming its attributes differently
#' cds <- read_gff_cds(paste0("https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists/current/",
#'                            "gff3/plasmodium_falciparum/",
#'                            "Plasmodium_falciparum.GCA000002765v3.63.gff3.gz"))
#' }
#' @export
read_gff_cds <- function(gff) {
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
    g <- utils::read.delim(text = txt, header = FALSE, comment.char = "#", quote = "",
                           stringsAsFactors = FALSE, col.names = cols)
  } else {
    g <- utils::read.delim(gff, header = FALSE, comment.char = "#", quote = "",
                           stringsAsFactors = FALSE, col.names = cols)
  }
  cds <- g[g$type == "CDS", , drop = FALSE]
  if (!nrow(cds)) stop("no CDS features in ", gff, call. = FALSE)

  # VEuPathDB writes `Parent=PF3D7_0709000.1` and a `gene_id`; Ensembl writes
  # `Parent=transcript:PF3D7_0709000.1` and no gene_id at all, so strip the type prefix and
  # fall back to the transcript id minus its `.N` suffix.
  tx <- sub("^(transcript|mRNA|rna|gene):", "", .gff_attr(cds$attr, "Parent"))
  gid <- .gff_attr(cds$attr, "gene_id")
  gid[is.na(gid)] <- sub("\\.[0-9]+$", "", tx[is.na(gid)])
  out <- data.frame(
    transcript_id = tx,
    gene_id = gid,
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
  tibble::as_tibble(out[order(out$transcript_id, out$start), , drop = FALSE])
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
#' [bed_intersect()]. `name` is `<transcript_id>-AA<aa_position>`. Those genomic bounds are
#' 0-based half-open like every other interval here, while the amino-acid position itself is
#' 1-based and `codon_positions` lists 1-based bases -- the deliberate exception described in
#' `?"plasgenomicsutilsR-coordinates"`.
#'
#' A codon can straddle an intron, in which case its three bases are not contiguous: `start`
#' and `end` then span the intron as well, and `spans_intron` flags it so the width is not
#' mistaken for three bases. `codon_positions` always lists the three base positions
#' themselves.
#'
#' @param positions A data frame with `transcript_id` and `aa_position`. `aa_position` is
#'   **1-based**, counting the initiator methionine as 1, matching how residues are numbered in
#'   the literature -- so `76` is the residue everyone calls 76. `transcript_id` may
#'   be a transcript id (`"PF3D7_0709000.1"`), a gene id (`"PF3D7_0709000"` -- every
#'   transcript of it is returned), or a gene symbol when `genes` is given (`"pfcrt"`).
#' @param gff A GFF path, or the result of [read_gff_cds()] (parse once, reuse).
#' @param genes Optional gene table with `name` and `gene_id` columns, so `transcript_id` can
#'   be a symbol; defaults to [PF3D7_GENES].
#' @param one_based_output Return 1-based inclusive coordinates instead of the package's
#'   0-based half-open convention. `FALSE` (default) keeps `start` 0-based so the output can
#'   be used as an interval table directly.
#' @return A tibble with `chr` (normalised), `chrom` (as the GFF spells it), `start`,
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
  tibble::as_tibble(out)
}

# Transcript-relative base index (1-based along the concatenated CDS, before phase) for a
# vector of genomic positions known to lie in exon `j` of `ex`.
.tx_offsets <- function(ex) {
  ex <- ex[order(ex$start, decreasing = ex$strand[1] == "-"), , drop = FALSE]
  len <- ex$end - ex$start + 1L
  ends <- cumsum(len)
  list(ex = ex, cum_start = ends - len + 1L, skip = ex$phase[1], total = sum(len))
}

#' The amino acid a SNP falls in
#'
#' The reverse of [aa_intervals()]: given SNP positions and a GFF, which codon of which
#' transcript each one sits in. That is how a variant gets talked about -- "*pfdhps* A437G" --
#' and it is the step between "these SNPs are in this gene" and "these are the amino acids they
#' change".
#'
#' `aa_position` is **1-based**, counting the initiator methionine as 1, because that is how the
#' literature numbers residues. `codon_base` says which of the codon's three bases the SNP is,
#' in transcript orientation, so on a minus-strand gene `codon_base == 1` is the *highest*
#' genomic coordinate of the three.
#'
#' A SNP outside any CDS gets `NA` (or is dropped by `keep = "hits"`) -- introns, UTRs and
#' intergenic space are all simply non-coding here. A SNP inside overlapping isoforms yields one
#' row per transcript, since the codon it hits can differ between them.
#'
#' Being in a codon says nothing about whether the residue actually changes: that needs the
#' alleles and the reading frame's other two bases, which are not part of this.
#'
#' @param snps A data frame with `snp_id` (`"chr:pos"`) or `chr` and `pos` columns; any other
#'   columns are carried through.
#' @param gff A GFF path, or the result of [read_gff_cds()] (parse once, reuse).
#' @param keep `"all"` (default) keeps non-coding SNPs with `NA` annotations; `"hits"` keeps
#'   only those in a CDS.
#' @param one_based_snps The positions in `snps` are 1-based (VCF `POS`). Genotype-matrix column
#'   names from [load_genotypes()] are, so set this when feeding those in; the package's own
#'   tables are 0-based, which is the default.
#' @return `snps` with `transcript_id`, `gene_id`, `aa_position` (1-based), `codon_base`
#'   (1/2/3 in transcript orientation), `strand` and `coding` added.
#' @seealso [aa_intervals()] for the other direction, [annotate_snps()] to first ask which
#'   gene a SNP is in.
#' @examples
#' \dontrun{
#' cds <- read_gff_cds("Pf3D7.gff")
#' # every genotyped SNP in pfdhps, and the residues they sit on
#' snps <- data.frame(snp_id = colnames(ps$genotype("full")))
#' snp_aa_positions(snps, cds, keep = "hits", one_based_snps = TRUE)
#' }
#' @export
snp_aa_positions <- function(snps, gff, keep = c("all", "hits"), one_based_snps = FALSE) {
  keep <- match.arg(keep)
  df <- as.data.frame(snps, stringsAsFactors = FALSE)
  if (!nrow(df)) return(snps)

  if (all(c("chr", "pos") %in% names(df))) {
    chr <- normalise_chr(df$chr)
    pos <- as.numeric(df$pos)
  } else if ("snp_id" %in% names(df)) {
    id <- as.character(df$snp_id)
    chr <- normalise_chr(sub(":[^:]*$", "", id))
    pos <- suppressWarnings(as.numeric(sub("^.*:", "", id)))
    if (anyNA(pos))
      stop("could not read a position out of `snp_id`; expected \"chr:pos\"", call. = FALSE)
  } else {
    stop("`snps` needs a `snp_id` column, or `chr` and `pos` columns", call. = FALSE)
  }
  # the GFF is 1-based inclusive, so compare in those terms
  if (!isTRUE(one_based_snps)) pos <- pos + 1

  cds <- if (is.data.frame(gff)) gff else read_gff_cds(gff)
  cds <- as.data.frame(cds, stringsAsFactors = FALSE)
  cds$.chr <- normalise_chr(cds$chrom)

  hits <- list()
  for (ch in intersect(unique(chr), unique(cds$.chr))) {
    ex_chr <- cds[cds$.chr == ch, , drop = FALSE]
    si <- which(chr == ch)
    for (tx in unique(ex_chr$transcript_id)) {
      o <- .tx_offsets(ex_chr[ex_chr$transcript_id == tx, , drop = FALSE])
      ex <- o$ex
      for (j in seq_len(nrow(ex))) {
        # SNPs of this chromosome inside this exon
        k <- si[pos[si] >= ex$start[j] & pos[si] <= ex$end[j]]
        if (!length(k)) next
        off <- if (ex$strand[j] == "-") ex$end[j] - pos[k] else pos[k] - ex$start[j]
        t_eff <- o$cum_start[j] + off - o$skip        # 1-based from the first whole codon
        aa <- ifelse(t_eff >= 1, ceiling(t_eff / 3), NA_real_)
        base <- ifelse(t_eff >= 1, ((t_eff - 1) %% 3) + 1, NA_real_)
        hits[[length(hits) + 1L]] <- data.frame(
          .row = k, transcript_id = tx, gene_id = ex$gene_id[j],
          aa_position = as.integer(aa), codon_base = as.integer(base),
          strand = ex$strand[j], stringsAsFactors = FALSE)
      }
    }
  }
  ann <- if (length(hits)) do.call(rbind, hits) else
    data.frame(.row = integer(0), transcript_id = character(0), gene_id = character(0),
               aa_position = integer(0), codon_base = integer(0), strand = character(0),
               stringsAsFactors = FALSE)

  if (keep == "hits") {
    if (!nrow(ann)) return(tibble::as_tibble(df[0, , drop = FALSE]))
    out <- cbind(df[ann$.row, , drop = FALSE], ann[, -1, drop = FALSE])
    out$coding <- TRUE
  } else {
    miss <- setdiff(seq_len(nrow(df)), ann$.row)
    # the NA filler cannot be recycled to zero rows, so only add it when something is missing
    if (length(miss))
      ann <- rbind(ann, data.frame(.row = miss, transcript_id = NA_character_,
                                   gene_id = NA_character_, aa_position = NA_integer_,
                                   codon_base = NA_integer_, strand = NA_character_,
                                   stringsAsFactors = FALSE))
    ann <- ann[order(ann$.row), , drop = FALSE]
    out <- cbind(df[ann$.row, , drop = FALSE], ann[, -1, drop = FALSE])
    out$coding <- !is.na(out$transcript_id)
  }
  rownames(out) <- NULL
  tibble::as_tibble(out)
}
