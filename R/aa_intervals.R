# Amino-acid positions -> the genomic interval of their codon.

# The standard genetic code. Written out rather than taken from Biostrings, which would put a
# Bioconductor package in the way of one lookup table.
.CODON_TABLE <- c(
  TTT = "F", TTC = "F", TTA = "L", TTG = "L", CTT = "L", CTC = "L", CTA = "L", CTG = "L",
  ATT = "I", ATC = "I", ATA = "I", ATG = "M", GTT = "V", GTC = "V", GTA = "V", GTG = "V",
  TCT = "S", TCC = "S", TCA = "S", TCG = "S", CCT = "P", CCC = "P", CCA = "P", CCG = "P",
  ACT = "T", ACC = "T", ACA = "T", ACG = "T", GCT = "A", GCC = "A", GCA = "A", GCG = "A",
  TAT = "Y", TAC = "Y", TAA = "*", TAG = "*", CAT = "H", CAC = "H", CAA = "Q", CAG = "Q",
  AAT = "N", AAC = "N", AAA = "K", AAG = "K", GAT = "D", GAC = "D", GAA = "E", GAG = "E",
  TGT = "C", TGC = "C", TGA = "*", TGG = "W", CGT = "R", CGC = "R", CGA = "R", CGG = "R",
  AGT = "S", AGC = "S", AGA = "R", AGG = "R", GGT = "G", GGC = "G", GGA = "G", GGG = "G")

# Complement single bases in place. The codon positions are already in transcript order, so a
# minus-strand codon needs complementing but not reversing.
.complement <- function(x) chartr("ACGTNacgtn", "TGCANtgcan", x)

# FASTA text -> named character vector of sequences, keyed by the first token of the header
# (">Pf3D7_01_v3 | organism=..." is one id plus a description).
.parse_fasta <- function(lines) {
  h <- which(startsWith(lines, ">"))
  if (!length(h)) return(NULL)
  ids <- sub("^>\\s*([^\\s|]+).*$", "\\1", lines[h], perl = TRUE)
  from <- h + 1L
  to <- c(h[-1] - 1L, length(lines))
  seqs <- vapply(seq_along(h), function(i) {
    if (from[i] > to[i]) "" else paste(lines[from[i]:to[i]], collapse = "")
  }, character(1))
  stats::setNames(toupper(gsub("[[:space:]]", "", seqs)), ids)
}

# Read a FASTA from a path or URL (plain or gzipped), or pass through a named vector that is
# already sequence.
.read_fasta <- function(fasta) {
  if (is.null(fasta)) return(NULL)
  # names are what separates "these are sequences" from "this is where to read them from" --
  # a one-sequence vector is still a vector of sequences
  if (is.character(fasta) && !is.null(names(fasta)))
    return(stats::setNames(toupper(gsub("[[:space:]]", "", fasta)), names(fasta)))
  if (!is.character(fasta) || length(fasta) != 1)
    stop("`fasta` must be one path or URL, or a named vector of sequences", call. = FALSE)
  if (grepl("^(https?|ftp)://", fasta)) {
    con <- if (grepl("\\.gz$", fasta)) gzcon(url(fasta, open = "rb")) else url(fasta)
    on.exit(try(close(con), silent = TRUE), add = TRUE)
    lines <- readLines(con, warn = FALSE)
  } else {
    if (!file.exists(fasta)) stop("no such file: ", fasta, call. = FALSE)
    lines <- readLines(fasta, warn = FALSE)
  }
  out <- .parse_fasta(lines)
  if (is.null(out)) stop("no FASTA records in ", fasta, call. = FALSE)
  out
}

# Sequence for the reference bases: an explicit `fasta` wins, otherwise whatever came embedded
# in the GFF, otherwise nothing and the reference columns are left off.
.resolve_sequence <- function(fasta, cds) {
  if (!is.null(fasta)) return(.read_fasta(fasta))
  attr(cds, "sequence", exact = TRUE)
}

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
#' # the CDS of six 3D7 drug-resistance genes ship with the package
#' cds <- read_gff_cds(system.file("extdata", "pf3d7_drug_gene_cds.gff",
#'                                 package = "plasgenomicsutilsR"))
#' cds
#'
#' \dontrun{
#' # a whole released annotation, read straight from the web
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
  out <- tibble::as_tibble(out[order(out$transcript_id, out$start), , drop = FALSE])
  # carried along so snp_aa_positions() can report reference codons without being handed the
  # genome a second time
  if (!is.null(seqs)) {
    attr(out, "sequence") <- seqs
    message(length(seqs), " sequence(s) read from the GFF's own ##FASTA section")
  }
  out
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
#' cds <- read_gff_cds(system.file("extdata", "pf3d7_drug_gene_cds.gff",
#'                                 package = "plasgenomicsutilsR"))
#'
#' # K76T and C580Y, at the positions they are reported at
#' aa <- data.frame(transcript_id = c("pfcrt", "pfkelch13"), aa_position = c(76, 580))
#' aa_intervals(aa, cds, one_based_output = TRUE)[, c("name", "chr", "start", "end", "strand")]
#'
#' # pfcrt has 13 coding exons, so four of its codons straddle an intron
#' every_codon <- aa_intervals(data.frame(transcript_id = "pfcrt", aa_position = 1:424), cds)
#' subset(every_codon, spans_intron)[, c("name", "codon_positions")]
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

# Genomic position of 1-based transcript coordinates (measured along the concatenated CDS from
# the first exon's start, phase included). Vectorised, and NA past either end of the CDS.
.tx_to_genomic <- function(o, t_full) {
  ex <- o$ex
  len <- ex$end - ex$start + 1L
  ends <- cumsum(len)
  starts <- ends - len + 1L
  out <- rep(NA_real_, length(t_full))
  ok <- !is.na(t_full) & t_full >= 1 & t_full <= ends[length(ends)]
  if (!any(ok)) return(out)
  j <- findInterval(t_full[ok], starts)
  off <- t_full[ok] - starts[j]
  out[ok] <- ifelse(ex$strand[j] == "-", ex$end[j] - off, ex$start[j] + off)
  out
}

# The reference codon and residue for rows whose codon positions are known. `seqs` is a named
# vector of chromosome sequences; anything it cannot cover comes back NA rather than guessed.
.ref_codons <- function(chrom, p1, p2, p3, strand, seqs) {
  n <- length(chrom)
  codon <- rep(NA_character_, n)
  key <- normalise_chr(names(seqs))
  want <- normalise_chr(chrom)
  base_at <- function(ch, p) {
    b <- rep(NA_character_, length(p))
    for (u in unique(ch[!is.na(ch)])) {
      s <- seqs[[which(key == u)[1]]]
      i <- which(ch == u & !is.na(p))
      if (!length(i) || is.null(s)) next
      inside <- p[i] >= 1 & p[i] <= nchar(s)
      b[i[inside]] <- substring(s, p[i][inside], p[i][inside])
    }
    b
  }
  have <- want %in% key & !is.na(p1) & !is.na(p2) & !is.na(p3)
  if (!any(have)) return(list(codon = codon, aa = rep(NA_character_, n)))
  i <- which(have)
  b <- vapply(list(p1, p2, p3), function(p) base_at(want[i], p[i]), character(length(i)))
  if (length(i) == 1L) b <- matrix(b, nrow = 1L)
  # positions are in transcript order already, so a minus-strand codon is complemented in place
  neg <- !is.na(strand[i]) & strand[i] == "-"
  b[neg, ] <- .complement(b[neg, , drop = FALSE])
  cod <- toupper(paste0(b[, 1], b[, 2], b[, 3]))
  cod[grepl("NA", cod, fixed = TRUE)] <- NA_character_
  codon[i] <- cod
  list(codon = codon, aa = unname(.CODON_TABLE[codon]))
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
#' alleles, which are not part of this. What reference sequence does buy you is the residue the
#' codon currently codes for -- see below.
#'
#' @section The reference residue:
#' Given sequence, the result also carries `ref_codon` (the codon's three bases, in transcript
#' orientation, complemented on the minus strand) and `ref_aa` (its residue, one letter, `*` for
#' a stop). Sequence can come from either place:
#'
#' * `fasta =` -- a path or URL to a genome FASTA (plain or gzipped), or a named vector of
#'   sequences you already have in hand.
#' * the GFF itself, when it ends with a `##FASTA` section. [read_gff_cds()] keeps those
#'   sequences, so nothing extra needs passing.
#'
#' Without either, both columns are simply absent: the reference base cannot be inferred from an
#' annotation alone. Some GFFs do carry a translated protein in a feature attribute, but too few
#' agree on how to make it worth reading, so it is not used.
#'
#' Names are matched through [normalise_chr()], so a FASTA headed `Pf3D7_07_v3` lines up with a
#' GFF spelling it the same way or differently. If nothing matches, both columns come back `NA`
#' with a warning naming what was compared -- silence there would look like a genome with no
#' coding SNPs.
#'
#' This is worth doing even when you think you know the answer. `ref_aa` is what the *reference*
#' carries, and a reference is one isolate's genome -- not a consensus, and not the ancestral or
#' wild-type sequence. Pf3D7 reads `G` at *pfdhps* 437, the residue A437G is named for changing
#' *to*, so a callset aligned to 3D7 shows no variant at that position precisely because the
#' reference is already the mutant. This is not a quirk of one locus: across *Plasmodium* species
#' the reference is sometimes the non-wild-type allele, so "REF" means "what this isolate has",
#' never "what came first".
#'
#' Two conventions here follow from that. `load_genotypes()` records which allele its dosages
#' count rather than letting REF stand in for a baseline, and [run_ihs()] and [plot_ehh()] default
#' to `polarized = FALSE`, treating the two states as simply the two states, since ancestral
#' versus derived cannot be read off REF and ALT.
#'
#' @param snps A data frame with `snp_id` (`"chr:pos"`) or `chr` and `pos` columns; any other
#'   columns are carried through.
#' @param gff A GFF path, or the result of [read_gff_cds()] (parse once, reuse).
#' @param keep `"all"` (default) keeps non-coding SNPs with `NA` annotations; `"hits"` keeps
#'   only those in a CDS.
#' @param one_based_snps The positions in `snps` are 1-based (VCF `POS`). Genotype-matrix column
#'   names from [load_genotypes()] are, so set this when feeding those in; the package's own
#'   tables are 0-based, which is the default.
#' @param fasta Optional reference sequence, for the `ref_codon` / `ref_aa` columns: a path or
#'   URL to a FASTA, or a named character vector of sequences. Defaults to whatever `gff`
#'   carried; `NULL` with a GFF holding no sequence leaves the two columns off.
#' @return `snps` with `transcript_id`, `gene_id`, `aa_position` (1-based), `codon_base`
#'   (1/2/3 in transcript orientation), `strand` and `coding` added, plus `ref_codon` and
#'   `ref_aa` when there is sequence to read them from.
#' @seealso [aa_intervals()] for the other direction, [annotate_snps()] to first ask which
#'   gene a SNP is in.
#' @examples
#' cds <- read_gff_cds(system.file("extdata", "pf3d7_drug_gene_cds.gff",
#'                                 package = "plasgenomicsutilsR"))
#'
#' # the three bases of pfkelch13 codon 580, on the minus strand: base 1 is the highest
#' snp_aa_positions(data.frame(chr = "Pf3D7_13_v3", pos = c(1725260, 1725259, 1725258)),
#'                  cds, keep = "hits", one_based_snps = TRUE)[
#'   , c("pos", "aa_position", "codon_base", "strand")]
#'
#' # `ref_codon` / `ref_aa` need sequence. A GFF ending in a `##FASTA` section carries its own,
#' # so nothing extra is passed:
#' gff <- tempfile(fileext = ".gff")
#' writeLines(c("##gff-version 3",
#'              "demo\t.\tCDS\t11\t25\t.\t+\t0\tID=c1;Parent=T.1;gene_id=T",
#'              "##FASTA", ">demo", "CCCCCCCCCCATGAAATTTGGGTAAC"), gff)
#' snp_aa_positions(data.frame(chr = "demo", pos = 11:22), read_gff_cds(gff),
#'                  keep = "hits", one_based_snps = TRUE)[
#'   , c("pos", "aa_position", "codon_base", "ref_codon", "ref_aa")]
#'
#' \dontrun{
#' # every genotyped SNP that is coding, and the residue it sits on. Genotype-matrix ids
#' # are 0-based, so no `one_based_snps` here.
#' snp_aa_positions(data.frame(snp_id = colnames(ps$genotype("full"))), cds, keep = "hits")
#'
#' # or point `fasta` at the released genome, read straight from the web like the GFF is.
#' # pfcrt codon 76 comes back "AAA" / "K".
#' genome <- paste0("https://plasmodb.org/common/downloads/Current_Release/",
#'                  "Pfalciparum3D7/fasta/data/PlasmoDB-68_Pfalciparum3D7_Genome.fasta")
#' snp_aa_positions(data.frame(chr = "Pf3D7_07_v3", pos = 403625), cds, keep = "hits",
#'                  one_based_snps = TRUE, fasta = genome)
#' }
#' @export
snp_aa_positions <- function(snps, gff, keep = c("all", "hits"), one_based_snps = FALSE,
                             fasta = NULL) {
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

  parsed <- if (is.data.frame(gff)) gff else read_gff_cds(gff)
  seqs <- .resolve_sequence(fasta, parsed)
  cds <- as.data.frame(parsed, stringsAsFactors = FALSE)
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
        # where the codon's three bases sit, so a reference codon can be read off the genome.
        # Back to full transcript coordinates (phase included) to reuse the exon walk.
        t1 <- (aa - 1) * 3 + 1 + o$skip
        hits[[length(hits) + 1L]] <- data.frame(
          .row = k, transcript_id = tx, gene_id = ex$gene_id[j],
          aa_position = as.integer(aa), codon_base = as.integer(base),
          strand = ex$strand[j], .chrom = ex$chrom[j],
          .p1 = .tx_to_genomic(o, t1), .p2 = .tx_to_genomic(o, t1 + 1),
          .p3 = .tx_to_genomic(o, t1 + 2), stringsAsFactors = FALSE)
      }
    }
  }
  ann <- if (length(hits)) do.call(rbind, hits) else
    data.frame(.row = integer(0), transcript_id = character(0), gene_id = character(0),
               aa_position = integer(0), codon_base = integer(0), strand = character(0),
               .chrom = character(0), .p1 = numeric(0), .p2 = numeric(0), .p3 = numeric(0),
               stringsAsFactors = FALSE)

  if (keep == "hits") {
    if (!nrow(ann)) return(tibble::as_tibble(df[0, , drop = FALSE]))
    # in the order the SNPs were given, not the order the chromosomes happened to be walked --
    # the rows line up with `snps` so columns can be carried alongside the result
    ann <- ann[order(ann$.row), , drop = FALSE]
    out <- cbind(df[ann$.row, , drop = FALSE], ann[, -1, drop = FALSE])
    out$coding <- TRUE
  } else {
    miss <- setdiff(seq_len(nrow(df)), ann$.row)
    # the NA filler cannot be recycled to zero rows, so only add it when something is missing
    if (length(miss))
      ann <- rbind(ann, data.frame(.row = miss, transcript_id = NA_character_,
                                   gene_id = NA_character_, aa_position = NA_integer_,
                                   codon_base = NA_integer_, strand = NA_character_,
                                   .chrom = NA_character_, .p1 = NA_real_, .p2 = NA_real_,
                                   .p3 = NA_real_, stringsAsFactors = FALSE))
    ann <- ann[order(ann$.row), , drop = FALSE]
    out <- cbind(df[ann$.row, , drop = FALSE], ann[, -1, drop = FALSE])
    out$coding <- !is.na(out$transcript_id)
  }
  # reference codon and residue, only when there is sequence to read them off
  if (!is.null(seqs) && nrow(out)) {
    rc <- .ref_codons(out$.chrom, out$.p1, out$.p2, out$.p3, out$strand, seqs)
    out$ref_codon <- rc$codon
    out$ref_aa <- rc$aa
    if (all(is.na(rc$codon[out$coding])) && any(out$coding))
      warning("no reference codon could be read: the sequence names (",
              paste(utils::head(names(seqs), 3), collapse = ", "),
              ") do not match the GFF's chromosomes (",
              paste(utils::head(unique(stats::na.omit(out$.chrom)), 3), collapse = ", "), ")",
              call. = FALSE)
  }
  out <- out[, setdiff(names(out), c(".chrom", ".p1", ".p2", ".p3")), drop = FALSE]
  rownames(out) <- NULL
  tibble::as_tibble(out)
}
