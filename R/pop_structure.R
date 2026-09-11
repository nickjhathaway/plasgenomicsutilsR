# Population-structure analysis: LD-pruned genotypes -> PCA / UMAP, and sNMF
# admixture, with plot_*() functions. Heavy analysis packages (SNPRelate, gdsfmt,
# uwot, LEA) are optional (Suggests), guarded at call time.

# ---- genotypes -------------------------------------------------------------

# The GDS builder parses a VCF as text, so a binary BCF reaches it as mojibake and fails deep
# inside the parser with "invalid multibyte string". Hand it a VCF instead: reuse one
# already sitting next to the BCF if there is one, otherwise convert with bcftools into
# `vcf_dir` (default: alongside the BCF). Reusing rather than re-converting is the point
# -- these files are large and otherwise accumulate one copy per analysis.
.as_text_vcf <- function(path, vcf_dir = NULL, refresh = "stale") {
  if (!grepl("\\.bcf$", path, ignore.case = TRUE)) return(path)
  base <- sub("\\.bcf$", "", path, ignore.case = TRUE)
  out <- if (is.null(vcf_dir)) paste0(base, ".vcf.gz")
         else file.path(vcf_dir, paste0(basename(base), ".vcf.gz"))

  found <- NULL
  for (cand in unique(c(out, paste0(base, ".vcf.gz"), paste0(base, ".vcf")))) {
    if (file.exists(cand)) { found <- cand; break }
  }
  if (!is.null(found) && !identical(refresh, "always")) {
    if (file.mtime(found) >= file.mtime(path)) {
      message("reusing the existing VCF: ", found)
      return(found)
    }
    # The BCF has moved on. Reusing the VCF here is not a stale cache in the harmless
    # sense: everything downstream -- the GDS, the pruned panel, every analysis -- is then
    # built from the old records, and the GDS looks fresh because it is newer than the
    # stale VCF it came from. So the default is to reconvert, and reusing it anyway has to
    # be asked for.
    if (identical(refresh, "never")) {
      warning(sprintf(paste("%s is older than the BCF beside it, and refresh = \"never\":",
                            "every result below is built from the older records."),
                      basename(found)), call. = FALSE)
      message("reusing the existing VCF: ", found)
      return(found)
    }
    message("the VCF beside the BCF is older than it; reconverting")
  }

  bcftools <- Sys.which("bcftools")
  if (!nzchar(bcftools))
    stop(sprintf(paste0("`%s` is a BCF, which SNPRelate cannot read, and bcftools is ",
                        "not on PATH to convert it.\n  bcftools view -Oz -o %s %s"),
                 basename(path), out, path), call. = FALSE)
  if (!is.null(vcf_dir)) dir.create(vcf_dir, showWarnings = FALSE, recursive = TRUE)
  message("converting BCF to ", out, " (the GDS builder reads VCF text only)")
  log <- system2(bcftools, c("view", "-Oz", "-o", shQuote(out), shQuote(path)),
                 stdout = TRUE, stderr = TRUE)
  st <- attr(log, "status")
  if ((!is.null(st) && st != 0L) || !file.exists(out))
    stop("bcftools failed to convert ", path, ":\n  ",
         paste(utils::tail(log, 3), collapse = "\n  "), call. = FALSE)
  out
}

# Which `variants` a cached GDS was built with, recorded as a node in the file itself so a
# reused GDS is never mistaken for one built the other way. A GDS from before this was written
# has no node and was necessarily biallelic-only.
.GDS_VARIANTS_NODE <- "plasgenomicsutilsR.variants"

.gds_variants <- function(gds) {
  if (!file.exists(gds)) return(NA_character_)
  f <- try(gdsfmt::openfn.gds(gds), silent = TRUE)
  if (inherits(f, "try-error")) return(NA_character_)
  on.exit(gdsfmt::closefn.gds(f), add = TRUE)
  n <- try(gdsfmt::index.gdsn(f, .GDS_VARIANTS_NODE), silent = TRUE)
  if (inherits(n, "try-error")) "biallelic_snvs" else as.character(gdsfmt::read.gdsn(n))[1]
}

.tag_gds_variants <- function(gds, variants) {
  f <- try(gdsfmt::openfn.gds(gds, readonly = FALSE), silent = TRUE)
  if (inherits(f, "try-error")) return(invisible(FALSE))
  on.exit(gdsfmt::closefn.gds(f), add = TRUE)
  old <- try(gdsfmt::index.gdsn(f, .GDS_VARIANTS_NODE), silent = TRUE)
  if (!inherits(old, "try-error")) gdsfmt::delete.gdsn(old)
  gdsfmt::add.gdsn(f, .GDS_VARIANTS_NODE, variants, closezip = TRUE)
  invisible(TRUE)
}
# SNPRelate reports a record's alleles as "REF/ALT1,ALT2". Splitting that is the whole of
# the allele-identity recovery: the strings are already in the GDS and were simply not read.
.parse_alleles <- function(x) {
  x <- as.character(x)
  ref <- ifelse(is.na(x), NA_character_, sub("/.*$", "", x))
  rest <- ifelse(is.na(x) | !grepl("/", x), "", sub("^[^/]*/", "", x))
  alt <- lapply(rest, function(z) {
    if (!nzchar(z)) return(character(0))
    a <- strsplit(z, ",", fixed = TRUE)[[1]]
    a[nzchar(a) & a != "."]
  })
  list(ref = ref, alt = alt, n_alt = vapply(alt, length, integer(1)))
}

# "1 multiallelic, 2 indel" -- only the kinds actually present, so a clean callset says so
# rather than listing three zeroes.
#
# The tally is exact: SeqArray has already read every record, so the classification the
# selection used is the classification reported. Under the SNPRelate backend this had to be
# taken by shelling out to bcftools, and went unstated when bcftools was absent.
.skipped_note <- function(counts) {
  if (is.null(counts) || !length(counts)) return("every record in the callset was read")
  paste0("not read: ",
         paste(sprintf("%d %s", unlist(counts), names(counts)), collapse = ", "))
}

# ---- SeqArray backend ------------------------------------------------------

# A SeqVarGDS holds allele *indices*, so it can carry a record with any number of alleles.
# A SNP-GDS holds a dosage -- how many copies of one allele a sample has -- which is two
# numbers where a multiallelic site needs k. Measured on 200 real multiallelic Pf7 records,
# SNPRelate's `biallelic.only` read **none** of them and `copy.num.of.ref` read all 200 and
# gave every alternate the same number. That is the whole reason the backend changed.
#
# One consequence worth knowing: the GDS no longer depends on which records you asked for.
# SeqArray reads everything and the selection happens in R, so one file serves every
# `variants`/`encoding` combination and the old per-`variants` cache tag is gone.
.GDS_BACKEND <- "seqarray"

# Which class each record belongs to, from its allele strings. Same reading as the Python
# package's `classify_record`, so the two agree about what a panel will hold: a record is one
# class, and a `*` disqualifies it -- part of the cohort has no base there to compare, so it
# is not a clean SNP site whatever its other alleles read.
.classify_alleles <- function(ref, alt) {
  vapply(seq_along(ref), function(i) {
    a <- alt[[i]]
    a <- a[nzchar(a) & a != "."]
    if (!length(a)) return("no_alt")
    if (any(a == "*")) return("spanning_del")
    kind <- vapply(a, function(x) {
      if (nchar(x) != nchar(ref[i])) return("indel")
      if (nchar(ref[i]) == 1L) return("snv")
      # Equal length and more than one base is not automatically an MNP: a single
      # substitution is often written with padding -- REF=TTATA ALT=CTATA differs only at
      # the first base -- and counting those as MNPs invents a population of them that is
      # not there. bcftools reads them as SNVs; so does the companion Python package's
      # `classify_record`, and so must this or a panel and its callset disagree.
      rr <- strsplit(ref[i], "")[[1]]; xx <- strsplit(x, "")[[1]]
      if (sum(rr != xx) == 1L) "snv" else "mnp"
    }, character(1))
    if (!all(kind == "snv")) return(if (length(unique(kind)) > 1L) "mixed" else
                                    if (kind[1] == "mnp") "mnp" else "indel")
    if (length(a) > 1L) return("multiallelic")
    "biallelic_snv"
  }, character(1))
}

# ploidy x sample x variant allele indices -> what the caller asked for.
#
# `dosage` counts alternate copies, which only means something when there is one alternate;
# `allele_index` names which allele a sample carries. A heterozygous call is a mixed
# infection rather than a diploid genotype -- the package's standing reading -- and a single
# index cannot name two clones, so it becomes missing under `allele_index`.
.encode_genotypes <- function(gt, encoding) {
  d <- dim(gt)
  ploidy <- d[1]
  out <- matrix(NA_integer_, d[2], d[3])
  if (identical(encoding, "dosage")) {
    alt <- apply(gt > 0L, c(2, 3), sum)
    called <- apply(!is.na(gt), c(2, 3), all)
    out[called] <- as.integer(alt[called]) * (if (ploidy == 1L) 2L else 1L)
    return(out)
  }
  first <- gt[1, , , drop = FALSE]
  dim(first) <- d[-1]
  if (ploidy == 1L) return(matrix(as.integer(first), d[2], d[3]))
  same <- apply(gt, c(2, 3), function(v) !anyNA(v) && length(unique(v)) == 1L)
  out[same] <- as.integer(first[same])
  out
}

# A sites table is only useful if it lines up with the matrix it describes, so a panel that
# holds a different set of columns gets NULL rather than a table that quietly disagrees.
.align_sites <- function(sites, cols) {
  if (is.null(sites) || !is.data.frame(sites) || is.null(sites$site_key)) return(NULL)
  k <- match(cols, sites$site_key)
  if (anyNA(k)) return(NULL)
  out <- sites[k, , drop = FALSE]
  rownames(out) <- NULL
  out
}

# One row per record kept, in the matrix's column order. Built from what SNPRelate already
# returns, so it costs no extra I/O and no new dependency.
.sites_table <- function(snpinfo, idx, site_key) {
  al <- .parse_alleles(if (is.null(snpinfo$allele)) NA_character_ else snpinfo$allele[idx])
  out <- data.frame(
    site_key = site_key,
    chr = as.character(snpinfo$chromosome[idx]),
    pos = as.numeric(snpinfo$position[idx]) - 1,   # 0-based, as everywhere else
    ref = al$ref,
    n_alt = al$n_alt,
    stringsAsFactors = FALSE
  )
  # assigned rather than passed to data.frame(), which would wrap it in AsIs and make every
  # downstream comparison against a plain list fail on class alone
  out$alt <- al$alt
  out[c("site_key", "chr", "pos", "ref", "alt", "n_alt")]
}

#' Load genotypes from a VCF, optionally LD-pruned
#'
#' Converts a VCF to GDS (only when needed) and returns the genotype matrix, LD-pruned by
#' default. The backend is \pkg{SeqArray}, whose SeqVarGDS stores allele **indices** and so
#' can hold a record with any number of alleles.
#'
#' Pruned or not depends on the question. Pruning is right for PCA, UMAP and admixture, where
#' correlated SNPs would let one locus dominate the structure. It is wrong wherever the
#' correlation between neighbouring SNPs *is* the signal -- differentiation ([pop_diff()])
#' and haplotypes ([plot_region_haplotypes()]) -- because it keeps one SNP out of each
#' correlated run and drops the rest. Holding both is cheap: the GDS is reused, so a second
#' call with `prune = FALSE` only re-reads it.
#'
#' @section Which records reach the panel:
#' `prune = FALSE` means unpruned, not every record: the default `variants = "biallelic_snvs"`
#' keeps only records whose single ALT is a substitution, so the panel is usually smaller than
#' the VCF's record count. Three kinds are skipped, and the message says **how many of each**:
#'
#' * sites with **no ALT allele** (`ALT="."`) -- the reference positions an all-sites caller
#'   emits. Usually the biggest share by far, and the easiest to miss, since nothing about
#'   them says "variant": `bcftools view --exclude-types indels` leaves every one in place,
#'   because they are not indels.
#' * **indels** and other non-SNV records.
#' * sites with **more than one ALT** -- which `variants = "all"` keeps.
#'
#' The tally is exact and costs nothing: SeqArray has already read every record by the time
#' the selection happens, so the classification that chose the panel is the one reported. A
#' record carrying a `*` is its own class and is skipped: part of the cohort has no base
#' there to compare, so `A > *,T` is not a clean SNP site however well T behaves. The
#' companion Python package's `spanning_del_filter` is the way to keep such a site -- it
#' recodes the deleted calls as missing and drops the allele, after which the record really
#' is an ordinary SNP.
#'
#' Sites that are **invariant across the loaded samples are kept**, as long as the record
#' lists an ALT. A site every sample calls `1/1`, or every sample calls `0/0`, comes through
#' -- which is why the functions needing variable sites ([parasite_haplotypes()],
#' [run_ihs()]) apply their own `maf` cutoff instead of trusting the panel.
#'
#' @section Keeping every variant:
#' `variants = "all"` keeps every record -- multiallelic sites, indels and `ALT="."` positions
#' included. Pair it with `encoding = "allele_index"` to keep the alternates apart; with the
#' default `encoding = "dosage"` it warns, because a dosage cannot say which alternate a
#' non-reference call carries:
#'
#' ```
#' encoding = "dosage"                  encoding = "allele_index"
#'   pos  alleles  s1 s2 s3               pos  alleles  s1 s2 s3
#'   400  C/T,G     0  2  2               400  C/T,G     0  1  2
#'                     ^  ^                                  ^  ^
#'                     1/1 and 2/2                           T and G, told apart
#'                     land on one number
#' ```
#'
#' One GDS serves every combination: SeqArray reads the whole callset and the selection
#' happens in R, so asking for a different `variants` or `encoding` re-reads the same file
#' rather than rebuilding it.
#'
#' @param vcf Path to a (bgzipped) VCF, or a **BCF** -- converted to VCF text first with
#'   `bcftools`, reusing any VCF already sitting next to it rather than making another copy.
#' @param gds Optional GDS path; derived from `vcf` if `NULL`.
#' @param refresh What to do with a derived file older than what it was built from -- the
#'   text VCF beside a BCF, and the GDS beside either. `"stale"` (default) rebuilds it,
#'   `"always"` rebuilds regardless, `"never"` reuses it and warns. The default is a
#'   rebuild because the alternative is silent: update the BCF and every result below comes
#'   from the old records, with the GDS looking current because it is newer than the stale
#'   VCF it was built from. A GDS left over from the older \pkg{SNPRelate} backend is a
#'   different format entirely and is rebuilt whatever `refresh` says.
#' @param prune LD-prune. Defaults to `TRUE` for `encoding = "dosage"` and `FALSE` otherwise,
#'   because pruning measures correlation between dosages and an allele index is not a count
#'   of anything. `FALSE` returns every record `variants` admits
#'   (see *Which records reach the panel*), unpruned -- use this for the genotype matrix fed
#'   to [pop_diff()] / [pop_diff_table()], since LD-pruning removes the very SNPs that carry
#'   the differentiation signal.
#' @param ld_threshold,slide_max_bp,slide_max_n,autosome_only Passed to
#'   [SNPRelate::snpgdsLDpruning()] (defaults 0.2 / 20000 / 200 / `FALSE`); ignored
#'   when `prune = FALSE`.
#' @param maf,missing_rate Optional MAF / per-SNP missing-rate cutoffs for pruning.
#' @param seed Random seed for the pruning.
#' @param vcf_dir Where to put the VCF converted from a BCF (default: alongside the BCF).
#'   Point it somewhere scratch to keep converted copies out of the data directory.
#' @param allele Which allele the returned dosage counts, `"alt"` (default) or `"ref"`. Only
#'   reported allele frequencies (and the arbitrary sign of a PCA axis) depend on this --
#'   every diversity, differentiation, LD and selection statistic here is symmetric in `p`
#'   and `1 - p`. Meaningless for `encoding = "allele_index"`, which names an allele rather
#'   than counting copies of one, and refused there.
#' @param variants Which records to read. `"biallelic_snvs"` (default) keeps biallelic SNVs
#'   only; `"all"` keeps every record, multiallelic sites and indels included.
#' @param encoding How `genotype` is coded. `"dosage"` (default) counts alternate copies,
#'   0/1/2 -- the classical matrix, and what PCA, admixture, diversity and differentiation
#'   read. `"allele_index"` names **which** allele each sample carries, 0 for the reference
#'   and 1..k for the alternates, which is the only one of the two that can carry a
#'   multiallelic site. See *Which encoding to ask for*.
#' @section Which encoding to ask for:
#' A **dosage** says how many copies of one allele a sample has. That is two numbers, and a
#' site with three alleles needs three, so a dosage cannot say *which* alternate a
#' non-reference call carries -- at a codon where D384A, D384G and D384Y arose separately,
#' every carrier lands on the same number. Ask for it when the panel is biallelic and the
#' consumer is PCA, UMAP, admixture, diversity or differentiation.
#'
#' An **allele index** names the allele: 0 for the reference, 1..k for the alternates in the
#' record's own ALT order (`$sites$alt` says which is which). It carries a multiallelic site
#' faithfully, and the dosage statistics refuse it rather than mangling it -- there is no
#' faithful conversion, so contrast one allele at a time instead.
#'
#' The parasite is haploid, so a heterozygous call is a **mixed infection** rather than a
#' diploid genotype, and one index cannot name two clones. Those calls are missing under
#' `allele_index`, which is the same reading `pop_diversity(het = "missing")` takes.
#'
#' @return A list with `genotype` (matrix; sample row names and `chr:pos0` column names --
#'   0-based, like every other position in the package),
#'   `sample.id`, `snp.id`, and the facts the matrix itself cannot carry: `allele` (which
#'   allele a dosage counts), `encoding`, `pruned`, `positions` and `variants`.
#'   [PopStructure] keeps them, so anything that names a call, refuses a non-dosage panel or
#'   warns about pruning can ask instead of assuming.
#'
#'   Also `sites`: one row per genotype column, in the same order, with `site_key`, `chr`,
#'   `pos` (0-based), `ref`, `alt`, `n_alt`, `n_alt_real` and `has_spanning_del`.
#'
#'   `alt` is a list column holding the record's **own** ALT list, in its own order and `*`
#'   included, because that is what an allele index names -- `alt[[i]][k]` is the allele a
#'   genotype of `k` refers to, and dropping anything from it would point the index at the
#'   wrong base. `n_alt` counts that list for the same reason.
#'
#'   Whether a record is *multiallelic* is a different question, since `*` is a missingness
#'   annotation rather than an allele: `n_alt_real` excludes it and is the one that answers,
#'   with `has_spanning_del` recording the `*` separately. So `A > *,T` is an ordinary
#'   biallelic SNV with `n_alt = 2` and `n_alt_real = 1`.
#'
#'   **Ask this table before reading a per-allele number off a dosage matrix**: a column with
#'   `n_alt_real > 1` has had its alternates collapsed onto one number and cannot answer one.
#' @seealso [PopStructure], [pop_structure()]
#' @examples
#' \dontrun{
#' # LD-pruned for PCA / UMAP / admixture
#' geno <- load_genotypes("clean.vcf.gz", gds = "clean.gds")
#'
#' # and the full panel for anything where SNP correlation IS the signal
#' full <- load_genotypes("clean.vcf.gz", gds = "clean.gds", prune = FALSE)
#' }
#' @export
load_genotypes <- function(vcf, gds = NULL, prune = NULL, ld_threshold = 0.2,
                         slide_max_bp = 20000, slide_max_n = 200, autosome_only = FALSE,
                         maf = NaN, missing_rate = NaN, seed = 42, vcf_dir = NULL,
                         allele = c("alt", "ref"),
                         variants = c("biallelic_snvs", "all"),
                         encoding = c("dosage", "allele_index"),
                         star = c("missing", "allele"),
                         refresh = c("stale", "never", "always")) {
  .need_package("SeqArray", "load_genotypes()")
  .need_package("gdsfmt", "load_genotypes()")
  allele <- match.arg(allele)
  variants <- match.arg(variants)
  encoding <- match.arg(encoding)
  star <- match.arg(star)
  refresh <- match.arg(refresh)
  # Pruning removes SNPs for correlating with a neighbour, which is what PCA and admixture
  # want and what they need dosages for. An allele-index panel is not headed there, so the
  # default follows the encoding rather than forcing the caller to switch it off.
  if (is.null(prune)) prune <- identical(encoding, "dosage")
  if (isTRUE(prune) && !identical(encoding, "dosage"))
    stop("`prune = TRUE` needs `encoding = \"dosage\"`: LD pruning measures correlation ",
         "between dosages, and an allele index is not a count of anything.", call. = FALSE)
  if (identical(allele, "ref") && !identical(encoding, "dosage"))
    stop("`allele = \"ref\"` needs `encoding = \"dosage\"`: an allele index names which ",
         "allele a sample carries rather than counting copies of one, so there is nothing ",
         "to flip.", call. = FALSE)
  if (!file.exists(vcf)) stop(sprintf("no such file: %s", vcf), call. = FALSE)
  vcf <- .as_text_vcf(vcf, vcf_dir, refresh)
  if (is.null(gds)) gds <- sub("\\.vcf(\\.gz)?$", ".gds", vcf, ignore.case = TRUE)
  if (identical(gds, vcf)) gds <- paste0(vcf, ".gds")

  # SeqArray reads every record, so one GDS serves every `variants`/`encoding` combination
  # and the selection happens below in R. A GDS left over from the SNPRelate backend holds a
  # different format entirely, so the tag also catches those and rebuilds.
  stale <- !file.exists(gds) || !identical(.gds_variants(gds), .GDS_BACKEND) ||
    identical(refresh, "always") ||
    (file.mtime(gds) < file.mtime(vcf) && !identical(refresh, "never"))
  if (file.exists(gds) && !stale && file.mtime(gds) < file.mtime(vcf))
    warning(sprintf(paste("%s is older than %s, and refresh = \"never\": it is being reused",
                          "rather than rebuilt."), basename(gds), basename(vcf)),
            call. = FALSE)
  if (stale) {
    SeqArray::seqVCF2GDS(vcf, gds, storage.option = "ZIP_RA", verbose = FALSE)
    .tag_gds_variants(gds, .GDS_BACKEND)
  }

  f <- SeqArray::seqOpen(gds)
  on.exit(SeqArray::seqClose(f), add = TRUE)
  # a filter set by an earlier call persists in the file, so start from everything
  SeqArray::seqResetFilter(f, verbose = FALSE)

  chrom <- as.character(SeqArray::seqGetData(f, "chromosome"))
  pos1 <- as.integer(SeqArray::seqGetData(f, "position"))
  al <- strsplit(SeqArray::seqGetData(f, "allele"), ",", fixed = TRUE)
  ref <- vapply(al, function(x) x[1], character(1))
  alt <- lapply(al, function(x) {
    a <- x[-1]
    a[nzchar(a) & a != "."]
  })
  kind <- .classify_alleles(ref, alt)
  keep <- if (identical(variants, "all")) rep(TRUE, length(kind))
          else kind == "biallelic_snv"
  if (!any(keep))
    stop("no records left after `variants = \"", variants, "\"` in ", basename(vcf),
         ": the callset holds ", paste(sprintf("%d %s", table(kind), names(table(kind))),
                                       collapse = ", "), ".", call. = FALSE)

  # Counted from the callset itself rather than asked of an external tool: SeqArray has
  # already read every record, so the tally is exact and costs nothing.
  skipped <- table(kind[!keep])
  skipped <- as.list(skipped[skipped > 0L])
  if (identical(variants, "all")) {
    message(sum(keep), " records in ", basename(vcf),
            " (every variant, including indels and multiallelic sites)")
  } else {
    message(sum(keep), " biallelic SNVs in ", basename(vcf),
            " (", .skipped_note(skipped), ")")
  }

  SeqArray::seqSetFilter(f, variant.sel = which(keep), verbose = FALSE)
  if (isTRUE(prune)) {
    .need_package("SNPRelate", "load_genotypes(prune = TRUE)")
    set.seed(seed)
    snpset <- SNPRelate::snpgdsLDpruning(
      f, autosome.only = autosome_only, ld.threshold = ld_threshold,
      slide.max.bp = slide_max_bp, slide.max.n = slide_max_n,
      maf = maf, missing.rate = missing_rate, verbose = FALSE)
    sel <- sort(unlist(snpset, use.names = FALSE))
    SeqArray::seqSetFilter(f, variant.id = sel, verbose = FALSE)
  }

  vid <- SeqArray::seqGetData(f, "variant.id")
  idx <- match(vid, seq_along(kind))
  samples <- as.character(SeqArray::seqGetData(f, "sample.id"))
  gt <- SeqArray::seqGetData(f, "genotype")
  mat <- .encode_genotypes(gt, encoding)
  # SeqArray reports 1-based VCF POS. Every SNP id and interval in both packages is 0-based
  # (`?"plasgenomicsutilsR-coordinates"`), so shift here, at the one place positions enter R
  # -- otherwise a scan built from these genotypes sits one base off every IBD table and
  # interval, which breaks exact joins and mis-assigns SNPs on a gene boundary.
  dimnames(mat) <- list(samples, paste0(chrom[idx], ":", pos1[idx] - 1L))
  # `*` is not a base. It says the sequence at this position is deleted on that haplotype --
  # a confident observation, but not one of the alleles being compared, and every k-allele
  # estimator in the package (`he`, pi, Jost's D, the one-hot expansion, `allele_states()`)
  # counts whatever distinct values it finds in this matrix. Left in, a site with 60% `*`
  # reads as a highly diverse site rather than a mostly-deleted one.
  #
  # So blank those calls by default: the haplotype has no base here, which is what NA means
  # everywhere else in the matrix. The allele *indices* are untouched -- `sites$alt` still
  # lists `*` in its own slot -- so index 2 still names the same base it did.
  #
  # `star = "allele"` keeps them, for the deliberate case of treating presence/absence as the
  # state (Pf dimorphic sequence, where the deletion IS the other haplotype).
  n_star <- 0L
  star_col <- 0L
  if (identical(star, "missing")) {
    for (j in seq_len(ncol(mat))) {
      k <- match("*", alt[[idx[j]]])
      if (is.na(k)) next
      # allele_index: exactly the calls that ARE the `*` allele. dosage: every alternate
      # call, because a dosage has already collapsed the alternates and cannot say whether
      # a 2 is the deleted haplotype or the base beside it -- blanking what cannot be
      # attributed is the same rule, applied to a coarser matrix.
      hit <- if (identical(encoding, "allele_index")) !is.na(mat[, j]) & mat[, j] == k
             else !is.na(mat[, j]) & mat[, j] > 0L
      if (!any(hit)) next
      n_star <- n_star + sum(hit)
      star_col <- star_col + 1L
      mat[hit, j] <- NA_integer_
    }
  }

  if (identical(encoding, "dosage") && identical(allele, "ref")) mat <- 2L - mat

  # `alt` is the record's own ALT list, `*` included and in its own order, because that is
  # what an allele index names -- dropping `*` would make index 1 point at the wrong base.
  # `n_alt` counts it faithfully for the same reason. Whether a record is *multiallelic* is a
  # separate question, since `*` is a missingness annotation rather than an allele, so
  # `n_alt_real` answers that one and the two must not be conflated.
  star_site <- vapply(alt[idx], function(a) any(a == "*"), logical(1))
  sites <- data.frame(site_key = colnames(mat), chr = chrom[idx],
                      pos = as.numeric(pos1[idx] - 1L), ref = ref[idx],
                      n_alt = vapply(alt[idx], length, integer(1)),
                      n_alt_real = vapply(alt[idx], function(a) sum(a != "*"), integer(1)),
                      has_spanning_del = star_site,
                      stringsAsFactors = FALSE)
  sites$alt <- alt[idx]
  sites <- sites[c("site_key", "chr", "pos", "ref", "alt", "n_alt", "n_alt_real",
                   "has_spanning_del")]

  n_multi <- sum(sites$n_alt_real > 1L, na.rm = TRUE)
  if (n_multi > 0L) {
    worst <- sites[which.max(sites$n_alt_real), ]
    if (identical(encoding, "dosage")) {
      warning(n_multi, " of ", nrow(sites), " records carry more than one ALT allele (most: ",
              worst$site_key, ", ", worst$n_alt_real, " alternates), and a dosage cannot say ",
              "which one a call carries -- every non-reference genotype there collapses to ",
              "the same number. Use `encoding = \"allele_index\"` to keep them apart.",
              call. = FALSE)
    } else {
      message(n_multi, " of ", nrow(sites), " records carry more than one ALT allele (most: ",
              worst$site_key, ", ", worst$n_alt_real, " alternates); `$genotype` holds allele ",
              "indices, so they stay distinct.")
    }
  }
  if (n_star > 0L)
    message(n_star, " call(s) of the `*` spanning-deletion allele across ", star_col,
            " record(s) set to missing -- `*` marks sequence that is not there rather than ",
            "a base, so counting it as an allele inflates every diversity estimate. Pass ",
            "`star = \"allele\"` to keep them as a state.")

  list(genotype = mat, sample.id = samples, snp.id = vid,
       allele = allele, pruned = isTRUE(prune), positions = "0-based",
       variants = variants, encoding = encoding, star = star, sites = sites)
}

#' Deprecated name for load_genotypes()
#'
#' `run_ld_prune()` was renamed to [load_genotypes()]: the old name reads oddly for what it
#' mostly does, and reads as a contradiction with `prune = FALSE`. Kept so existing scripts
#' keep working.
#'
#' @param ... Passed to [load_genotypes()].
#' @return See [load_genotypes()].
#' @examples
#' \dontrun{
#' # deprecated: use load_genotypes()
#' geno <- run_ld_prune("clean.vcf.gz", gds = "clean.gds")
#' }
#' @export
run_ld_prune <- function(...) {
  warning("`run_ld_prune()` has been renamed to `load_genotypes()`", call. = FALSE)
  load_genotypes(...)
}

# mean-impute missing genotypes per SNP; drop all-missing columns
.impute_geno <- function(mat, what = "pop_structure()") {
  .require_dosage(mat, what)
  mat <- as.matrix(mat)
  cm <- colMeans(mat, na.rm = TRUE)
  keep <- !is.nan(cm)
  mat <- mat[, keep, drop = FALSE]
  cm <- cm[keep]
  na <- which(is.na(mat), arr.ind = TRUE)
  if (nrow(na)) mat[na] <- cm[na[, "col"]]
  mat
}

# ---- PCA + UMAP container --------------------------------------------------

#' Compute PCA and UMAP from a genotype matrix
#'
#' Mean-imputes missing genotypes, runs PCA ([stats::prcomp()]) and, optionally, a
#' UMAP embedding (\pkg{uwot}) with PCA initialisation. Returns a `pop_structure`
#' object the `plot_*()` functions read.
#'
#' @param geno A genotype matrix (samples x SNPs, 0/1/2, `NA` allowed) or the list
#'   returned by [load_genotypes()].
#' @param samples Sample ids (defaults to the genotype list's `sample.id`, or matrix
#'   row names).
#' @param meta Optional per-sample metadata: a data frame with a `sample` column
#'   (plus e.g. `region`, `country`) used for colouring.
#' @param n_pcs Number of PCs to summarise (variance explained).
#' @param umap Compute a UMAP embedding.
#' @param umap_pca,n_neighbors,min_dist UMAP parameters (defaults 30 / 15 / 0.1);
#'   `n_neighbors` is clamped to the sample count.
#' @param seed Random seed for UMAP.
#' @return A `pop_structure` object (a list with `samples`, `pca` scores,
#'   `pca_var`, `umap`, `meta`).
#' @examples
#' # PCA (and optionally UMAP) from a genotype matrix
#' G <- matrix(stats::rbinom(30 * 80, 2, 0.4), 30, 80)
#' rownames(G) <- paste0("s", 1:30)
#' ps <- pop_structure(G, umap = FALSE)
#' dim(ps$pca)
#' head(ps$pca_var)
#' @export
pop_structure <- function(geno, samples = NULL, meta = NULL, n_pcs = 50, umap = TRUE,
                          umap_pca = 30, n_neighbors = 15, min_dist = 0.1, seed = 42) {
  meta <- .normalise_meta(meta)
  if (is.list(geno) && !is.null(geno$genotype)) {
    # asked while the list is still here: once it is a bare matrix, a dosage and an allele
    # index cannot be told apart
    .require_dosage(geno, "pop_structure()")
    if (is.null(samples)) samples <- geno$sample.id
    mat <- geno$genotype
  } else {
    mat <- as.matrix(geno)
    if (is.null(samples)) samples <- rownames(mat)
  }
  if (is.null(samples)) samples <- as.character(seq_len(nrow(mat)))
  mat <- .impute_geno(mat)

  pca <- stats::prcomp(mat, center = TRUE, scale. = FALSE)
  ve <- pca$sdev^2 / sum(pca$sdev^2)
  npc <- min(n_pcs, length(ve))
  pca_var <- data.frame(
    PC = seq_len(npc),
    var_explained = round(ve[seq_len(npc)] * 100, 2),
    cumulative = round(cumsum(ve)[seq_len(npc)] * 100, 2))

  umap_df <- NULL
  if (umap) {
    .need_package("uwot", "the UMAP embedding in pop_structure()")
    nn <- max(2L, min(n_neighbors, nrow(mat) - 1L))
    # umap_pca may be a count (>= 1) or a variance fraction (0 < x < 1)
    npca <- if (umap_pca > 0 && umap_pca < 1) n_pcs_for_variance(pca, umap_pca) else as.integer(umap_pca)
    npca <- min(npca, ncol(mat), nrow(mat) - 1L)
    set.seed(seed)
    emb <- uwot::umap(mat, pca = npca, pca_center = TRUE, n_neighbors = nn,
                      min_dist = min_dist)
    umap_df <- data.frame(sample = samples, UMAP1 = emb[, 1], UMAP2 = emb[, 2],
                          stringsAsFactors = FALSE)
  }
  structure(list(samples = samples, pca = pca$x, prcomp = pca, pca_var = pca_var,
                 umap = umap_df, meta = meta), class = "pop_structure")
}

#' @export
print.pop_structure <- function(x, ...) {
  cat("<pop_structure>", length(x$samples), "samples,",
      ncol(x$pca), "PCs\n")
  cat("  PC1-2 variance:", x$pca_var$var_explained[1], "% /",
      x$pca_var$var_explained[2], "%\n")
  cat("  UMAP:", if (is.null(x$umap)) "-" else "yes",
      "  metadata:", if (is.null(x$meta)) "-" else paste(ncol(x$meta), "cols"), "\n")
  invisible(x)
}

.ps_meta_join <- function(df, meta) {
  if (is.null(meta) || !"sample" %in% names(meta)) return(df)
  merge(df, meta, by = "sample", all.x = TRUE, sort = FALSE)
}

# Grouping values for a set of samples as a factor, carrying the metadata column's own
# level order through. Every plot derives its group ordering from here, so setting levels
# on the metadata (or via PopStructure$set_levels()) drives them all the same way. A
# column that is not a factor is natural-sorted (see .levels_of()).
.group_factor <- function(meta, group, samples = NULL) {
  v <- meta[[group]]
  if (!is.null(samples)) v <- v[match(samples, meta$sample)]
  .as_group_factor(v)
}

# The levels actually present, in the column's level order (empty levels dropped).
.group_order <- function(v) {
  f <- .as_group_factor(v)
  levels(f)[levels(f) %in% as.character(f[!is.na(f)])]
}

# A factor keeps its own level order; anything else gets one from .levels_of().
.as_group_factor <- function(v) {
  if (is.factor(v)) v else factor(as.character(v), levels = .levels_of(v))
}

# ---- plots -----------------------------------------------------------------

#' PCA scatter plot
#'
#' @param x A `pop_structure` object or a [PopStructure] R6 object.
#' @param pcs Which two PCs to plot (default `c(1, 2)`).
#' @param colour,color Metadata column to colour points by (needs `meta`).
#' @param colors,colours Optional named `level -> colour` vector for the colour scale
#'   (e.g. from [meta_colors()]); a `PopStructure` supplies its stored map.
#' @param point_size,point_alpha Point aesthetics.
#' @param legend_point_size Size of the coloured dots in the legend (default `3`, larger
#'   than the plotted points so the key is easy to read); `NULL` leaves it as-is.
#' @return A ggplot object.
#' @examples
#' ps <- example_pop_structure(umap = FALSE)
#' plot_pca(ps, colour = "country")
#' plot_pca(ps, colour = "country", pcs = c(2, 3))
#' @export
plot_pca <- function(x, pcs = c(1, 2), colour = NULL, colors = NULL,
                     point_size = 1.6, point_alpha = 0.8, legend_point_size = 3,
                     color = NULL, colours = NULL) {
  colour <- .alias_arg("colour", "color")
  colors <- .alias_arg("colors", "colours")
  .need_package("ggplot2", "plot_pca()")
  if (inherits(x, "PopStructure")) {
    if (is.null(colors) && !is.null(colour)) colors <- x$get_colors()[[colour]]
    x <- x$as_ps()
  }
  a <- pcs[1]; b <- pcs[2]
  df <- data.frame(sample = x$samples, .x = x$pca[, a], .y = x$pca[, b],
                   stringsAsFactors = FALSE)
  df <- .ps_meta_join(df, x$meta)
  ve <- x$pca_var$var_explained
  p <- ggplot2::ggplot(df, ggplot2::aes(.data$.x, .data$.y)) +
    ggplot2::labs(x = sprintf("PC%d (%.1f%%)", a, ve[a]),
                  y = sprintf("PC%d (%.1f%%)", b, ve[b])) +
    ggplot2::theme_minimal(base_size = 11)
  p + .scatter_points(colour, point_size, point_alpha, colors, legend_point_size)
}

#' UMAP scatter plot
#'
#' @param x A `pop_structure` object (built with `umap = TRUE`) or a [PopStructure].
#' @param colour,color Metadata column to colour points by (needs `meta`).
#' @param colors,colours Optional named `level -> colour` vector for the colour scale
#'   (e.g. from [meta_colors()]); a `PopStructure` supplies its stored map.
#' @param point_size,point_alpha Point aesthetics.
#' @param legend_point_size Size of the coloured dots in the legend (default `3`, larger
#'   than the plotted points so the key is easy to read); `NULL` leaves it as-is.
#' @param point_border Outline colour for each point (e.g. `"black"`), or `NULL`
#'   (default) for unoutlined points. Outlined points use shape 21, so the
#'   categories drive `fill` rather than `colour`.
#' @param point_stroke Outline width (default `0.3`); ignored when `point_border`
#'   is `NULL`.
#' @param legend_title Title for the colour legend (default: the `colour` column name).
#' @return A ggplot object.
#' @examples
#' ps <- example_pop_structure(umap = FALSE)
#' ps$run_umap(pca_components = 5)
#' plot_umap(ps, colour = "country")
#' @export
plot_umap <- function(x, colour = NULL, colors = NULL, point_size = 1.6,
                      point_alpha = 0.8, legend_point_size = 3,
                      point_border = NULL, point_stroke = 0.3,
                      legend_title = NULL,
                      color = NULL, colours = NULL) {
  colour <- .alias_arg("colour", "color")
  colors <- .alias_arg("colors", "colours")
  .need_package("ggplot2", "plot_umap()")
  if (inherits(x, "PopStructure")) {
    if (is.null(colors) && !is.null(colour)) colors <- x$get_colors()[[colour]]
    x <- x$as_ps()
  }
  if (is.null(x$umap)) stop("this pop_structure has no UMAP (build with umap = TRUE)",
                            call. = FALSE)
  df <- .ps_meta_join(x$umap, x$meta)
  df$.x <- df$UMAP1; df$.y <- df$UMAP2
  ggplot2::ggplot(df, ggplot2::aes(.data$.x, .data$.y)) +
    ggplot2::labs(x = "UMAP 1", y = "UMAP 2") +
    ggplot2::coord_equal() +
    ggplot2::theme_minimal(base_size = 11) +
    .scatter_points(colour, point_size, point_alpha, colors, legend_point_size,
                    point_border, point_stroke, legend_title)
}

# shared point + colour layers for the scatter plots
#
# `border` outlines each point: that needs shape 21, where the category drives `fill` and
# the outline owns `colour`, so the whole scale/guide/label set moves from colour to fill.
# With `border = NULL` (the default) nothing changes -- same shape 19 layer as before.
.scatter_points <- function(colour, point_size, point_alpha, colors = NULL,
                            legend_point_size = NULL, border = NULL, stroke = 0.3,
                            legend_title = NULL) {
  bordered <- !is.null(border)
  if (is.null(colour)) {
    return(list(if (bordered)
      ggplot2::geom_point(size = point_size, alpha = point_alpha, shape = 21,
                          colour = border, stroke = stroke)
      else ggplot2::geom_point(size = point_size, alpha = point_alpha)))
  }
  aes_key <- if (bordered) "fill" else "colour"
  out <- list(
    if (bordered)
      ggplot2::geom_point(ggplot2::aes(fill = .data[[colour]]), shape = 21,
                          colour = border, stroke = stroke,
                          size = point_size, alpha = point_alpha)
    else
      ggplot2::geom_point(ggplot2::aes(colour = .data[[colour]]),
                          size = point_size, alpha = point_alpha),
    # the key the category is mapped to differs with `border`, so name the legend by it
    ggplot2::labs(!!aes_key := if (is.null(legend_title)) colour else legend_title))
  if (!is.null(colors))
    out <- c(out, list(if (bordered) ggplot2::scale_fill_manual(values = colors)
                       else ggplot2::scale_colour_manual(values = colors)))
  if (!is.null(legend_point_size))
    out <- c(out, list(ggplot2::guides(
      !!aes_key := ggplot2::guide_legend(override.aes = list(size = legend_point_size)))))
  out
}

# ---- sNMF admixture --------------------------------------------------------

#' Run sNMF (LEA) admixture on a genotype matrix
#'
#' Writes a `.geno` file (missing coded as 9) and runs [LEA::snmf()] over a range of
#' `K`, returning the sNMF project. Cross-entropy is computed so [snmf_best_k()] can
#' pick `K`.
#'
#' @param geno A genotype matrix (samples x SNPs, 0/1/2, `NA` allowed) or the list
#'   from [load_genotypes()].
#' @param K Integer vector of ancestral-population counts to fit (default `1:10`).
#' @param rep Repetitions per `K`.
#' @param alpha Regularisation.
#' @param seed Random seed.
#' @param cpu CPU cores.
#' @param cache Reuse a previously computed project when the genotypes and all
#'   parameters match (sNMF is slow); set `FALSE` to always recompute.
#' @param cache_dir Directory for the cached `.geno`/project (default a per-session
#'   temp dir). Point it at a persistent path to reuse across sessions.
#' @param verbose Show LEA's (voluminous) console output. Default `FALSE` runs it
#'   quietly; set `log_file` to capture it instead of discarding it.
#' @param log_file Optional file to write LEA's output to when `verbose = FALSE`.
#' @return An `snmf_fit` object: a list with the LEA `project`, the fitted `K` range,
#'   and `samples`.
#' @examples
#' \dontrun{
#' fit <- run_snmf(geno, K = 1:10, rep = 10, cache_dir = "snmf_cache")
#' snmf_best_k(fit)
#' }
#' @export
run_snmf <- function(geno, K = 1:10, rep = 10, alpha = 10, seed = 42, cpu = 1,
                     cache = TRUE, cache_dir = NULL, verbose = FALSE, log_file = NULL) {
  .need_package("LEA", "run_snmf()")
  if (is.list(geno) && !is.null(geno$genotype)) {
    samples <- geno$sample.id
    mat <- geno$genotype
  } else {
    mat <- as.matrix(geno)
    samples <- rownames(mat)
  }
  K <- as.integer(K)
  if (is.null(cache_dir)) cache_dir <- file.path(tempdir(), "plas_snmf")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  key <- substr(rlang::hash(list(mat, K, rep, alpha, seed)), 1, 16)
  geno_file <- file.path(cache_dir, paste0("snmf_", key, ".geno"))
  proj_file <- sub("\\.geno$", ".snmfProject", geno_file)
  samp_file <- file.path(cache_dir, paste0("snmf_", key, ".samples.rds"))

  if (cache && file.exists(proj_file)) {
    project <- .run_quiet(function() LEA::load.snmfProject(proj_file), verbose, log_file)
    if (file.exists(samp_file)) samples <- readRDS(samp_file)
  } else {
    .require_dosage(mat, "run_snmf()")
    mat[is.na(mat)] <- 9L
    project <- .run_quiet(function() {
      LEA::write.geno(mat, output.file = geno_file)
      LEA::snmf(geno_file, K = K, repetitions = rep, alpha = alpha, entropy = TRUE,
                project = "new", seed = seed, CPU = cpu)
    }, verbose, log_file)
    if (!is.null(samples)) saveRDS(samples, samp_file)
  }
  structure(list(project = project, K = K, samples = samples, geno_file = geno_file),
            class = "snmf_fit")
}

# run a thunk quietly (LEA prints a lot): discard stdout, or send it to `log_file`.
.run_quiet <- function(fn, verbose = FALSE, log_file = NULL) {
  if (verbose) return(fn())
  con <- if (is.null(log_file)) nullfile() else log_file
  out <- NULL
  utils::capture.output(out <- fn(), file = con, type = "output")
  out
}

#' Number of PCs explaining a target cumulative variance
#'
#' Handy for setting how many principal components feed the UMAP: pass the fraction of
#' variance you want the PCA step to capture.
#'
#' @param x A [stats::prcomp()] result or a numeric vector of eigenvalues/variances.
#' @param target Cumulative variance to reach, as a proportion (`0.1`) or a percent
#'   (`10`).
#' @return The smallest number of leading PCs whose cumulative variance >= `target`.
#' @examples
#' ps <- example_pop_structure(umap = FALSE)
#' # how many PCs to feed UMAP for a given share of the variance
#' n_pcs_for_variance(ps$prcomp(), 0.5)
#' n_pcs_for_variance(ps$prcomp(), 50)     # percent is accepted too
#' @export
n_pcs_for_variance <- function(x, target = 0.8) {
  ve <- if (inherits(x, "prcomp")) x$sdev^2 / sum(x$sdev^2)
        else if (is.numeric(x)) x / sum(x)
        else stop("`x` must be a prcomp result or a numeric variance vector", call. = FALSE)
  if (target > 1) target <- target / 100
  w <- which(cumsum(ve) >= target)
  if (length(w)) w[1] else length(ve)
}

#' @export
print.snmf_fit <- function(x, ...) {
  cat("<snmf_fit>", length(x$samples), "samples,  K =",
      paste(range(x$K), collapse = "-"), "\n")
  invisible(x)
}

.snmf_project <- function(x) if (inherits(x, "snmf_fit")) x$project else x

#' Pick the best K from an sNMF fit by cross-entropy
#'
#' sNMF fits `rep` replicates at each K, so a K has a spread of cross-entropies rather than
#' one. They are summarised by their **minimum** across replicates, which is what
#' \pkg{LEA}'s own `plot()` of an sNMF project shows, and the criterion its vignette reads a
#' best K off. It is also the replicate [snmf_q()] returns: sNMF's objective is non-convex,
#' replicates land in different local optima, and the best-fitting one is the model whose
#' ancestry you go on to plot -- so comparing minima compares the models actually used at
#' each K.
#'
#' `stat = "mean"` averages the replicates instead. Reach for it when the replicate count
#' differs between K values (`n_runs` in [snmf_cross_entropy()]), since a minimum over more
#' replicates is expected to be smaller whether or not that K fits better; averaging is not
#' sensitive to how many were run.
#'
#' @param x An [run_snmf()] result (or a raw LEA project, with `K` given).
#' @param K Candidate K values (defaults to the fitted range from [run_snmf()]).
#' @param stat Combine replicates by `"min"` (default) or `"mean"` cross-entropy.
#' @return The K minimising the summarised cross-entropy.
#' @seealso [snmf_cross_entropy()] for the per-K spread, [plot_snmf_cross_entropy()] for the
#'   elbow.
#' @examples
#' \dontrun{
#' snmf_best_k(fit)                  # minimum cross-entropy, as LEA plots it
#' snmf_best_k(fit, stat = "mean")   # for an uneven number of replicates per K
#' }
#' @export
snmf_best_k <- function(x, K = NULL, stat = c("min", "mean")) {
  .need_package("LEA", "snmf_best_k()")
  stat <- match.arg(stat)
  project <- .snmf_project(x)
  if (is.null(K)) {
    if (inherits(x, "snmf_fit")) K <- x$K
    else stop("pass K (the fitted range) when giving a raw sNMF project", call. = FALSE)
  }
  fn <- if (stat == "min") min else mean
  vals <- vapply(K, function(k) fn(LEA::cross.entropy(project, K = k)), numeric(1))
  K[which.min(vals)]
}

#' Cross-entropy of every sNMF replicate, summarised per K
#'
#' sNMF fits `rep` independent replicates at each K and scores each by cross-entropy
#' (lower is better). [snmf_q()] and [plot_admixture()] use the **minimum**-cross-entropy
#' replicate, so the `min` column is the one that describes the ancestry actually plotted;
#' `mean` and `max` show how much the replicates disagreed, which is worth a look before
#' trusting a K.
#'
#' `min` is also what K is chosen on, by [snmf_best_k()] and by the
#' [plot_snmf_cross_entropy()] elbow, so the K you settle on and the ancestry you draw at it
#' are scored by the same number. A curve that keeps falling or is flat means the data do not
#' support a well-defined K, whatever [snmf_best_k()] returns; `n_runs` is worth a glance
#' first, since a K whose replicates mostly failed has its minimum taken over fewer of them.
#'
#' @param x An [run_snmf()] result (or a raw LEA project, with `K` given).
#' @param K Candidate K values (defaults to the fitted range).
#' @return A tibble, one row per K: `K`, `n_runs`, `min`, `mean`, `max`, and `best_run`
#'   (the replicate index attaining `min`, i.e. the one [snmf_q()] returns).
#' @seealso [plot_snmf_cross_entropy()], [snmf_best_k()], [snmf_q()]
#' @examples
#' \dontrun{
#' snmf_cross_entropy(fit)           # one row per K: min, mean, and n_runs
#' }
#' @export
snmf_cross_entropy <- function(x, K = NULL) {
  .need_package("LEA", "snmf_cross_entropy()")
  project <- .snmf_project(x)
  if (is.null(K)) {
    if (inherits(x, "snmf_fit")) K <- x$K
    else stop("pass K (the fitted range) when giving a raw sNMF project", call. = FALSE)
  }
  rows <- lapply(K, function(k) {
    ce <- as.numeric(LEA::cross.entropy(project, K = k))
    ce <- ce[is.finite(ce)]
    if (!length(ce)) {
      return(data.frame(K = k, n_runs = 0L, min = NA_real_, mean = NA_real_,
                        max = NA_real_, best_run = NA_integer_))
    }
    data.frame(K = k, n_runs = length(ce), min = min(ce), mean = mean(ce),
               max = max(ce), best_run = which.min(ce))
  })
  tibble::as_tibble(do.call(rbind, rows))
}

#' Cross-entropy elbow plot for choosing K
#'
#' Cross-entropy against K, so the elbow (or the absence of one) is visible. The line
#' follows `stat`, and `show_range = TRUE` adds the replicate min-max band -- a wide band
#' means the replicates disagreed and that K is not reproducible.
#'
#' @param x An [run_snmf()] result, or a table from [snmf_cross_entropy()].
#' @param K Candidate K values (defaults to the fitted range).
#' @param stat Which summary the line follows, and which the red marker minimises:
#'   `"min"` (default, as in [snmf_best_k()]) or `"mean"`. See [snmf_best_k()] for why
#'   `"min"` is the default and when `"mean"` is worth asking for.
#' @param show_range Draw the replicate min-max band (default `TRUE`).
#' @param best_k K to mark in red; `NULL` (default) marks the K minimising `stat`, `NA`
#'   marks none.
#' @param point_size,line_width Point and line sizes.
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' fit <- run_snmf(geno, K = 1:10, rep = 10)
#' snmf_cross_entropy(fit)          # the numbers
#' plot_snmf_cross_entropy(fit)     # the elbow
#' }
#' @export
plot_snmf_cross_entropy <- function(x, K = NULL, stat = c("min", "mean"),
                                    show_range = TRUE, best_k = NULL,
                                    point_size = 2.4, line_width = 0.6) {
  .need_package("ggplot2", "plot_snmf_cross_entropy()")
  stat <- match.arg(stat)
  ce <- if (is.data.frame(x) && all(c("K", "min", "mean") %in% names(x))) x
        else snmf_cross_entropy(x, K)
  ce <- ce[is.finite(ce[[stat]]), , drop = FALSE]
  if (!nrow(ce)) stop("no finite cross-entropy values to plot", call. = FALSE)
  if (is.null(best_k)) best_k <- ce$K[which.min(ce[[stat]])]
  ce$.y <- ce[[stat]]
  ce$.best <- !is.na(best_k) & ce$K == best_k

  p <- ggplot2::ggplot(ce, ggplot2::aes(.data$K, .data$.y))
  if (show_range && any(is.finite(ce$max))) {
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$min, ymax = .data$max),
                                  fill = "grey80", alpha = 0.5)
  }
  p <- p +
    ggplot2::geom_line(linewidth = line_width, colour = "grey30") +
    ggplot2::geom_point(ggplot2::aes(colour = .data$.best), size = point_size) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey25"),
                                 guide = "none") +
    ggplot2::scale_x_continuous(breaks = ce$K) +
    ggplot2::labs(x = "K", y = sprintf("cross-entropy (%s of replicates)", stat),
                  title = "sNMF cross-entropy by K") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  if (!is.na(best_k) && any(ce$.best)) {
    p <- p + ggplot2::annotate("text", x = best_k, y = ce$.y[ce$.best],
                               label = paste0("  best K = ", best_k),
                               hjust = 0, vjust = -0.9, colour = "firebrick", size = 3.4)
  }
  p
}

#' Admixture bar plots across every K, as pages
#'
#' One [plot_admixture()] per K, returned as a named list -- hand it straight to
#' [save_plot()] for a multi-page PDF, one K per page. Optionally leads with the
#' cross-entropy elbow ([plot_snmf_cross_entropy()]) so the page that tells you which K to
#' believe comes first, with the best K marked in red.
#'
#' Bars stay in the same place from page to page when a shared `sample_order` is used,
#' which is what makes the pages comparable: a sample sits at the same x on every K.
#'
#' @param x A [PopStructure] with a fitted sNMF, or an [run_snmf()] result (then give
#'   `meta` and `samples`).
#' @param K Which K values to draw (default: every fitted K with a Q matrix, K >= 2 --
#'   K = 1 is a single block and carries no information).
#' @param group Metadata column to facet by (e.g. `"region"`).
#' @param cross_entropy_first Lead with the cross-entropy elbow page (default `TRUE`).
#' @param sample_order Explicit sample order shared by every page. Overrides
#'   `sample_order_best_k`.
#' @param sample_order_best_k Derive one shared sample order from the best K's Q and use
#'   it on every page (default `TRUE`). `FALSE` lets each page cluster its own samples,
#'   so bars move between pages.
#' @param best_k The K to treat as best; `NULL` (default) uses [snmf_best_k()].
#' @param stat How [snmf_best_k()] and the elbow combine replicates: `"min"` (default)
#'   or `"mean"`.
#' @param meta,samples Metadata and sample ids, needed only when `x` is a raw
#'   [run_snmf()] result.
#' @param ... Passed to [plot_admixture()] (e.g. `group_bar`, `border`, `colours`).
#' @return A named list of ggplots: `"cross_entropy"` (when asked for) then `"K=2"`,
#'   `"K=3"`, ... The best K's page is titled as such.
#' @examples
#' \dontrun{
#' ps$run_snmf(K = 1:12)
#' pages <- ps$plot_admixture_multi_k(group = "region")
#' save_plot("admixture_all_k.pdf", pages)     # one K per page
#' }
#' @seealso [plot_snmf_cross_entropy()], [plot_admixture()], [save_plot()]
#' @export
plot_admixture_multi_k <- function(x, K = NULL, group = NULL,
                                   cross_entropy_first = TRUE,
                                   sample_order = NULL, sample_order_best_k = TRUE,
                                   best_k = NULL, stat = c("min", "mean"),
                                   meta = NULL, samples = NULL, ...) {
  meta <- .normalise_meta(meta)
  .need_package("ggplot2", "plot_admixture_multi_k()")
  stat <- match.arg(stat)
  is_ps <- inherits(x, "PopStructure")
  fit <- if (is_ps) x$get_snmf_fit() else x
  if (is.null(fit)) stop("run_snmf() first", call. = FALSE)
  if (is.null(meta)) meta <- if (is_ps) x$get_meta() else NULL
  if (is.null(samples)) samples <- if (is_ps) x$get_samples() else fit$samples

  fitted_k <- if (inherits(fit, "snmf_fit")) fit$K else K
  if (is.null(K)) K <- fitted_k[fitted_k >= 2]
  K <- as.integer(K)
  if (!length(K)) stop("no K values to draw (K = 1 alone carries no structure)", call. = FALSE)
  if (is.null(best_k)) best_k <- snmf_best_k(fit, K = fitted_k, stat = stat)

  q_at <- function(k) {
    qq <- snmf_q(fit, K = k)
    if (!is.null(samples)) qq <- qq[rownames(qq) %in% samples, , drop = FALSE]
    qq
  }
  # one shared order keeps a sample at the same x on every page
  if (is.null(sample_order) && isTRUE(sample_order_best_k)) {
    sample_order <- admixture_order(q_at(best_k), meta = meta, group = group)
  }

  pages <- list()
  if (isTRUE(cross_entropy_first)) {
    pages[["cross_entropy"]] <- plot_snmf_cross_entropy(fit, K = fitted_k, stat = stat,
                                                        best_k = best_k)
  }
  for (k in K) {
    ttl <- if (!is.na(best_k) && k == best_k) sprintf("K = %d  (best by %s cross-entropy)", k, stat)
           else sprintf("K = %d", k)
    p <- plot_admixture(q_at(k), meta = meta, group = group,
                        sample_order = sample_order, ...) +
      ggplot2::labs(title = ttl)
    pages[[sprintf("K=%d", k)]] <- p
  }
  pages
}

#' Best-run Q (ancestry proportion) matrix for a given K
#'
#' @param x An [run_snmf()] result (or a raw LEA project).
#' @param K The number of ancestral populations.
#' @param run Which replicate; defaults to the lowest cross-entropy run.
#' @return A samples-by-K matrix of ancestry proportions, with sample ids as row names
#'   when available.
#' @examples
#' \dontrun{
#' q <- snmf_q(fit, K = 5)           # the best-fitting replicate at K = 5
#' dim(q)
#' }
#' @export
snmf_q <- function(x, K, run = NULL) {
  .need_package("LEA", "snmf_q()")
  project <- .snmf_project(x)
  if (is.null(run)) run <- which.min(LEA::cross.entropy(project, K = K))
  q <- LEA::Q(project, K = K, run = run)
  colnames(q) <- paste0("K", seq_len(ncol(q)))
  if (inherits(x, "snmf_fit") && !is.null(x$samples)) rownames(q) <- x$samples
  q
}

# A legend of `n` keys, wrapped so it does not run off the canvas: down the side it grows
# into extra columns, along the bottom into extra rows. `rows` caps the keys per column
# (or per row when horizontal); NULL uses a default that keeps a right-hand legend inside a
# typical page.
.LEGEND_MAX_KEYS <- 10L

.LEGEND_MAX_ROWS_HORIZONTAL <- 2L

.legend_wrap <- function(n, rows = NULL, position = "right", order = 0) {
  if (identical(position, "none")) return("none")
  if (position %in% c("bottom", "top")) {
    # along the bottom, keep it shallow so it does not eat the panel's height
    nr <- if (is.null(rows)) .LEGEND_MAX_ROWS_HORIZONTAL else max(1L, as.integer(rows))
    return(ggplot2::guide_legend(order = order, nrow = min(n, nr)))
  }
  per <- if (is.null(rows)) .LEGEND_MAX_KEYS else max(1L, as.integer(rows))
  ggplot2::guide_legend(order = order, ncol = max(1L, ceiling(n / per)))
}

#' Sample order for admixture bars
#'
#' Orders samples within each `group` by hierarchical clustering of their ancestry
#' vectors. Compute it once (e.g. at the best K) and pass it as `sample_order` to
#' [plot_admixture()] for every K, so bars stay in the same position across K.
#'
#' @param q A samples-by-K ancestry matrix.
#' @param samples Sample ids (default `rownames(q)`).
#' @param meta,group Optional metadata + the column to order within.
#' @return An ordered character vector of sample ids.
#' @examples
#' q <- matrix(c(0.9, 0.1, 0.2, 0.8, 0.85, 0.15, 0.1, 0.9), ncol = 2, byrow = TRUE)
#' rownames(q) <- paste0("s", 1:4)
#' meta <- data.frame(sample = rownames(q), region = c("A", "A", "B", "B"))
#' admixture_order(q, meta = meta, group = "region")
#' @export
admixture_order <- function(q, samples = NULL, meta = NULL, group = NULL) {
  meta <- .normalise_meta(meta)
  q <- as.matrix(q)
  if (is.null(samples)) samples <- rownames(q)
  if (is.null(samples)) samples <- as.character(seq_len(nrow(q)))
  grp <- factor(rep("all", length(samples)))
  if (!is.null(meta) && !is.null(group) && "sample" %in% names(meta) && group %in% names(meta)) {
    grp <- .group_factor(meta, group, samples)
  }
  ord <- character(0)
  # groups are swept in the column's level order, so samples sit in the requested order
  for (g in .group_order(grp)) {
    idx <- which(!is.na(grp) & grp == g)
    if (length(idx) > 2) {
      hc <- stats::hclust(stats::dist(q[idx, , drop = FALSE]), method = "ward.D2")
      idx <- idx[hc$order]
    }
    ord <- c(ord, samples[idx])
  }
  # samples whose group is not one of the levels are left out (set_levels() drops unlisted
  # levels), which silently shrinks the plot -- so say how many went
  if (anyNA(grp)) {
    warning(sum(is.na(grp)), " sample(s) have no level in the grouping column and are ",
            "not plotted; add their level with set_levels() to keep them", call. = FALSE)
  }
  ord
}

#' Admixture (STRUCTURE) bar plot from a Q matrix
#'
#' Stacked ancestry-proportion bars, one per sample, optionally faceted by a grouping
#' column, with an optional group colour strip whose colours can match another plot
#' (e.g. UMAP points).
#'
#' @param q A samples-by-K ancestry matrix (e.g. from [snmf_q()]).
#' @param samples Sample ids (defaults to `rownames(q)`).
#' @param meta Optional metadata data frame with a `sample` column.
#' @param group Optional metadata column to facet by (e.g. `"region"`).
#' @param order_within Order samples within each group by clustering (ignored when
#'   `sample_order` is supplied).
#' @param sample_order Optional explicit sample order (see [admixture_order()]); keeps
#'   bars in the same place across different K.
#' @param colours,colors Optional fill colours for the K clusters.
#' @param group_bar Draw a coloured strip above the bars keyed by `group`.
#' @param group_colours,group_colors Named `level -> colour` vector for the group strip (see
#'   [meta_colors()]); pass the same mapping you use for the UMAP to match colours.
#' @param border Outline each sample's bar (default `TRUE`, matching
#'   [plot_structure_figure()]) so neighbours with nearly identical ancestry stay distinct;
#'   `FALSE` for borderless bars. With very many samples in a narrow render the outlines
#'   can swamp the fills -- either render wider ([save_plot()] uses the attached width) or
#'   set `border = FALSE` / a thinner `border_linewidth`.
#' @param border_colour,border_color,border_linewidth Colour and width of the per-sample outline.
#' @param legend_position Where the legends go: `"right"` (default), `"bottom"`, `"top"`,
#'   `"left"`, or `"none"`. A large K makes for a tall legend stack, so `"bottom"` often
#'   fits better on a wide, short admixture panel.
#' @param legend_rows Keys per legend column (or per row when the legend is horizontal)
#'   before wrapping into another column/row. `NULL` (default) wraps a side legend every
#'   `r plasgenomicsutilsR:::.LEGEND_MAX_KEYS` keys and splits a horizontal one over two
#'   rows, which keeps `K` = 15 plus a group strip on the page. The suggested output height
#'   accounts for whatever this works out to.
#' @param cluster_label Legend title for the ancestry fills. The K components are what
#'   sNMF calls clusters, but next to a UMAP -- where the visible groupings are also
#'   clusters -- a legend reading "cluster" invites reading the two as the same thing.
#'   The default wraps over two lines so the longer wording costs no legend width.
#' @return A ggplot object.
#' @examples
#' q <- matrix(c(0.9, 0.1, 0.2, 0.8, 0.85, 0.15, 0.1, 0.9), ncol = 2, byrow = TRUE)
#' rownames(q) <- paste0("s", 1:4)
#' meta <- data.frame(sample = rownames(q), region = c("A", "A", "B", "B"))
#' plot_admixture(q, meta = meta, group = "region")
#' @export
plot_admixture <- function(q, samples = NULL, meta = NULL, group = NULL,
                           order_within = TRUE, sample_order = NULL, colours = NULL,
                           group_bar = FALSE, group_colours = NULL,
                           border = TRUE, border_colour = "black",
                           border_linewidth = 0.15,
                           legend_position = c("right", "bottom", "top", "left", "none"),
                           legend_rows = NULL, cluster_label = "Ancestry\ncomponent",
                           colors = NULL, group_colors = NULL, border_color = NULL) {
  meta <- .normalise_meta(meta)
  colours <- .alias_arg("colours", "colors")
  group_colours <- .alias_arg("group_colours", "group_colors")
  border_colour <- .alias_arg("border_colour", "border_color")
  .need_package("ggplot2", "plot_admixture()")
  legend_position <- match.arg(legend_position)
  q <- as.matrix(q)
  if (is.null(colnames(q))) colnames(q) <- paste0("K", seq_len(ncol(q)))
  if (is.null(samples)) samples <- rownames(q)
  if (is.null(samples)) samples <- as.character(seq_len(nrow(q)))
  df <- data.frame(sample = samples, q, check.names = FALSE, stringsAsFactors = FALSE)
  df <- .ps_meta_join(df, meta)

  if (is.null(sample_order)) {
    sample_order <- if (order_within) admixture_order(q, samples, meta, group) else samples
  }
  df <- df[match(sample_order, df$sample), , drop = FALSE]
  df$sample <- factor(df$sample, levels = sample_order)

  long <- stats::reshape(
    df[, c("sample", if (!is.null(group)) group, colnames(q))],
    direction = "long", varying = colnames(q), v.names = "q",
    times = colnames(q), timevar = "cluster", idvar = "sample")
  long$sample <- factor(long$sample, levels = sample_order)
  # reshape() leaves `cluster` a character, so K10 would sort before K2. Keep the Q
  # matrix's own column order, which is K1..K<K>.
  long$cluster <- factor(long$cluster, levels = colnames(q))

  cl_cols <- if (is.null(colours)) .pick_palette(ncol(q)) else colours
  p <- ggplot2::ggplot(long, ggplot2::aes(.data$sample, .data$q, fill = .data$cluster)) +
    ggplot2::geom_col(width = 1, position = "stack",
                      colour = if (border) border_colour else NA,
                      linewidth = border_linewidth) +
    ggplot2::scale_fill_manual(values = cl_cols, drop = FALSE,
                               guide = .legend_wrap(ncol(q), legend_rows, legend_position,
                                                    order = 1)) +
    ggplot2::labs(x = NULL, y = "ancestry proportion", fill = cluster_label) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank())

  if (group_bar && !is.null(group) && group %in% names(df)) {
    .need_package("ggnewscale", "the group colour bar")
    # the strip's facet column has to stay the same factor as the bars': a character copy
    # here makes ggplot combine the two layers' facet values alphabetically
    gdf <- data.frame(sample = factor(df$sample, levels = sample_order),
                      grp = .as_group_factor(df[[group]]), stringsAsFactors = FALSE)
    gdf[[group]] <- gdf$grp   # carry the facet variable so facet_grid subsets the strip
    if (is.null(group_colours)) group_colours <- meta_colors(df, cols = group)[[group]]
    p <- p + ggnewscale::new_scale_fill() +
      ggplot2::geom_tile(data = gdf, ggplot2::aes(.data$sample, y = 1.06, fill = .data$grp),
                         height = 0.05, inherit.aes = FALSE) +
      ggplot2::scale_fill_manual(values = group_colours, name = group,
                                 guide = .legend_wrap(length(group_colours), legend_rows,
                                                      legend_position, order = 2)) +
      ggplot2::coord_cartesian(clip = "off")
  }
  p <- p + ggplot2::theme(legend.position = legend_position,
                          legend.box = if (legend_position %in% c("bottom", "top"))
                            "horizontal" else "vertical")
  n_groups <- if (!is.null(group) && group %in% names(df)) length(unique(df[[group]])) else 1L
  if (!is.null(group) && group %in% names(df)) {
    p <- p + ggplot2::facet_grid(stats::as.formula(paste("~", group)),
                                 scales = "free_x", space = "free")
  }
  # Width scales with the number of bars so save_plot() makes it wide enough to read.
  # Height has to clear a side legend as well, or a large K clips its own keys -- measure
  # the built legend rather than estimating from the key count, since key size, text size
  # and how many columns the guide wrapped into all contribute.
  # Measured on a throwaway device: measuring on the current one leaves a stray Rplots.pdf
  # in a script and an empty figure in a notebook chunk (see .with_null_device()).
  panel_h <- 4
  height <- .with_null_device({
    gt <- try(ggplot2::ggplotGrob(p), silent = TRUE)
    if (inherits(gt, "try-error")) panel_h
    else if (legend_position %in% c("right", "left")) {
      # a side legend spans the canvas, so the canvas has to be at least that tall
      max(panel_h, .guide_box_height(gt, c("guide-box-right", "guide-box-left")) + 0.25)
    } else {
      # a legend along the bottom/top is a layout row: it adds to the height instead of
      # competing with it, so give the panel its space on top of the legend's
      panel_h + .guide_box_height(gt, c("guide-box-bottom", "guide-box-top"))
    }
  })
  attr(p, "plasgenomics_dims") <- c(
    width = round(max(6, min(24, nrow(q) * 0.06 + n_groups * 0.3)), 1),
    height = round(height, 1))
  p
}

# ---- PopStructure R6 wrapper ----------------------------------------------

#' Population-structure workspace (PCA + UMAP + admixture)
#'
#' An R6 object that bundles a genotype matrix, its PCA (the full [stats::prcomp()]
#' result), an optional UMAP embedding, per-sample metadata, a shared metadata colour
#' map, and an sNMF admixture fit. Because it keeps the fitted objects, you can save it
#' (`saveRDS`) and re-plot without recomputing, colour PCA/UMAP/admixture consistently
#' (same `level -> colour` map), reorder admixture bars once and reuse the order across
#' K, and sub-select to a set of samples (or a metadata match) for output.
#'
#' @examples
#' ps <- example_pop_structure()
#' ps$run_umap(pca_components = 0.5)      # PCs covering 50% of variance
#' \dontrun{
#' ps$run_snmf(K = 1:6)                   # cached + quiet
#' ps$plot_admixture(K = ps$best_k(), group = "population")
#' west <- ps$subset(population = c("PopA", "PopB"))
#' }
#' @export
PopStructure <- R6::R6Class("PopStructure",
  public = list(
    #' @description Build from a genotype matrix (or a [load_genotypes()] list).
    #' @param geno Genotype matrix (samples x SNPs, 0/1/2, `NA`) or `load_genotypes()` list.
    #' @param samples Sample ids (default row names / the list's `sample.id`).
    #' @param meta Optional metadata (data frame with a `sample` column).
    #' @param n_pcs Number of PCs to summarise.
    #' @param colors Optional named list of colour maps (see [meta_colors()]).
    #' @param one_based The genotype column names carry 1-based VCF positions (a matrix built
    #'   by an older `load_genotypes()`, or read straight off a VCF); shift them to the 0-based
    #'   convention on the way in.
    #' @param pruned Whether the SNPs were LD-pruned. Taken from a [load_genotypes()] list
    #'   when it says so. Pruning is right for PCA / admixture and wrong for looking at
    #'   haplotypes, since it drops SNPs precisely for being correlated with a neighbour.
    #' @param allele Which allele the dosages count, `"alt"` or `"ref"`. Taken from a
    #'   [load_genotypes()] list when it says so; a bare matrix cannot say, and the two codings
    #'   are indistinguishable afterwards, so anything that names the calls has to be told.
    #' @param pca_panel Which derived view PCA and UMAP run on when the panel holds allele
    #'   indices. `"onehot"` gives every ALT its own indicator column, so a multiallelic site
    #'   still reaches the ordination; `"dosage"` drops those sites, which is what the object
    #'   did before. `"auto"` (default) takes one-hot when the panel holds a multiallelic site
    #'   and dosage otherwise -- so a wholly biallelic panel's PCA is unchanged, because on a
    #'   biallelic site the indicator column *is* the alt-dosage column.
    initialize = function(geno, samples = NULL, meta = NULL, n_pcs = 50, colors = NULL,
                          allele = NULL, pruned = NULL, full = NULL,
                          one_based = FALSE, colours = NULL,
                          pca_panel = c("auto", "onehot", "dosage")) {
      pca_panel <- match.arg(pca_panel)
      colors <- .alias_arg("colors", "colours")
      meta <- .normalise_meta(meta)
      positions <- NULL
      sites <- NULL
      enc <- NULL
      if (is.list(geno) && !is.null(geno$genotype)) {
        if (is.null(samples)) samples <- geno$sample.id
        if (is.null(allele)) allele <- geno$allele
        if (is.null(pruned)) pruned <- geno$pruned
        positions <- geno$positions
        sites <- geno$sites
        enc <- geno$encoding
        geno <- geno$genotype
      }
      # A matrix built before positions were shifted (or straight off a VCF) carries 1-based
      # ids; move them once here so nothing downstream has to ask.
      if (isTRUE(one_based) && !is.null(colnames(geno))) {
        ids <- .parse_snp_ids(colnames(geno))
        colnames(geno) <- paste0(ids$chr, ":", as.integer(ids$pos) - 1L)
        positions <- "0-based"
      }
      private$snp_positions <- positions %||% if (isTRUE(one_based)) "0-based" else NULL
      # shifted to 0-based above if it needed it, so the keys still match the matrix
      if (!is.null(sites) && isTRUE(one_based)) sites$site_key <- colnames(geno)
      private$snp_sites <- sites
      private$snp_encoding <- enc %||% "dosage"
      private$was_pruned <- if (is.null(pruned)) NULL else isTRUE(pruned)
      private$allele_counted <- if (is.null(allele)) NULL else
        match.arg(allele, c("alt", "ref"))
      mat <- as.matrix(geno)
      if (is.null(samples)) samples <- rownames(mat)
      if (is.null(samples)) samples <- as.character(seq_len(nrow(mat)))
      rownames(mat) <- samples
      private$geno_mat <- mat
      private$sample_ids <- samples
      private$active_ids <- samples
      # The primary panel is named for what it is, so `$genotype("full")` finds it when the
      # object was built unpruned and nothing else has to know which call made it.
      private$primary <- if (isTRUE(pruned)) "pruned" else
        if (isFALSE(pruned)) "full" else "genotypes"
      private$panel_list <- stats::setNames(
        list(list(genotype = mat, allele = private$allele_counted,
                  pruned = private$was_pruned, encoding = private$snp_encoding,
                  sites = .align_sites(private$snp_sites, colnames(mat)))),
        private$primary)
      private$n_pcs <- n_pcs
      private$pca_panel <- pca_panel
      private$ps <- pop_structure(private$pca_matrix(), samples = samples, n_pcs = n_pcs,
                                  umap = FALSE)
      private$colors <- if (is.null(colors)) list() else colors
      if (!is.null(full)) self$add_panel("full", full, pruned = FALSE)
      if (!is.null(meta)) self$add_meta(meta)
      invisible(self)
    },

    #' @description Register another genotype panel under a name, so one object can hold the
    #'   pruned SNPs for PCA / admixture and the full set for the analyses where the
    #'   correlation between neighbouring SNPs is the signal.
    #' @param name Panel name (`"full"` and `"pruned"` are the ones other functions ask for).
    #' @param geno Genotype matrix, or a [load_genotypes()] list.
    #' @param allele,pruned What this panel is; taken from a `load_genotypes()` list when it
    #'   says so.
    add_panel = function(name, geno, allele = NULL, pruned = NULL, sites = NULL,
                         encoding = NULL) {
      private$ensure_panels()
      if (!is.character(name) || length(name) != 1 || !nzchar(name))
        stop("`name` must be a single non-empty string", call. = FALSE)
      if (is.list(geno) && !is.null(geno$genotype)) {
        if (is.null(allele)) allele <- geno$allele
        if (is.null(pruned)) pruned <- geno$pruned
        if (is.null(sites)) sites <- geno$sites
        if (is.null(encoding)) encoding <- geno$encoding
        if (is.null(rownames(geno$genotype)) && !is.null(geno$sample.id))
          rownames(geno$genotype) <- geno$sample.id
        geno <- geno$genotype
      }
      mat <- as.matrix(geno)
      if (is.null(rownames(mat)) || is.null(colnames(mat)))
        stop("a panel needs sample row names and `chr:pos` column names", call. = FALSE)
      gone <- setdiff(private$sample_ids, rownames(mat))
      if (length(gone))
        stop("panel \"", name, "\" is missing ", length(gone), " of this object's ",
             length(private$sample_ids), " samples", call. = FALSE)
      private$panel_list[[name]] <- list(
        genotype = mat[private$sample_ids, , drop = FALSE],
        allele = if (is.null(allele)) NULL else match.arg(allele, c("alt", "ref")),
        pruned = if (is.null(pruned)) NULL else isTRUE(pruned),
        encoding = encoding %||% "dosage",
        sites = .align_sites(sites, colnames(mat)))
      invisible(self)
    },

    #' @description The panels this object holds, primary first.
    panels = function() {
      private$ensure_panels()
      unique(c(private$primary, names(private$panel_list)))
    },

    #' @description Attach/replace metadata; auto-assigns colours for new columns.
    #' @param meta Data frame with a `sample` column.
    #' @param colors Optional colour overrides (`column -> (level -> colour)`).
    add_meta = function(meta, colors = NULL, colours = NULL) {
      colors <- .alias_arg("colors", "colours")
      meta <- .normalise_meta(meta)
      meta <- as.data.frame(meta, stringsAsFactors = FALSE)
      if (!"sample" %in% names(meta)) stop("`meta` needs a `sample` column", call. = FALSE)
      private$meta_df <- meta
      private$ps$meta <- meta
      auto <- meta_colors(meta)
      for (nm in names(auto)) if (is.null(private$colors[[nm]])) private$colors[[nm]] <- auto[[nm]]
      if (!is.null(colors)) self$set_colors(colors)
      invisible(self)
    },

    #' @description Set/override colour maps for metadata columns.
    #' @param colors Named list `column -> (level -> colour)`.
    set_colors = function(colors, colours = NULL) {
      colors <- .alias_arg("colors", "colours")
      for (nm in names(colors)) private$colors[[nm]] <- colors[[nm]]
      invisible(self)
    },

    #' @description Fix the level order of a metadata column. Because every plot reads
    #'   the shared metadata and colour map, this one order flows through the legends
    #'   (PCA / UMAP), the admixture facet order, and the colour strips. Existing colours
    #'   follow their level (only the order changes); unlisted levels are dropped.
    #' @param column Metadata column name.
    #' @param levels The desired level order.
    set_levels = function(column, levels) {
      m <- private$meta_df
      if (is.null(m) || !column %in% names(m))
        stop("no metadata column '", column, "'", call. = FALSE)
      m[[column]] <- factor(as.character(m[[column]]), levels = levels)
      private$meta_df <- m
      private$ps$meta <- m
      cur <- private$colors[[column]]
      private$colors[[column]] <-
        if (!is.null(cur) && all(levels %in% names(cur))) cur[levels]
        else meta_colors(m, cols = column)[[column]]
      invisible(self)
    },

    #' @description Compute a UMAP embedding.
    #' @param pca_components PCs feeding UMAP: a count (`>= 1`) or a variance fraction
    #'   (`0 < x < 1`, e.g. `0.1` uses enough PCs for 10% of variance, via
    #'   [n_pcs_for_variance()]).
    #' @param n_neighbors,min_dist UMAP parameters.
    #' @param seed Random seed.
    run_umap = function(pca_components = 30, n_neighbors = 15, min_dist = 0.1, seed = 42) {
      private$collapse_to_active()
      # UMAP runs on PCA scores, so it takes the same matrix the object's own PCA does --
      # not `geno_mat`, which on an allele-index panel would have a call of allele 2 read as
      # two copies of the alternate and separate samples by *which* alternate they carry as
      # though that were a distance.
      mat <- private$pca_matrix()
      private$ps <- pop_structure(mat, samples = rownames(mat),
                                  meta = private$meta_df, n_pcs = private$n_pcs,
                                  umap = TRUE, umap_pca = pca_components,
                                  n_neighbors = n_neighbors, min_dist = min_dist, seed = seed)
      invisible(self)
    },

    #' @description Fit sNMF admixture (cached and quiet by default).
    #' @param K,rep,alpha,seed,cpu,cache,cache_dir,verbose,log_file Passed to [run_snmf()].
    run_snmf = function(K = 1:10, rep = 10, alpha = 10, seed = 42, cpu = 1,
                        cache = TRUE, cache_dir = NULL, verbose = FALSE, log_file = NULL) {
      private$collapse_to_active()
      # sNMF's `.geno` alphabet is 0/1/2/9 -- one number per (sample, locus) counting copies
      # of one allele -- so a multiallelic locus has no representation in it. One-hot would
      # not rescue it either: sNMF would read the indicator columns as independent loci when
      # they are perfectly anti-correlated, inflating the locus count and distorting the
      # ancestry estimate. So admixture takes the biallelic subset, which the object derives
      # and says the size of. A genome-wide summary is the place that costs least.
      mat <- self$genotype(needs = "dosage")
      private$snmf_fit <- run_snmf(
        list(genotype = mat, sample.id = rownames(mat)),
        K = K, rep = rep, alpha = alpha, seed = seed, cpu = cpu,
        cache = cache, cache_dir = cache_dir, verbose = verbose, log_file = log_file)
      invisible(self)
    },

    #' @description Best K (cross-entropy) from the fitted sNMF.
    #' @param stat Combine replicates by `"min"` or `"mean"`; see [snmf_best_k()].
    best_k = function(stat = c("min", "mean")) snmf_best_k(private$snmf_fit, stat = match.arg(stat)),

    #' @description Per-K cross-entropy summary of the sNMF replicates
    #'   (see [snmf_cross_entropy()]).
    #' @param ... Passed to [snmf_cross_entropy()].
    cross_entropy = function(...) {
      private$require_snmf()
      snmf_cross_entropy(private$snmf_fit, ...)
    },

    #' @description Cross-entropy elbow plot for choosing K
    #'   (see [plot_snmf_cross_entropy()]).
    #' @param ... Passed to [plot_snmf_cross_entropy()].
    plot_cross_entropy = function(...) {
      private$require_snmf()
      plot_snmf_cross_entropy(private$snmf_fit, ...)
    },

    #' @description The fitted sNMF result (`NULL` before `run_snmf()`).
    get_snmf_fit = function() private$snmf_fit,

    #' @description One admixture plot per K as pages, for a multi-page PDF
    #'   (see [plot_admixture_multi_k()]).
    #' @param ... Passed to [plot_admixture_multi_k()].
    plot_admixture_multi_k = function(...) {
      private$require_snmf()
      plot_admixture_multi_k(self, ...)
    },

    #' @description Q (ancestry) matrix at K, restricted to the active samples.
    #' @param K Number of ancestral populations (default the best K).
    #' @param run Replicate (default the lowest cross-entropy run).
    q = function(K = NULL, run = NULL) {
      if (is.null(private$snmf_fit)) stop("run_snmf() first", call. = FALSE)
      if (is.null(K)) K <- self$best_k()
      qq <- snmf_q(private$snmf_fit, K = K, run = run)
      qq[rownames(qq) %in% private$active_ids, , drop = FALSE]
    },

    #' @description Restrict to a set of samples in place (used by `subset()`).
    #' @param samples Sample ids to keep.
    restrict = function(samples) {
      private$active_ids <- intersect(private$sample_ids, samples)
      invisible(self)
    },

    #' @description A new `PopStructure` limited to given samples and/or metadata
    #'   matches (does not recompute PCA/UMAP/sNMF -- the embeddings are shared and
    #'   simply filtered).
    #' @param samples Sample ids to keep.
    #' @param drop_samples Sample ids to remove, applied after `samples` and the metadata
    #'   filters. For taking out the one or two samples an otherwise good subset should not
    #'   include -- a contaminant, a duplicate -- without having to spell out the keep list.
    #'   Ids that are not in the object are reported rather than ignored, since a typo would
    #'   otherwise look exactly like a sample that was already gone.
    #' @param ... `column = value(s)` metadata filters (e.g. `region = "West Africa"`).
    subset = function(samples = NULL, drop_samples = NULL, ...) {
      keep <- private$active_ids
      if (!is.null(samples)) keep <- intersect(keep, samples)
      filt <- list(...)
      if (length(filt) && !is.null(private$meta_df)) {
        m <- private$meta_df
        ok <- rep(TRUE, nrow(m))
        for (nm in names(filt)) if (nm %in% names(m)) ok <- ok & (m[[nm]] %in% filt[[nm]])
        keep <- intersect(keep, m$sample[ok])
      }
      if (!is.null(drop_samples)) {
        drop_samples <- as.character(drop_samples)
        unknown <- setdiff(drop_samples, private$sample_ids)
        if (length(unknown))
          warning("drop_samples not in this object: ", paste(unknown, collapse = ", "),
                  call. = FALSE)
        keep <- setdiff(keep, drop_samples)
      }
      if (!length(keep)) stop("that leaves no samples", call. = FALSE)
      new <- self$clone(deep = FALSE)
      new$restrict(keep)
      new
    },

    #' @description The genotype matrix for the active samples (samples x SNPs).
    #' @param panel Panel to return by name; the primary one when `NULL`.
    #' @param prefer Panel to use *if the object has it*, falling back to the primary one --
    #'   how an analysis asks for the panel it wants without requiring it.
    #' @param needs What the caller can read: `"dosage"`, `"allele_index"`, or `"onehot"`
    #'   (one indicator column per ALT, which is how a multiallelic site reaches PCA). The
    #'   object serves or derives the panel that answers it, so one object covers every
    #'   analysis. `NULL` takes the panel as it is.
    genotype = function(panel = NULL, prefer = NULL, needs = NULL) {
      nm <- private$pick_panel(panel, prefer)
      if (!is.null(needs)) nm <- private$panel_for(nm, needs)   # may add a derived panel
      p <- private$panel_list[[nm]]
      p$genotype[private$active_ids, , drop = FALSE]
    },
    #' @description How a panel's `$genotype()` is coded: `"dosage"` (alternate copies) or
    #'   `"allele_index"` (which allele each sample carries). An object built from a bare
    #'   matrix reports `"dosage"`, since that is what one has always been.
    #' @param panel Which panel to report on; the primary one when `NULL`.
    encoding = function(panel = NULL) private$panel_field(panel, "encoding") %||% "dosage",
    #' @description PCA scores for the active samples.
    pca_scores = function() private$ps$pca[private$idx(), , drop = FALSE],
    #' @description PCA variance-explained table.
    pca_variance = function() private$ps$pca_var,
    #' @description The full [stats::prcomp()] object (all samples).
    prcomp = function() private$ps$prcomp,
    #' @description UMAP data frame for the active samples (or `NULL`).
    umap_df = function() if (is.null(private$ps$umap)) NULL else
      private$ps$umap[private$ps$umap$sample %in% private$active_ids, , drop = FALSE],
    #' @description Metadata for the active samples.
    get_meta = function() if (is.null(private$meta_df)) NULL else
      private$meta_df[private$meta_df$sample %in% private$active_ids, , drop = FALSE],
    #' @description The shared colour maps.
    get_colors = function() private$colors,

    #' @description Alias for `get_colors()`.
    get_colours = function() private$colors,
    #' @description Which allele the dosages count (`"alt"` / `"ref"`), or `NULL` when the
    #'   object does not record it (built from a bare matrix, or saved by an older version).
    #' @param panel Which panel to report on; the primary one when `NULL`.
    allele = function(panel = NULL) private$panel_field(panel, "allele"),
    #' @description The per-record allele table for a panel, or `NULL` when the panel was
    #'   built from a bare matrix that never carried one. One row per genotype column, in the
    #'   same order; see [load_genotypes()] for the columns. Ask it before reading a
    #'   per-allele number off a dosage matrix -- a column with `n_alt_real > 1` has had its
    #'   alternates collapsed onto one number and cannot answer one.
    #' @param panel Which panel to describe; the primary one when `NULL`.
    sites = function(panel = NULL) private$panel_field(panel, "sites"),
    #' @description Which convention the SNP positions follow, or `NULL` when the object does
    #'   not record it (built before this was tracked -- rebuild it, or pass `one_based`).
    positions = function() private$snp_positions,
    #' @description Whether the SNPs were LD-pruned (`TRUE` / `FALSE`), or `NULL` when the
    #'   object does not record it.
    #' @param panel Which panel to report on; the primary one when `NULL`.
    pruned = function(panel = NULL) private$panel_field(panel, "pruned"),
    #' @description Active sample ids.
    get_samples = function() private$active_ids,

    #' @description A `pop_structure` S3 view (active samples) for the `plot_*()` fns.
    as_ps = function() {
      idx <- private$idx()
      structure(list(samples = private$active_ids,
                     pca = private$ps$pca[idx, , drop = FALSE],
                     prcomp = private$ps$prcomp,
                     pca_var = private$ps$pca_var,
                     umap = self$umap_df(),
                     meta = self$get_meta()),
                class = "pop_structure")
    },

    #' @description PCA scatter coloured by a metadata column (shared colours).
    #' @param colour,color Metadata column to colour by (either spelling).
    #' @param pcs Which two PCs.
    #' @param ... Passed to [plot_pca()].
    plot_pca = function(colour = NULL, pcs = c(1, 2), ..., color = NULL) {
      colour <- .alias_arg("colour", "color")
      private$check_colour_column(colour, "plot_pca")
      plot_pca(self, pcs = pcs, colour = colour, ...)
    },

    #' @description UMAP scatter coloured by a metadata column (shared colours).
    #' @param colour,color Metadata column to colour by (either spelling).
    #' @param ... Passed to [plot_umap()].
    plot_umap = function(colour = NULL, ..., color = NULL) {
      colour <- .alias_arg("colour", "color")
      private$check_colour_column(colour, "plot_umap")
      plot_umap(self, colour = colour, ...)
    },

    #' @description Admixture bars; the group strip reuses the shared colour map, so it
    #'   matches the UMAP/PCA colouring.
    #' @param K Number of ancestral populations (default best K).
    #' @param group Metadata column to facet by and colour the strip with.
    #' @param colour,color Metadata column for the strip colours (default `group`, either
    #'   spelling). To recolour the ancestry clusters themselves, pass `colours` / `colors`
    #'   (a palette), which goes to [plot_admixture()].
    #' @param sample_order Explicit sample order (see [admixture_order()]).
    #' @param group_bar Draw the group colour strip.
    #' @param ... Passed to [plot_admixture()].
    plot_admixture = function(K = NULL, group = NULL, colour = group, sample_order = NULL,
                              group_bar = !is.null(group), ..., color = NULL) {
      colour <- .alias_arg("colour", "color")
      private$check_colour_column(colour, "plot_admixture")
      qq <- self$q(K)
      plot_admixture(qq, meta = self$get_meta(), group = group, sample_order = sample_order,
                     group_bar = group_bar,
                     group_colours = if (is.null(colour)) NULL else private$colors[[colour]], ...)
    },

    #' @description Combined UMAP + admixture figure (see [plot_structure_figure()]).
    #' @param group Metadata column to facet/colour the admixture by.
    #' @param ... Passed to [plot_structure_figure()].
    plot_figure = function(group = NULL, ...) plot_structure_figure(self, group = group, ...),

    #' @description Per-SNP population differentiation between the levels of a metadata
    #'   column (see [pop_diff()]); uses the object's genotype matrix (pass
    #'   `genotype = load_genotypes(vcf, prune = FALSE)` to run on the full unpruned set).
    #' @param group Metadata column defining the groups.
    #' @param ... Passed to [pop_diff()] (e.g. `statistic = "fst"`, `genotype = `).
    pop_diff = function(group = NULL, ...) pop_diff(self, group = group, ...),

    #' @description Per-SNP Jost's D between the levels of a metadata column ([jost_d()]).
    #' @param group Metadata column defining the groups.
    #' @param ... Passed to [jost_d()].
    jost_d = function(group = NULL, ...) jost_d(self, group = group, ...),

    #' @description Group-pair differentiation table across all statistics
    #'   (see [pop_diff_table()]).
    #' @param group Metadata column defining the groups.
    #' @param ... Passed to [pop_diff_table()].
    pop_diff_table = function(group = NULL, ...) pop_diff_table(self, group = group, ...),

    #' @description Group x group differentiation triangle heatmap (see
    #'   [plot_diff_heatmap()]); metadata annotation is resolved against this object's
    #'   metadata automatically.
    #' @param group Metadata column defining the groups.
    #' @param ... Passed to [plot_diff_heatmap()], plus `statistic` (`"jost_d"` default,
    #'   `"gst_hedrick"`, `"fst"`) selecting the measure, and `genotype` (a full/unpruned
    #'   override, see [pop_diff()]).
    plot_diff_heatmap = function(group = NULL, ...) {
      dots <- list(...)
      statistic <- if (is.null(dots$statistic)) "jost_d" else dots$statistic
      genotype  <- dots$genotype
      dots$statistic <- NULL; dots$genotype <- NULL
      # annotation strips read the shared colour map, as the UMAP, the PCA and the admixture
      # group bar do, so a level drawn here is the colour it is everywhere else
      if (!is.null(dots$annotate) &&
          is.null(dots$annotate_colours) && is.null(dots$annotate_colors) &&
          length(private$colors))
        dots$annotate_colours <- private$colors
      do.call(plot_diff_heatmap,
              c(list(pop_diff(self, group = group, statistic = statistic, genotype = genotype),
                     meta = private$meta_df), dots))
    },

    #' @description Alias of `plot_diff_heatmap()` for Jost's D.
    #' @param group Metadata column defining the groups.
    #' @param ... Passed to [plot_diff_heatmap()].
    plot_jost_d_heatmap = function(group = NULL, ...)
      plot_diff_heatmap(jost_d(self, group = group), meta = private$meta_df, ...),

    #' @description Per-SNP differentiation in long form (see [pop_diff_snps()]).
    #' @param group Metadata column defining the groups.
    #' @param ... Passed to [pop_diff()] (e.g. `statistic = `, `genotype = `).
    pop_diff_snps = function(group = NULL, ...) pop_diff_snps(pop_diff(self, group = group, ...)),

    #' @description Genome-wide differentiation Manhattan (see [plot_diff_manhattan()]).
    #' @param group Metadata column defining the groups.
    #' @param ... Passed to [plot_diff_manhattan()]; `statistic` / `genotype` go to [pop_diff()].
    plot_diff_manhattan = function(group = NULL, ...) {
      dots <- list(...)
      statistic <- if (is.null(dots$statistic)) "jost_d" else dots$statistic
      genotype  <- dots$genotype
      dots$statistic <- NULL; dots$genotype <- NULL
      do.call(plot_diff_manhattan,
              c(list(pop_diff(self, group = group, statistic = statistic, genotype = genotype)),
                dots))
    },

    #' @description Within-group diversity (see [pop_diversity()]).
    #' @param group Metadata column defining the groups.
    #' @param ... Passed to [pop_diversity()] (e.g. `by = `, `accessible = `).
    diversity = function(group = NULL, ...) pop_diversity(self, group = group, ...),

    #' @description Multilocus index of association (see [ld_index()]).
    #' @param group Metadata column defining the groups.
    #' @param ... Passed to [ld_index()].
    ld_index = function(group = NULL, ...) ld_index(self, group = group, ...),

    #' @description Beta scores for balancing selection (see [beta_score()]).
    #' @param group Metadata column defining the groups.
    #' @param ... Passed to [beta_score()].
    beta_score = function(group = NULL, ...) beta_score(self, group = group, ...),

    #' @description Phased haplotypes for a haplotype-homozygosity scan
    #'   (see [parasite_haplotypes()]).
    #' @param ... Passed to [parasite_haplotypes()] (e.g. `fws = `, `maf = `).
    haplotypes = function(...) parasite_haplotypes(self, ...),

    #' @description Genotype heatmap over one region, samples clustered
    #'   (see [plot_region_haplotypes()]).
    #' @param region The interval to draw.
    #' @param ... Passed to [plot_region_haplotypes()] (e.g. `split = `, `spacing = `).
    plot_region_haplotypes = function(region, ...)
      plot_region_haplotypes(self, region, ...),

    #' @description Integrated haplotype score (see [run_ihs()]); builds the haplotypes
    #'   first unless one is supplied.
    #' @param group Metadata column defining the groups.
    #' @param hap Optional [parasite_haplotypes()] object to reuse.
    #' @param ... Passed to [run_ihs()].
    ihs = function(group = NULL, hap = NULL, ...) {
      run_ihs(if (is.null(hap)) parasite_haplotypes(self) else hap, group = group, ...)
    },

    #' @description Save the whole workspace (genotype, PCA, UMAP, metadata, colours,
    #'   and sNMF fit) to an `.rds` file so it can be reloaded without recomputing.
    #'   Note: an sNMF fit references LEA project files on disk -- run `run_snmf()` with
    #'   a persistent `cache_dir` if you want the admixture to survive a reload.
    #' @param file Destination path.
    #' @param compress Passed to [saveRDS()] (default `"xz"` for a compact file).
    save = function(file, compress = "xz") {
      saveRDS(self, file, compress = compress)
      invisible(self)
    },

    #' @description Compact summary.
    #' @param ... Ignored.
    print = function(...) {
      cat("<PopStructure>", length(private$active_ids), "of",
          length(private$sample_ids), "samples,", ncol(private$ps$pca), "PCs\n")
      cat("  UMAP:", if (is.null(private$ps$umap)) "-" else "yes",
          "  sNMF:", if (is.null(private$snmf_fit)) "-"
                     else paste0("K ", paste(range(private$snmf_fit$K), collapse = "-")),
          "  meta:", if (is.null(private$meta_df)) "-"
                     else paste(setdiff(names(private$meta_df), "sample"), collapse = ", "),
          "\n")
      invisible(self)
    }
  ),
  private = list(
    pca_panel = "auto",

    # The matrix PCA and UMAP are actually run on.
    #
    # PCA needs a numeric matrix whose columns each count copies of ONE allele. An
    # allele-index panel is not that -- allele 2 is a different base, not two copies -- so
    # something has to be derived. There are two honest derivations and they answer different
    # questions:
    #
    #   "dosage"  drops every multiallelic site and keeps the biallelic ones. That is what
    #             the object did before, and on the pfpx1 384 kind of locus it drops exactly
    #             the site the analysis is about.
    #   "onehot"  keeps every site and gives each ALT its own indicator column. Nothing is
    #             dropped, and on a biallelic site the indicator column IS the alt-dosage
    #             column -- so a wholly biallelic panel's PCA does not move by one digit.
    #
    # "auto" takes one-hot whenever the panel holds a multiallelic site, because dropping a
    # site loses more than weighting one slightly heavily: a k-alternate site contributes k
    # anti-correlated columns, so it carries more of the total variance than a biallelic site
    # does. On a panel that is a handful of multiallelic sites in tens of thousands that is
    # unmeasurable; on a panel that is mostly multiallelic it is uniform. It is the middle
    # case that is uneven, and `pca_panel = "dosage"` is the way back.
    pca_matrix = function() {
      enc <- private$snp_encoding %||% "dosage"
      if (identical(enc, "dosage")) return(private$geno_mat)
      want <- private$pca_panel %||% "auto"
      if (identical(want, "auto")) {
        st <- private$panel_list[[private$primary]]$sites
        multi <- !is.null(st) && any(st$n_alt_real > 1L, na.rm = TRUE)
        want <- if (multi) "onehot" else "dosage"
      }
      # two statements on purpose: `panel_for()` *adds* the derived panel, and R evaluates
      # `private$panel_list` before the index expression that creates it, so the one-liner
      # reads the list as it was and comes back NULL
      nm <- private$panel_for(private$primary, want)
      private$panel_list[[nm]]$genotype
    },

    geno_mat = NULL, sample_ids = NULL, active_ids = NULL, meta_df = NULL,
    allele_counted = NULL, was_pruned = NULL, snp_positions = NULL, snp_sites = NULL,
    snp_encoding = NULL,
    panel_list = NULL, primary = NULL,
    told = character(0),

    # Resolve a panel name. `panel` is a requirement (absent -> error); `prefer` is a wish
    # (absent -> the primary panel, with one note when that means handing pruned SNPs to an
    # analysis that asked for the full set).
    # An object saved before panels existed has `geno_mat` but no panel list, so build one on
    # first use rather than letting every accessor fail on a loaded object.
    ensure_panels = function() {
      if (is.null(private$panel_list) || !length(private$panel_list)) {
        private$primary <- if (isTRUE(private$was_pruned)) "pruned" else
          if (isFALSE(private$was_pruned)) "full" else "genotypes"
        private$panel_list <- stats::setNames(
          list(list(genotype = private$geno_mat, allele = private$allele_counted,
                    sites = .align_sites(private$snp_sites, colnames(private$geno_mat)),
                    encoding = private$snp_encoding %||% "dosage",
                    pruned = private$was_pruned)), private$primary)
      }
      invisible(TRUE)
    },

    # The panel that can answer what an analysis needs, deriving one if that is possible.
    #
    # A dosage view of the biallelic sites is derivable from an allele-index panel: drop the
    # multiallelic columns, which a dosage cannot carry, and read index 0 as 0 copies and
    # index 1 as 2 (the haploid reading the rest of the package takes). The reverse is not
    # derivable -- a dosage never recorded which alternate a call named -- so an object built
    # from dosages says so rather than guessing.
    #
    # The derived view is cached as its own panel, so the second analysis that needs it does
    # not rebuild it and `$panels()` shows where it came from.
    panel_for = function(nm, needs) {
      needs <- match.arg(needs, c("dosage", "allele_index", "onehot"))
      have <- private$panel_list[[nm]]$encoding %||% "dosage"
      if (identical(have, needs)) return(nm)
      # One-hot is how a multiallelic site still reaches PCA and UMAP. One indicator column
      # per ALT, reference column dropped -- which on a biallelic site *is* the alt-dosage
      # column, so it is a strict generalisation and a biallelic panel's PCA does not move.
      # A dosage panel is therefore already one-hot and is handed back as it is.
      if (identical(needs, "onehot")) {
        if (identical(have, "dosage")) return(nm)
        return(private$onehot_panel(nm))
      }
      if (identical(needs, "allele_index")) {
        # A dosage *can* be read as an allele index where the site has one alternate: 0 copies
        # is allele 0 and 2 copies is allele 1. That is only safe when something says every
        # site is biallelic, because a dosage of 2 at a collapsed multiallelic site could be
        # any of the alternates -- which is the information a dosage never recorded.
        p <- private$panel_list[[nm]]
        if (!is.null(p$sites) && all(p$sites$n_alt_real == 1L, na.rm = TRUE) &&
            !anyNA(p$sites$n_alt_real))
          return(private$index_from_dosage(nm))
        stop("this analysis needs allele indices and panel \"", nm, "\" holds dosages. ",
             "Which alternate a call carries cannot be recovered from a dosage -- it was ",
             "never recorded. Rebuild with ",
             "`load_genotypes(vcf, variants = \"all\", encoding = \"allele_index\")`.",
             call. = FALSE)
      }

      derived <- paste0("biallelic_dosage", if (identical(nm, private$primary)) "" else
                        paste0("_", nm))
      if (!is.null(private$panel_list[[derived]])) return(derived)
      p <- private$panel_list[[nm]]
      sites <- p$sites
      if (is.null(sites))
        stop("panel \"", nm, "\" holds allele indices but carries no `sites` table, so ",
             "there is no way to tell which of its columns are biallelic. Rebuild it with ",
             "`load_genotypes()`.", call. = FALSE)
      keep <- which(sites$n_alt_real == 1L)
      if (!length(keep))
        stop("panel \"", nm, "\" has no biallelic site in it, so no dosage view exists.",
             call. = FALSE)
      g <- p$genotype[, keep, drop = FALSE]
      d <- matrix(NA_integer_, nrow(g), ncol(g), dimnames = dimnames(g))
      d[!is.na(g) & g == 0L] <- 0L
      d[!is.na(g) & g == 1L] <- 2L
      private$panel_list[[derived]] <- list(
        genotype = d, allele = "alt", pruned = p$pruned, encoding = "dosage",
        sites = sites[keep, , drop = FALSE])
      if (!derived %in% private$told) {
        private$told <- c(private$told, derived)
        message("derived a biallelic dosage panel (\"", derived, "\") from \"", nm,
                "\": ", length(keep), " of ", ncol(p$genotype), " sites, the ",
                ncol(p$genotype) - length(keep),
                " multiallelic ones left out because a dosage cannot carry them")
      }
      derived
    },

    # A biallelic dosage read as allele indices: 0 copies is allele 0, 2 copies is allele 1.
    # Only reachable when the panel's sites table says every site has one alternate.
    index_from_dosage = function(nm) {
      derived <- paste0("allele_index", if (identical(nm, private$primary)) "" else
                        paste0("_", nm))
      if (!is.null(private$panel_list[[derived]])) return(derived)
      p <- private$panel_list[[nm]]
      g <- p$genotype
      idx <- matrix(NA_integer_, nrow(g), ncol(g), dimnames = dimnames(g))
      idx[!is.na(g) & g == 0L] <- 0L
      idx[!is.na(g) & g == 2L] <- 1L
      private$panel_list[[derived]] <- list(
        genotype = idx, allele = p$allele, pruned = p$pruned,
        encoding = "allele_index", sites = p$sites)
      derived
    },

    # One column per (site, ALT), value 2 where the sample carries that allele and 0 where it
    # carries another, NA where it has no call. Values are 0/2 rather than 0/1 so the scale
    # matches the dosage matrix a biallelic panel would give.
    #
    # Worth knowing before reading a PCA off it: a site with k alternates contributes k
    # columns, so it carries more of the total variance than a biallelic one. That is a real
    # property of a more informative site rather than an artefact, but it means a panel with
    # a few very multiallelic loci is not weighted the way the same panel of biallelic loci
    # would be.
    onehot_panel = function(nm) {
      derived <- paste0("onehot", if (identical(nm, private$primary)) "" else paste0("_", nm))
      if (!is.null(private$panel_list[[derived]])) return(derived)
      p <- private$panel_list[[nm]]
      sites <- p$sites
      if (is.null(sites))
        stop("panel \"", nm, "\" holds allele indices but carries no `sites` table, so its ",
             "columns cannot be named by allele. Rebuild it with `load_genotypes()`.",
             call. = FALSE)
      g <- p$genotype
      cols <- list(); keys <- character(0); rows <- list()
      for (j in seq_len(ncol(g))) {
        alts <- sites$alt[[j]]
        real <- which(alts != "*")
        if (!length(real)) next
        for (k in real) {
          v <- rep(NA_integer_, nrow(g))
          ok <- !is.na(g[, j])
          v[ok] <- ifelse(g[ok, j] == k, 2L, 0L)
          cols[[length(cols) + 1L]] <- v
          # A site with one real ALT keeps its plain `chr:pos` name, because that column *is*
          # the alt-dosage column and everything downstream keys on `chr:pos`. Only a site
          # that needed splitting gets the widened `chr:pos:allele` key, which is the one
          # place a column name has to name an allele.
          key <- if (length(real) == 1L) sites$site_key[j]
                 else paste0(sites$site_key[j], ":", alts[k])
          keys <- c(keys, key)
          # the row describes THIS column, so its `site_key` has to be this column's name.
          # Leaving the parent `chr:pos` there makes the table look aligned while any join
          # on `site_key` silently merges a split site's alleles back together -- the exact
          # collapse the one-hot expansion exists to undo. `parent_key` keeps the link.
          r <- sites[j, , drop = FALSE]
          r$parent_key <- sites$site_key[j]
          r$allele <- alts[k]
          r$site_key <- key
          rows[[length(rows) + 1L]] <- r
        }
      }
      if (!length(cols)) stop("panel \"", nm, "\" has no real ALT allele to expand.",
                              call. = FALSE)
      m <- do.call(cbind, cols)
      dimnames(m) <- list(rownames(g), keys)
      private$panel_list[[derived]] <- list(
        genotype = m, allele = "alt", pruned = p$pruned, encoding = "dosage",
        sites = do.call(rbind, rows))
      if (!derived %in% private$told) {
        private$told <- c(private$told, derived)
        message("expanded \"", nm, "\" to a one-hot panel (\"", derived, "\"): ",
                ncol(g), " sites -> ", ncol(m), " allele columns, one per ALT. On a ",
                "biallelic site that is the alt-dosage column, so nothing biallelic moves.")
      }
      derived
    },

    # One field off one panel. Two statements on purpose: `pick_panel()` may *build*
    # `panel_list` (a legacy or freshly-loaded object has none), and R evaluates the
    # object of `[[` before the index, so `panel_list[[pick_panel()]]` reads the list as
    # it was *before* the build and silently returns NULL.
    panel_field = function(panel, field) {
      nm <- private$pick_panel(panel)
      private$panel_list[[nm]][[field]]
    },
    pick_panel = function(panel = NULL, prefer = NULL) {
      private$ensure_panels()
      if (!is.null(panel)) {
        if (!panel %in% names(private$panel_list))
          stop("no panel called \"", panel, "\"; this object has: ",
               paste(names(private$panel_list), collapse = ", "), call. = FALSE)
        return(panel)
      }
      if (!is.null(prefer) && prefer %in% names(private$panel_list)) return(prefer)
      if (identical(prefer, "full") &&
          isTRUE(private$panel_list[[private$primary]]$pruned) &&
          !"full" %in% private$told) {
        private$told <- c(private$told, "full")
        message("this reads best on the full SNP set, but the object only holds a pruned ",
                "panel; add one with ",
                "`$add_panel(\"full\", load_genotypes(vcf, prune = FALSE))`")
      }
      private$primary
    },
    colors = NULL, ps = NULL, snmf_fit = NULL, n_pcs = NULL,
    idx = function() match(private$active_ids, private$sample_ids),
    # `colour` names a metadata column whose shared colour map is reused; a palette passed
    # here reaches `private$colors[[colour]]` and fails on the subscript instead of on the
    # argument. The near-miss worth naming is `colours`, one letter away and a palette.
    check_colour_column = function(colour, fn) {
      if (is.null(colour)) return(invisible(NULL))
      cols <- if (is.null(private$meta_df)) character() else
        setdiff(names(private$meta_df), "sample")
      if (!is.character(colour) || length(colour) != 1L)
        stop("`colour` in ", fn, "() names one metadata column to colour by, not a palette.",
             if (identical(fn, "plot_admixture"))
               " To recolour the ancestry clusters, pass `colours = `." else "",
             "\n  columns available: ",
             if (length(cols)) paste(cols, collapse = ", ") else "none (add_meta() first)",
             call. = FALSE)
      if (!colour %in% cols)
        stop("no metadata column `", colour, "` to colour by.\n  columns available: ",
             if (length(cols)) paste(cols, collapse = ", ") else "none (add_meta() first)",
             call. = FALSE)
      invisible(NULL)
    },
    # Drop the samples a subset() excluded from the data itself, not just from the
    # accessors. Called before anything that refits: PCA, UMAP and sNMF are joint fits, so
    # a refit that still saw the excluded samples would place the kept ones by reference to
    # samples the caller asked to leave out.
    collapse_to_active = function() {
      if (identical(private$sample_ids, private$active_ids)) return(invisible(NULL))
      keep <- private$active_ids
      i <- private$idx()
      private$geno_mat <- private$geno_mat[i, , drop = FALSE]
      if (!is.null(private$panel_list))
        private$panel_list <- lapply(private$panel_list, function(p) {
          p$genotype <- p$genotype[keep, , drop = FALSE]
          p
        })
      if (!is.null(private$meta_df))
        private$meta_df <- private$meta_df[private$meta_df$sample %in% keep, , drop = FALSE]
      if (!is.null(private$ps)) {
        private$ps$pca <- private$ps$pca[i, , drop = FALSE]
        if (!is.null(private$ps$umap))
          private$ps$umap <-
            private$ps$umap[private$ps$umap$sample %in% keep, , drop = FALSE]
        if (!is.null(private$ps$prcomp$x))
          private$ps$prcomp$x <- private$ps$prcomp$x[i, , drop = FALSE]
        private$ps$meta <- private$meta_df
      }
      private$sample_ids <- keep
      invisible(NULL)
    },
    require_snmf = function() {
      if (is.null(private$snmf_fit)) stop("run_snmf() first", call. = FALSE)
      invisible(NULL)
    }
  )
)

# An R6 object carries its own copy of the methods it was built with, so one saved by an
# earlier version of the package comes back without anything added since. Move the saved state
# into an instance built by the class as it stands now; a field the saved object never had
# keeps the current default.
.refresh_pop_structure <- function(obj) {
  old <- obj$.__enclos_env__$private
  fresh <- PopStructure$new(matrix(c(0L, 1L, 1L, 0L), nrow = 2,
                                   dimnames = list(c("a", "b"), NULL)))
  new <- fresh$.__enclos_env__$private
  for (nm in names(PopStructure$private_fields)) {
    v <- old[[nm]]
    if (!is.null(v)) new[[nm]] <- v
  }
  # Only non-NULL fields are copied, so an object saved before panels existed would keep the
  # placeholder's panel list -- which describes the two-sample dummy, not the genotypes just
  # copied over. Drop it and let it be rebuilt from those on first use.
  if (is.null(old$panel_list)) { new$panel_list <- NULL; new$primary <- NULL }
  # Positions moved to 0-based at the one place they enter R. An object saved before that has
  # 1-based ids, which is a silent one-base offset against every IBD table and interval rather
  # than an error, so say so once.
  if (is.null(new$snp_positions))
    message("this object does not record whether its SNP positions are 0-based; objects saved ",
            "before that change carry 1-based ids, one base off every interval and IBD table. ",
            "Rebuild it, or reload the matrix with `PopStructure$new(..., one_based = TRUE)`")
  fresh
}

#' Load a saved PopStructure workspace
#'
#' Reads an `.rds` written by `PopStructure$save()` (or plain [saveRDS()]). The workspace is
#' re-bound to the installed version of the class, so a file written by an older version of
#' the package gains the methods added since.
#'
#' @param file Path to the `.rds` file.
#' @return The [PopStructure] object.
#' @examples
#' ps <- example_pop_structure(umap = FALSE)
#' f <- tempfile(fileext = ".rds")
#' ps$save(f)
#' load_pop_structure(f)
#' @export
load_pop_structure <- function(file) {
  obj <- readRDS(file)
  if (!inherits(obj, "PopStructure"))
    stop("`file` does not contain a PopStructure object", call. = FALSE)
  .refresh_pop_structure(obj)
}

# ---- example ---------------------------------------------------------------

#' Public population-structure example datasets
#'
#' Builds a [PopStructure] object from a bundled **public** genotype matrix. Two datasets
#' ship:
#' * `"ghana_cambodia"` (default) -- 60 Pf7 samples (30 Ghana, 30 Cambodia) at 49
#'   biallelic SNPs; two well-separated populations, a clean minimal demo.
#' * `"africa"` -- 258 published East/Central-African samples (`country`: DRC, Kenya,
#'   Tanzania, Uganda; finer `site` sub-regions; a macro `region`) at the 2,000 SNPs that
#'   most differentiate the sites (top Jost's D; see [top_differentiating_snps()]), so the
#'   regional structure is clear -- a richer, multi-region demo for the combined UMAP +
#'   admixture figure ([plot_structure_figure()]).
#'
#' @param dataset Which bundled dataset to load.
#' @param umap Also compute a UMAP embedding (needs \pkg{uwot}); skipped with a message
#'   if the package is missing.
#' @param seed Random seed for the UMAP.
#' @return A [PopStructure] object with a `meta` data frame.
#' @examples
#' ps <- example_pop_structure(umap = FALSE)
#' ps
#' @export
example_pop_structure <- function(dataset = c("ghana_cambodia", "africa", "multiallelic"),
                                  umap = TRUE, seed = 42) {
  dataset <- match.arg(dataset)
  file <- switch(dataset,
                 ghana_cambodia = "pop_structure_ghana_cambodia.rds",
                 africa = "pop_structure_africa.rds",
                 multiallelic = "pop_structure_multiallelic.rds")
  f <- system.file("extdata", file, package = "plasgenomicsutilsR")
  if (!nzchar(f)) stop("example genotype data not found in the installed package",
                       call. = FALSE)
  d <- readRDS(f)
  if (identical(dataset, "multiallelic")) {
    # The one fixture with multiallelic SNPs in it, and the only one built by a script
    # (`data-raw/pop_structure_multiallelic.R`). Its primary panel holds allele indices, so
    # a dosage analysis derives its own biallelic view on demand rather than being handed a
    # panel that quietly cannot carry the sites it was asked about.
    ps <- PopStructure$new(d$index, meta = d$meta)
    if (umap) message("the multiallelic fixture is for encoding tests; ",
                      "UMAP wants the dosage panel, so ask for it with ",
                      "`$genotype(needs = \"dosage\")`")
    return(ps)
  }
  # The other two fixtures are alt dosage (load_genotypes()'s default). `ghana_cambodia` also
  # carries a `full` panel: the same sparse genome-wide SNPs plus every biallelic SNP around
  # pfcrt / pfdhps / pfkelch13, so the locus, haplotype and EHH examples have real density
  # while PCA / UMAP / admixture keep reading the thinned set they were tuned on.
  ps <- PopStructure$new(list(genotype = d$genotype, allele = d$allele %||% "alt",
                              pruned = d$pruned, positions = d$positions %||% "0-based"),
                         meta = d$meta, full = d$full)
  if (umap) {
    if (requireNamespace("uwot", quietly = TRUE)) {
      # params that spread each dataset nicely out of the box
      if (dataset == "africa")
        ps$run_umap(pca_components = 0.8, n_neighbors = 25, min_dist = 0.4, seed = seed)
      else ps$run_umap(seed = seed)
    } else message("install 'uwot' to add a UMAP embedding to the example")
  }
  ps
}

#' Merge loaded genotype sets into one
#'
#' Column-binds the genotype matrices from two or more [load_genotypes()] results, keyed on
#' the `chr:pos` names the loader assigns.
#'
#' Merging here rather than upstream is the point. Two callers write different INFO and
#' FORMAT fields, and reconciling them so `bcftools merge` will accept both is real work
#' that changes nothing about the answer -- by the time a callset is a dosage matrix, all
#' that survives is the calls themselves, and those combine by position. So a variant that
#' had to be re-called separately, because the original callset arrived already filtered and
#' the caller's own output was gone, can be added without rebuilding the VCF it came from.
#'
#' What is checked, because these are the ways a merge is silently wrong:
#'
#' * **`allele`** must agree. One set counting alternate alleles and another counting
#'   reference alleles are indistinguishable after the fact, and mixing them inverts the
#'   dosages of whichever half disagrees.
#' * **Samples** must be the same set, in whatever order; the columns are realigned to the
#'   first set's order. `samples = "common"` intersects instead, saying how many it dropped.
#' * **Positions** must not collide. The same `chr:pos` from two callers is two answers to
#'   one question, and picking silently would hide the disagreement -- name which one wins
#'   with `on_overlap`.
#'
#' The result is sorted by chromosome and position, so it reads like a callset rather than
#' like the order the files happened to be given in. `snp.id` is renumbered: SNPRelate's ids
#' are per-file integers starting at 1, so they collide across files and mean nothing once
#' merged -- `chr:pos` is the identity that survives.
#'
#' @param ... Two or more [load_genotypes()] results.
#' @param samples `"identical"` (default) requires the same sample set in every input;
#'   `"common"` keeps the intersection.
#' @param on_overlap What to do when the same `chr:pos` appears in more than one input:
#'   `"error"` (default), or `"first"` / `"last"` to keep that input's calls.
#' @return A list shaped like [load_genotypes()]'s: `genotype`, `sample.id`, `snp.id`,
#'   `allele`, `pruned`, `positions`, `variants`. `pruned` is `TRUE` only when every input
#'   was pruned -- a SNP added to a pruned set never faced pruning itself, which is usually
#'   the reason for adding it.
#' @examples
#' \dontrun{
#' main <- load_genotypes("cohort.vcf.gz", gds = "cohort.gds", prune = TRUE)
#' extra <- load_genotypes("one_recalled_snp.vcf.gz", gds = "extra.gds", prune = FALSE)
#' both <- merge_genotypes(main, extra)
#' }
#' @seealso [load_genotypes()]
#' @export
merge_genotypes <- function(..., samples = c("identical", "common"),
                            on_overlap = c("error", "first", "last")) {
  sets <- list(...)
  samples <- match.arg(samples)
  on_overlap <- match.arg(on_overlap)
  if (length(sets) == 1 && is.list(sets[[1]]) && is.null(sets[[1]]$genotype))
    sets <- sets[[1]]                      # a list of sets, rather than several arguments
  if (length(sets) < 2) stop("give at least two genotype sets to merge", call. = FALSE)
  for (i in seq_along(sets)) {
    g <- sets[[i]]
    if (!is.list(g) || is.null(g$genotype) || is.null(g$sample.id))
      stop(sprintf("input %d is not a load_genotypes() result", i), call. = FALSE)
  }

  # `allele` says whether a 2 means two reference copies or two alternate ones, and nothing
  # in the matrix itself distinguishes them, so a mismatch here is unrecoverable later
  al <- unique(vapply(sets, function(g) g$allele %||% NA_character_, character(1)))
  if (length(al) > 1)
    stop("these sets count different alleles (", paste(al, collapse = ", "),
         "); reload them with the same `allele` before merging", call. = FALSE)
  pos <- unique(vapply(sets, function(g) g$positions %||% NA_character_, character(1)))
  if (length(pos) > 1)
    stop("these sets use different position conventions (", paste(pos, collapse = ", "), ")",
         call. = FALSE)
  vr <- unique(vapply(sets, function(g) g$variants %||% NA_character_, character(1)))
  if (length(vr) > 1)
    warning("merging sets loaded with different `variants` (", paste(vr, collapse = ", "),
            "); a `variants = \"all\"` matrix counts reference copies at multiallelic sites ",
            "and holds record types the other does not", call. = FALSE)

  keep <- Reduce(if (identical(samples, "common")) intersect else union,
                 lapply(sets, `[[`, "sample.id"))
  if (identical(samples, "identical")) {
    for (i in seq_along(sets)) {
      miss <- setdiff(keep, sets[[i]]$sample.id)
      if (length(miss))
        stop(sprintf(paste("input %d is missing %d of the %d samples (e.g. %s);",
                           "use samples = \"common\" to merge on the overlap"),
                     i, length(miss), length(keep),
                     paste(utils::head(miss, 3), collapse = ", ")), call. = FALSE)
    }
  } else {
    keep <- sets[[1]]$sample.id[sets[[1]]$sample.id %in% keep]   # first set's order
    if (!length(keep)) stop("the sets share no samples", call. = FALSE)
    dropped <- length(unique(unlist(lapply(sets, `[[`, "sample.id")))) - length(keep)
    if (dropped) message("merging on ", length(keep), " shared sample(s); ", dropped,
                         " not in every set were dropped")
  }
  keep <- sets[[1]]$sample.id[sets[[1]]$sample.id %in% keep]

  mats <- lapply(sets, function(g) {
    m <- g$genotype
    if (is.null(rownames(m))) rownames(m) <- g$sample.id
    m[keep, , drop = FALSE]
  })
  ids <- lapply(mats, colnames)
  dup <- unique(unlist(ids)[duplicated(unlist(ids))])
  if (length(dup)) {
    if (identical(on_overlap, "error"))
      stop(sprintf(paste("%d position(s) appear in more than one set (e.g. %s). Two callers",
                         "answering for one position is a disagreement, not a merge -- pass",
                         "on_overlap = \"first\" or \"last\" to say which set wins."),
                   length(dup), paste(utils::head(dup, 3), collapse = ", ")), call. = FALSE)
    order_seen <- if (identical(on_overlap, "first")) seq_along(mats) else rev(seq_along(mats))
    taken <- character(0)
    for (i in order_seen) {
      drop_here <- intersect(colnames(mats[[i]]), taken)
      if (length(drop_here))
        mats[[i]] <- mats[[i]][, setdiff(colnames(mats[[i]]), drop_here), drop = FALSE]
      taken <- c(taken, colnames(mats[[i]]))
    }
    message(length(dup), " overlapping position(s) resolved by keeping the ", on_overlap,
            " set")
  }

  out <- do.call(cbind, mats)
  loc <- .parse_snp_ids(colnames(out))
  # natural order where the chromosome is a number, alphabetical where it is not, so
  # chr2 sorts before chr10 rather than after it
  chr_num <- suppressWarnings(as.numeric(loc$chr))
  ord <- order(is.na(chr_num), chr_num, loc$chr, loc$pos)
  out <- out[, ord, drop = FALSE]

  list(genotype = out, sample.id = keep, snp.id = seq_len(ncol(out)),
       allele = sets[[1]]$allele, pruned = all(vapply(sets, function(g)
         isTRUE(g$pruned), logical(1))),
       positions = sets[[1]]$positions, variants = sets[[1]]$variants)
}
