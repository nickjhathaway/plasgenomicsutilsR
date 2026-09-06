# Ensembl Protists download paths. The FTP layout is regular enough to build a URL from a
# species name, but the assembly token in every filename is not: it tracks whatever assembly
# that release carried (Pf3D7 is `ASM276v2` up to release 60 and `GCA000002765v3` from 61), so
# it is looked up per release rather than assumed.

.ENSEMBL_PROTISTS_FTP <- "https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists"

#' Default Ensembl Protists release
#'
#' The release [ensembl_gff_url()], [ensembl_genome_url()] and [ensembl_species()] use when
#' none is given. The species table for this release ships with the package, so the default
#' needs no download.
#'
#' @export
ENSEMBL_PROTISTS_RELEASE <- 63L

# Species, assembly and collection for every Plasmodium genome of the default release, taken
# from that release's species_EnsemblProtists.txt. `collection` is "" for the genomes Ensembl
# gives a top-level directory and non-empty for those filed under a collection.
.ENSEMBL_PLASMODIUM <- local({
  v <- c(
    "plasmodium_berghei", "GCA900002375v2", "",
    "plasmodium_berghei_gca_900044335", "PbK173", "protists_alveolata1_collection",
    "plasmodium_berghei_gca_900088445", "PbNK65E", "protists_alveolata1_collection",
    "plasmodium_berghei_gca_900095585", "PbSP11RLL", "protists_alveolata1_collection",
    "plasmodium_berghei_gca_900095635", "PbSP11A", "protists_alveolata1_collection",
    "plasmodium_chabaudi", "PCHAS01", "",
    "plasmodium_chabaudi_adami_gca_900095565", "PchDK", "protists_alveolata1_collection",
    "plasmodium_chabaudi_chabaudi_gca_900095605", "PchCB", "protists_alveolata1_collection",
    "plasmodium_coatneyi_gca_001680005", "ASM168000v1", "protists_alveolata1_collection",
    "plasmodium_cynomolgi_strain_b_gca_000321355", "PcynB_1.0", "protists_alveolata1_collection",
    "plasmodium_falciparum", "GCA000002765v3", "",
    "plasmodium_falciparum_7g8_gca_000150435", "Plas_falc_7G8_V1", "protists_alveolata1_collection",
    "plasmodium_falciparum_camp_malaysia_gca_000521115", "Plas_falc_Malayan_Camp_V1", "protists_alveolata1_collection",
    "plasmodium_falciparum_dd2_gca_000149795", "ASM14979v1", "protists_alveolata1_collection",
    "plasmodium_falciparum_fch_4_gca_000521155", "Plas_falc_FCH_4_V2", "protists_alveolata1_collection",
    "plasmodium_falciparum_hb3_gca_000149665", "ASM14966v2", "protists_alveolata1_collection",
    "plasmodium_falciparum_igh_cr14_gca_000186055", "ASM18605v2", "protists_alveolata1_collection",
    "plasmodium_falciparum_malips096_e11_gca_000521035", "Plas_falc_MaliPS096_E11_V1", "protists_alveolata1_collection",
    "plasmodium_falciparum_nf135_5_c10_gca_000521075", "Plas_falc_NF135_5_C10_V1", "protists_alveolata1_collection",
    "plasmodium_falciparum_nf54_gca_000401695", "Plas_falc_NF4_V1", "protists_alveolata1_collection",
    "plasmodium_falciparum_nf54_gca_002831795", "ASM283179v1", "protists_alveolata1_collection",
    "plasmodium_falciparum_palo_alto_uganda_gca_000521095", "Plas_falc_Uganda_Palo-Alto_FUP_H_V1", "protists_alveolata1_collection",
    "plasmodium_falciparum_raj116_gca_000186025", "ASM18602v2", "protists_alveolata1_collection",
    "plasmodium_falciparum_santa_lucia_gca_000150455", "Plas_falc_Santa_Lucia_Salvador_I_V1", "protists_alveolata1_collection",
    "plasmodium_falciparum_tanzania_2000708__gca_000521055", "Plas_falc_02000708_V1", "protists_alveolata1_collection",
    "plasmodium_falciparum_ugt5_1_gca_000401715", "Plas_falc_5_1_V1", "protists_alveolata1_collection",
    "plasmodium_falciparum_vietnam_oak_knoll_fvo__gca_000521015", "Plas_falc_Vietnam_Oak-Knoll_FVO_V1", "protists_alveolata1_collection",
    "plasmodium_fragile_gca_000956335", "Plas_frag_nilgiri_V1", "protists_alveolata1_collection",
    "plasmodium_gaboni_gca_001602025", "ASM160202v1", "protists_alveolata1_collection",
    "plasmodium_gallinaceum_gca_900005855", "PGAL8A", "protists_alveolata1_collection",
    "plasmodium_gonderi_gca_002157705", "Pgonderi_assembly01", "protists_alveolata1_collection",
    "plasmodium_inui_san_antonio_1_gca_000524495", "Plas_inui_San_Antonio_1_V1", "protists_alveolata1_collection",
    "plasmodium_knowlesi", "GCA000006355v3", "",
    "plasmodium_knowlesi_gca_002140095", "PKNOHv1", "protists_alveolata1_collection",
    "plasmodium_knowlesi_strain_h_gca_900004885", "PKNA1-C.2", "protists_alveolata1_collection",
    "plasmodium_malariae_gca_900088575", "Pmal", "protists_alveolata1_collection",
    "plasmodium_malariae_gca_900090045", "PmUG01", "protists_alveolata1_collection",
    "plasmodium_ovale_curtisi_gca_900088555", "Poc1", "protists_alveolata1_collection",
    "plasmodium_ovale_curtisi_gca_900088565", "Poc2", "protists_alveolata1_collection",
    "plasmodium_ovale_gca_900090025", "PowCR01", "protists_alveolata1_collection",
    "plasmodium_ovale_wallikeri_gca_900088485", "Pow1", "protists_alveolata1_collection",
    "plasmodium_ovale_wallikeri_gca_900088545", "Pow2", "protists_alveolata1_collection",
    "plasmodium_reichenowi_gca_001601855", "ASM160185v1", "protists_alveolata1_collection",
    "plasmodium_reichenowi_gca_900097025", "PRG01", "protists_alveolata1_collection",
    "plasmodium_relictum_gca_900005765", "PRELSG", "protists_alveolata1_collection",
    "plasmodium_sp_drc_itaito_gca_900240055", "PGABG02", "protists_alveolata1_collection",
    "plasmodium_sp_gca900257145", "GCA900257145v2", "",
    "plasmodium_sp_gorilla_clade_g2_gca_900097015", "PADLG01", "protists_alveolata1_collection",
    "plasmodium_vinckei_petteri_gca_000524515", "Plas_vinc_pett_CR_V1", "protists_alveolata1_collection",
    "plasmodium_vinckei_vinckei_gca_000709005", "Plas_vinc_vinckei_V1", "protists_alveolata1_collection",
    "plasmodium_vivax_brazil_i_gca_000320645", "PV_Brazil_I_V1", "protists_alveolata1_collection",
    "plasmodium_vivax_gca900093555", "GCA900093555v2", "",
    "plasmodium_vivax_india_vii_gca_000320625", "PV_India_VII_V2", "protists_alveolata1_collection",
    "plasmodium_vivax_mauritania_i_gca_000320665", "PV_Mauritania_I_V1", "protists_alveolata1_collection",
    "plasmodium_vivax_north_korean_gca_000320685", "PV_N_Korean_V1", "protists_alveolata1_collection",
    "plasmodium_yoelii", "GCA900002385v2", "",
    "plasmodium_yoelii_17x_gca_000505035", "Plas_yoel_17X_V2", "protists_alveolata1_collection",
    "plasmodium_yoelii_gca_900002395", "PYYM01", "protists_alveolata1_collection",
    "plasmodium_yoelii_yoelii_gca_000003085", "ASM308v2", "protists_alveolata1_collection")
  m <- matrix(v, ncol = 3, byrow = TRUE,
              dimnames = list(NULL, c("species", "assembly", "collection")))
  data.frame(m, stringsAsFactors = FALSE)
})

# Two-letter short names, expanded to the species epithet and then resolved against the
# release like any other bare epithet.
.ENSEMBL_ALIASES <- c(pf = "falciparum", pv = "vivax", pk = "knowlesi", pb = "berghei",
                      pc = "chabaudi", py = "yoelii", pm = "malariae", po = "ovale")

# One session's worth of fetched release tables, so asking for several URLs of a non-default
# release costs one download rather than one per call.
.ensembl_cache <- new.env(parent = emptyenv())

# The Plasmodium rows of a release's species table: species, assembly and collection.
.ensembl_release_table <- function(release) {
  key <- as.character(release)
  if (!is.null(.ensembl_cache[[key]])) return(.ensembl_cache[[key]])
  if (identical(as.integer(release), ENSEMBL_PROTISTS_RELEASE)) {
    .ensembl_cache[[key]] <- .ENSEMBL_PLASMODIUM
    return(.ENSEMBL_PLASMODIUM)
  }
  url <- sprintf("%s/release-%d/species_EnsemblProtists.txt", .ENSEMBL_PROTISTS_FTP, release)
  txt <- tryCatch(readLines(url, warn = FALSE), error = function(e) NULL)
  if (is.null(txt))
    stop("could not read the species list for Ensembl Protists release ", release,
         " from ", url, " -- check the release number, or pass `assembly` (and `collection`) ",
         "to build the URL without a lookup", call. = FALSE)
  f <- utils::read.delim(text = txt, header = FALSE, comment.char = "#", quote = "",
                         stringsAsFactors = FALSE)
  keep <- grepl("^plasmodium", f[[2]])
  coll <- sub("_core_.*", "", f[[14]][keep])
  sp <- f[[2]][keep]
  out <- data.frame(species = sp, assembly = f[[5]][keep],
                    collection = ifelse(coll == sp, "", coll), stringsAsFactors = FALSE)
  out <- out[order(out$species), , drop = FALSE]
  rownames(out) <- NULL
  .ensembl_cache[[key]] <- out
  out
}

.ensembl_release <- function(release) {
  r <- suppressWarnings(as.integer(release))
  if (length(r) != 1 || is.na(r) || r <= 0)
    stop("`release` must be one positive Ensembl Protists release number, e.g. ",
         ENSEMBL_PROTISTS_RELEASE, call. = FALSE)
  r
}

# "P. falciparum", "Pf", "falciparum" and "plasmodium_falciparum" all name the same genome.
.ensembl_resolve <- function(species, release) {
  if (length(species) != 1 || !is.character(species) || is.na(species))
    stop("`species` must be one species name", call. = FALSE)
  tab <- .ensembl_release_table(release)
  raw <- tolower(trimws(species))
  # An Ensembl name is taken as given -- a few carry a doubled underscore that tidying the
  # punctuation of the friendlier spellings below would flatten.
  hit <- match(raw, tab$species)
  if (!is.na(hit)) return(tab[hit, , drop = FALSE])
  q <- gsub("[^a-z0-9]+", "_", sub("^p\\.?[ _]+", "plasmodium_", raw))
  q <- sub("_+$", "", q)
  if (!is.na(.ENSEMBL_ALIASES[q])) q <- paste0("plasmodium_", .ENSEMBL_ALIASES[q])
  if (!startsWith(q, "plasmodium_")) q <- paste0("plasmodium_", q)
  hit <- match(q, tab$species)
  if (!is.na(hit)) return(tab[hit, , drop = FALSE])
  # A bare epithet names the genome Ensembl treats as that species' reference -- the one it
  # gives a top-level directory rather than filing under a collection. Asked of the release,
  # since both the name and which genome holds that spot change between releases.
  near <- tab[startsWith(tab$species, q), , drop = FALSE]
  ref <- near[!nzchar(near$collection), , drop = FALSE]
  if (nrow(ref) == 1) return(ref)
  stop("no single Ensembl Protists genome named \"", species, "\" in release ", release,
       if (nrow(near))
         paste0(" -- did you mean ", paste0("\"", near$species, "\"", collapse = ", "), "?")
       else " -- see ensembl_species() for the names available",
       call. = FALSE)
}

# Ensembl capitalises the species slug to make every filename's stem.
.ensembl_stem <- function(species) sub("^(.)", "\\U\\1", species, perl = TRUE)

.ensembl_dir <- function(release, kind, row) {
  parts <- c(.ENSEMBL_PROTISTS_FTP, sprintf("release-%d", release), kind,
             if (nzchar(row$collection)) row$collection, row$species)
  paste(parts, collapse = "/")
}

.ensembl_row <- function(species, release, assembly, collection) {
  if (!is.null(assembly)) {
    if (length(assembly) != 1 || !is.character(assembly) || !nzchar(assembly))
      stop("`assembly` must be one assembly name, e.g. \"GCA000002765v3\"", call. = FALSE)
    if (is.null(collection)) collection <- ""
    return(data.frame(species = tolower(species), assembly = assembly,
                      collection = collection, stringsAsFactors = FALSE))
  }
  row <- .ensembl_resolve(species, release)
  if (!is.null(collection)) row$collection <- collection
  row
}

#' Plasmodium genomes in an Ensembl Protists release
#'
#' The names [ensembl_gff_url()] and [ensembl_genome_url()] accept, with the assembly each
#' carries in that release. The default release's table ships with the package; any other
#' release is read from the Ensembl FTP site once per session.
#'
#' @param release Ensembl Protists release number. Defaults to [ENSEMBL_PROTISTS_RELEASE].
#' @param pattern Optional regular expression, matched against the species name -- e.g.
#'   `"vivax"` for the *P. vivax* genomes.
#' @return A tibble of `species` (the name to pass on), `assembly` (the token that appears in
#'   the filenames) and `collection` (the FTP collection directory the genome sits in, `""`
#'   for the genomes Ensembl gives a top-level directory).
#' @seealso [ensembl_gff_url()], [ensembl_genome_url()]
#' @examples
#' ensembl_species()
#'
#' # the several assemblies of one species
#' ensembl_species(pattern = "ovale")
#' @export
ensembl_species <- function(release = ENSEMBL_PROTISTS_RELEASE, pattern = NULL) {
  tab <- .ensembl_release_table(.ensembl_release(release))
  if (!is.null(pattern)) tab <- tab[grepl(pattern, tab$species), , drop = FALSE]
  tibble::as_tibble(tab)
}

#' Ensembl Protists download URLs
#'
#' The GFF3 annotation and the genome FASTA of a *Plasmodium* genome, so a released
#' annotation can be handed straight to [read_gff_cds()] or [snp_aa_positions()] without
#' looking the path up by hand. Nothing is downloaded -- these return the URL.
#'
#' `species` may be the Ensembl name (`"plasmodium_falciparum"`), the species alone
#' (`"falciparum"`), `"P. falciparum"`, or a two-letter short name (`"pf"`, `"pv"`, `"pk"`,
#' `"pb"`, `"pc"`, `"py"`). The short forms name the genome Ensembl treats as that species'
#' reference -- the one Ensembl gives a top-level directory. Species with no such genome,
#' *P. malariae* and *P. ovale* among them, have several competing assemblies and so must be
#' named in full; the error lists them, and [ensembl_species()] lists them all.
#'
#' The assembly token in every Ensembl filename follows whatever assembly the release
#' carried, and does change (Pf3D7 is `ASM276v2` through release 60 and `GCA000002765v3`
#' from release 61), so it is looked up for the release asked for rather than assumed. Pass
#' `assembly` to skip that lookup, which also builds URLs for a genome outside the table.
#'
#' @param species Species name; see details.
#' @param release Ensembl Protists release number. Defaults to [ENSEMBL_PROTISTS_RELEASE].
#' @param assembly Assembly token to use, skipping the species lookup. Optional.
#' @param collection FTP collection directory, or `""` for none. Optional; only needed
#'   alongside `assembly` for a genome outside the table.
#' @param masking Repeat masking of the genome FASTA: `"none"` (the default), `"soft"`
#'   (repeats in lower case) or `"hard"` (repeats as `N`).
#' @return One URL, as a string.
#' @seealso [ensembl_species()], [read_gff_cds()], [snp_aa_positions()]
#' @examples
#' ensembl_gff_url("falciparum")
#' ensembl_genome_url("falciparum")
#'
#' # any of these name the same genome
#' ensembl_gff_url("pf")
#' ensembl_gff_url("P. falciparum")
#'
#' # other species, and the assembly of an earlier release
#' ensembl_gff_url("vivax")
#' ensembl_genome_url("plasmodium_malariae_gca_900090045", masking = "soft")
#'
#' \dontrun{
#' # a whole released annotation, read straight from the web
#' cds <- read_gff_cds(ensembl_gff_url("falciparum"))
#' snp_aa_positions(snps, cds, fasta = ensembl_genome_url("falciparum"))
#' }
#' @export
ensembl_gff_url <- function(species, release = ENSEMBL_PROTISTS_RELEASE, assembly = NULL,
                            collection = NULL) {
  release <- .ensembl_release(release)
  row <- .ensembl_row(species, release, assembly, collection)
  sprintf("%s/%s.%s.%d.gff3.gz", .ensembl_dir(release, "gff3", row),
          .ensembl_stem(row$species), row$assembly, release)
}

#' @rdname ensembl_gff_url
#' @export
ensembl_genome_url <- function(species, release = ENSEMBL_PROTISTS_RELEASE,
                               masking = c("none", "soft", "hard"), assembly = NULL,
                               collection = NULL) {
  masking <- match.arg(masking)
  release <- .ensembl_release(release)
  row <- .ensembl_row(species, release, assembly, collection)
  dna <- c(none = "dna", soft = "dna_sm", hard = "dna_rm")[[masking]]
  sprintf("%s/dna/%s.%s.%s.toplevel.fa.gz", .ensembl_dir(release, "fasta", row),
          .ensembl_stem(row$species), row$assembly, dna)
}
