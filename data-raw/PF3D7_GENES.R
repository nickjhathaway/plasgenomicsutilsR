# Build the bundled Pf3D7 gene-coordinate datasets from the VEuPathDB/PlasmoDB GFF.
#
# PF3D7_GENES        : every protein-coding gene, coordinates = the CDS span (min CDS
#                      start .. max CDS end over every isoform), i.e. the translated
#                      extent, excluding UTRs. Genes with no CDS feature fall back to the
#                      mRNA span and then to the `protein_coding_gene` bounds.
#                      Coordinates are 0-based half-open [start, end), matching BED and
#                      the rest of the package; the GFF is 1-based inclusive, so `start`
#                      is shifted down by one and `end` is carried over unchanged.
# PF_EXAMPLE_DRUG_GENES : the small drug-resistance / selection subset used in docs
#                      and examples.
#
# Source GFF: Pf3D7, release version 2020-09-01 (VEuPathDB / PlasmoDB). Re-run this
# script (adjusting `gff`) to regenerate the .rda files under data/.

gff <- "/Users/nhathaway/Documents/tank/data/genomes/plasmodium/genomes/pf/info/gff/Pf3D7.gff"
PF3D7_SOURCE_VERSION <- "2020-09-01"

g <- read.delim(gff, header = FALSE, comment.char = "#", quote = "",
                stringsAsFactors = FALSE,
                col.names = c("seqid", "source", "type", "start", "end",
                              "score", "strand", "phase", "attr"))

# length-preserving GFF attribute extractor (NA when the key is absent)
attr_get <- function(a, key) {
  has <- grepl(paste0("(^|;)", key, "="), a)
  val <- sub(paste0(".*(^|;)", key, "=([^;]*).*"), "\\2", a)
  ifelse(has, val, NA_character_)
}

gn <- g[g$type == "protein_coding_gene", ]
genes <- data.frame(
  gene_id     = attr_get(gn$attr, "ID"),
  Pf3D7_chrom = gn$seqid,
  Name        = attr_get(gn$attr, "Name"),
  gene_start  = gn$start,
  gene_end    = gn$end,
  stringsAsFactors = FALSE
)

# coordinates = CDS span (translated extent), collapsing every isoform to the outer bounds
cds <- g[g$type == "CDS", ]
cds_gid <- attr_get(cds$attr, "gene_id")                   # CDS rows carry gene_id directly
cstart <- tapply(cds$start, cds_gid, min)
cend   <- tapply(cds$end,   cds_gid, max)
genes$start <- as.integer(cstart[genes$gene_id])
genes$end   <- as.integer(cend[genes$gene_id])

# fall back to the mRNA span, then the gene bounds, for anything with no CDS feature
mr     <- g[g$type == "mRNA", ]
parent <- sub("\\.[0-9]+$", "", attr_get(mr$attr, "ID"))   # PF3D7_0709000.1 -> PF3D7_0709000
pstart <- tapply(mr$start, parent, min)
pend   <- tapply(mr$end,   parent, max)
no_cds <- is.na(genes$start)
genes$start[no_cds] <- as.integer(pstart[genes$gene_id[no_cds]])
genes$end[no_cds]   <- as.integer(pend[genes$gene_id[no_cds]])
miss <- is.na(genes$start)
genes$start[miss] <- genes$gene_start[miss]
genes$end[miss]   <- genes$gene_end[miss]
message(sum(no_cds), " genes had no CDS feature (fell back to mRNA / gene bounds)")

# GFF 1-based inclusive -> 0-based half-open [start, end)
genes$start <- as.integer(genes$start) - 1L
genes$end   <- as.integer(genes$end)

# display name: "pf" + lowercase(Name), or the gene_id when the gene has no Name
genes$name <- ifelse(is.na(genes$Name) | !nzchar(genes$Name),
                     genes$gene_id, paste0("pf", tolower(genes$Name)))
# surgical overrides: the folate genes carry long compound Names in the GFF
genes$name[genes$gene_id == "PF3D7_0417200"] <- "pfdhfr"   # Name = DHFR-TS
genes$name[genes$gene_id == "PF3D7_0810800"] <- "pfdhps"   # Name = PPPK-DHPS

# numeric chromosome (1..14); apicoplast / mitochondrion keep API / MIT
chrom <- sub("^Pf3D7_", "", sub("_v3$", "", genes$Pf3D7_chrom))
num   <- suppressWarnings(as.integer(chrom))
genes$chrom <- ifelse(is.na(num), chrom, as.character(num))

PF3D7_GENES <- genes[order(match(genes$chrom, c(1:14, "API", "MIT")), genes$start),
                     c("Pf3D7_chrom", "start", "end", "chrom", "gene_id", "name")]
rownames(PF3D7_GENES) <- NULL
attr(PF3D7_GENES, "source_version") <- PF3D7_SOURCE_VERSION

# the documented drug-resistance / selection example subset
.drug_ids <- c("PF3D7_0709000", "PF3D7_0417200", "PF3D7_0523000", "PF3D7_0810800",
               "PF3D7_1343700", "PF3D7_0629500", "PF3D7_1224000", "PF3D7_0720700")
PF_EXAMPLE_DRUG_GENES <- PF3D7_GENES[match(.drug_ids, PF3D7_GENES$gene_id), ]
rownames(PF_EXAMPLE_DRUG_GENES) <- NULL
attr(PF_EXAMPLE_DRUG_GENES, "source_version") <- PF3D7_SOURCE_VERSION

save(PF3D7_GENES, file = "data/PF3D7_GENES.rda", compress = "xz")
save(PF_EXAMPLE_DRUG_GENES, file = "data/PF_EXAMPLE_DRUG_GENES.rda", compress = "xz")
message(nrow(PF3D7_GENES), " genes -> data/PF3D7_GENES.rda (",
        file.size("data/PF3D7_GENES.rda"), " bytes)")
