# Build the bundled Pf3D7 gene-coordinate datasets from the VEuPathDB/PlasmoDB GFF.
#
# PF3D7_GENES        : every protein-coding gene, coordinates = the mRNA (transcript)
#                      span (min transcript start .. max transcript end), i.e. CDS+UTR
#                      transcript extent rather than the wider `gene` feature.
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

# coordinates = transcript (mRNA) span, collapsing multiple isoforms to the outer bounds
mr     <- g[g$type == "mRNA", ]
parent <- sub("\\.[0-9]+$", "", attr_get(mr$attr, "ID"))   # PF3D7_0709000.1 -> PF3D7_0709000
pstart <- tapply(mr$start, parent, min)
pend   <- tapply(mr$end,   parent, max)
genes$start <- as.integer(pstart[genes$gene_id])
genes$end   <- as.integer(pend[genes$gene_id])
miss <- is.na(genes$start)                                 # genes with no mRNA: use gene bounds
genes$start[miss] <- genes$gene_start[miss]
genes$end[miss]   <- genes$gene_end[miss]

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
