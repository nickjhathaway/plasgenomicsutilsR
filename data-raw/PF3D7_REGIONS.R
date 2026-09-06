# Build the bundled Pf3D7 genome-region datasets from the BED assets shipped with the
# companion Python package (plasgenomicsutils/src/plasgenomicsutils/assets). Re-run to
# regenerate data/PF3D7_CORE_REGIONS.rda and data/PF3D7_PARALOG_GENES.rda.

assets <- "../plasgenomicsutils/src/plasgenomicsutils/assets"   # sibling repo in the homebase

# Strand is not in the BED assets, so it comes from the same GFF PF3D7_GENES.R reads. Keep the
# two paths in step: a strand from one release against coordinates from another is the kind of
# disagreement nothing downstream would notice.
gff <- "/Users/nhathaway/Documents/tank/data/genomes/plasmodium/genomes/pf/info/gff/Pf3D7.gff"

norm_chrom <- function(seqid) {
  chrom <- sub("^Pf3D7_", "", sub("_v3$", "", seqid))
  num <- suppressWarnings(as.integer(chrom))
  ifelse(is.na(num), chrom, as.character(num))
}

# core regions: 3-column BED (chr, start, end); everything else on a core chromosome is
# subtelomeric / internally hypervariable
core <- read.delim(file.path(assets, "pf3d7_core_regions.bed"), header = FALSE,
                   comment.char = "#", stringsAsFactors = FALSE,
                   col.names = c("Pf3D7_chrom", "start", "end"))
core$chrom <- norm_chrom(core$Pf3D7_chrom)
PF3D7_CORE_REGIONS <- core[order(match(core$chrom, c(1:14, "API", "MIT")), core$start),
                          c("Pf3D7_chrom", "start", "end", "chrom")]
rownames(PF3D7_CORE_REGIONS) <- NULL

# paralog / hypervariable gene families: 4-column BED, 4th field "gene_id|description"
par <- read.delim(file.path(assets, "pf3d7_paralog_genes.bed"), header = FALSE,
                  comment.char = "#", quote = "", stringsAsFactors = FALSE,
                  col.names = c("Pf3D7_chrom", "start", "end", "info"))
par$gene_id     <- sub("\\|.*$", "", par$info)
par$description <- sub("^[^|]*\\|", "", par$info)
par$chrom <- norm_chrom(par$Pf3D7_chrom)
# strand, joined on gene_id. A quarter of these families are pseudogenes -- var and rifin
# pseudogenes especially -- so the lookup has to span every gene-level feature type in the GFF
# and not just `protein_coding_gene`, or those come back NA.
g <- read.delim(gff, header = FALSE, comment.char = "#", quote = "",
                stringsAsFactors = FALSE,
                col.names = c("seqid", "source", "type", "start", "end",
                              "score", "strand", "phase", "attr"))
gene_level <- g[g$type %in% c("protein_coding_gene", "pseudogene", "ncRNA_gene"), ]
gene_level$ID <- sub(".*(^|;)ID=([^;]*).*", "\\2", gene_level$attr)
par$strand <- gene_level$strand[match(par$gene_id, gene_level$ID)]
stopifnot(all(par$strand %in% c("+", "-")))

PF3D7_PARALOG_GENES <- par[order(match(par$chrom, c(1:14, "API", "MIT")), par$start),
                          c("Pf3D7_chrom", "start", "end", "chrom", "gene_id", "description",
                            "strand")]
rownames(PF3D7_PARALOG_GENES) <- NULL

save(PF3D7_CORE_REGIONS,  file = "data/PF3D7_CORE_REGIONS.rda",  compress = "xz")
save(PF3D7_PARALOG_GENES, file = "data/PF3D7_PARALOG_GENES.rda", compress = "xz")
message(nrow(PF3D7_CORE_REGIONS), " core regions, ", nrow(PF3D7_PARALOG_GENES), " paralog genes")
