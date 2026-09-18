#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 01_build_phyloseq.R
# Build a filtered phyloseq object from the wheat-rhizosphere 16S amplicon data.
#
# Data: Garrido-Sanz & Keel, "Sequential propagation of a reproducible wheat
#       rhizosphere microbiome" (Keel lab, UNIL).
#       Zenodo: https://doi.org/10.5281/zenodo.14514438  (CC-BY-4.0)
#
# Inputs  (in data/):
#   asv_counts.txt   samples x ASVs count table, ";"-separated, 1st col = sample id
#   taxonomy.txt     ASV x taxonomic-rank table, ";"-separated, 1st col = ASV sequence
#   metadata.txt     sample metadata, TAB-separated, key column = id_samples
# Output (in outputs/):
#   phyloseq.rds     a filtered phyloseq object ready for analysis
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(phyloseq)
  library(data.table)
  library(tibble)
})

data_dir <- "data"
out_dir  <- "outputs"
dir.create(out_dir, showWarnings = FALSE)

message("== Loading count table ==")
seqtab <- fread(file.path(data_dir, "asv_counts.txt"), sep = ";", header = TRUE)
seqtab <- column_to_rownames(as.data.frame(seqtab), "V1")
rownames(seqtab) <- gsub("_$", "", rownames(seqtab))   # tidy sample ids
message(sprintf("   %d samples x %d ASVs", nrow(seqtab), ncol(seqtab)))

message("== Loading taxonomy ==")
taxa <- read.table(file.path(data_dir, "taxonomy.txt"), sep = ";", header = TRUE)
taxa <- column_to_rownames(as.data.frame(taxa), "X")

message("== Loading sample metadata ==")
meta <- read.table(file.path(data_dir, "metadata.txt"), sep = "\t", header = TRUE)
rownames(meta) <- meta$id_samples
# fix the propagation order and drop stray quotes on the colour column
meta$Name   <- factor(meta$Name, levels = meta$Name[order(meta$Order2)][!duplicated(meta$Name[order(meta$Order2)])])
meta$RGBcol <- gsub('"', "", meta$RGBcol)

# keep only ASVs present in both count and taxonomy tables
common_asv <- intersect(colnames(seqtab), rownames(taxa))
seqtab <- seqtab[, common_asv]
taxa   <- taxa[common_asv, ]

message("== Assembling phyloseq object ==")
ps <- phyloseq(
  otu_table(as.matrix(seqtab), taxa_are_rows = FALSE),
  tax_table(as.matrix(taxa)),
  sample_data(meta)
)

# ---- Prevalence / abundance filter -----------------------------------------
# Keep ASVs seen in >= 2 samples with a total count >= 10 (removes noise/singletons)
prev  <- apply(otu_table(ps), 2, function(x) sum(x > 0))
total <- taxa_sums(ps)
ps <- prune_taxa(prev >= 2 & total >= 10, ps)

message(sprintf("== Filtered phyloseq: %d samples, %d ASVs ==",
                nsamples(ps), ntaxa(ps)))
saveRDS(ps, file.path(out_dir, "phyloseq.rds"))
message("Saved outputs/phyloseq.rds")
