#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 02_diversity_ordination.R
# The reproducibility "money figure": alpha diversity + Bray-Curtis PCoA of a
# wheat rhizosphere microbiome propagated across sequential cycles.
#
# The 4 replicates per group cluster tightly (reproducible biology) and the
# groups trace a clear assembly trajectory (soil -> wash -> seed-borne ->
# cycles -> recovered). Re-run this anywhere and you get the same result -
# that is the reproducibility Renku is built to guarantee.
#
# Input : outputs/phyloseq.rds        (from 01_build_phyloseq.R)
# Output: outputs/figure_diversity_ordination.png
#         outputs/alpha_diversity.csv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(patchwork)
})

out_dir <- "outputs"
ps <- readRDS(file.path(out_dir, "phyloseq.rds"))

# named palette straight from the authors' metadata, ordered by propagation stage
pal <- sample_data(ps) %>% data.frame() %>%
  distinct(Name, RGBcol, Order2) %>% arrange(Order2)
palette <- setNames(pal$RGBcol, pal$Name)

# ---- 1. Alpha diversity -----------------------------------------------------
alpha <- estimate_richness(ps, measures = c("Observed", "Shannon"))
alpha <- alpha %>% rownames_to_column("id_samples") %>%
  mutate(id_samples = gsub("^X", "", id_samples)) %>%
  left_join(data.frame(sample_data(ps)) %>% rownames_to_column("id_samples"),
            by = "id_samples")

write.csv(alpha[, c("id_samples", "Name", "Observed", "Shannon")],
          file.path(out_dir, "alpha_diversity.csv"), row.names = FALSE)

p_alpha <- ggplot(alpha, aes(Name, Shannon, fill = Name)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.12, size = 1.4, colour = "grey20") +
  scale_fill_manual(values = palette) +
  labs(title = "Alpha diversity collapses as the community is propagated",
       x = NULL, y = "Shannon diversity") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# ---- 2. Beta diversity: Bray-Curtis PCoA ------------------------------------
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
ord    <- ordinate(ps_rel, method = "PCoA", distance = "bray")

ord_df <- plot_ordination(ps_rel, ord, justDF = TRUE)
centroids <- ord_df %>% group_by(Name) %>%
  summarise(Axis.1 = mean(Axis.1), Axis.2 = mean(Axis.2),
            Order2 = first(Order2), .groups = "drop") %>% arrange(Order2)
ev <- round(100 * ord$values$Relative_eig[1:2], 1)

p_beta <- ggplot(ord_df, aes(Axis.1, Axis.2, colour = Name)) +
  geom_path(data = centroids, mapping = aes(Axis.1, Axis.2, group = 1),
            colour = "grey50", linetype = "dashed", inherit.aes = FALSE) +
  geom_point(size = 2.6, alpha = 0.9) +
  scale_colour_manual(values = palette) +
  labs(title = "Bray-Curtis PCoA: reproducible assembly trajectory",
       x = sprintf("PCo1 (%.1f%%)", ev[1]),
       y = sprintf("PCo2 (%.1f%%)", ev[2]), colour = NULL) +
  theme_bw(base_size = 11)

fig <- p_alpha + p_beta +
  plot_annotation(
    title = "Wheat rhizosphere microbiome propagation (Garrido-Sanz & Keel, UNIL - NCCR WP3)",
    caption = "Data: Zenodo 10.5281/zenodo.14514438 (CC-BY-4.0) | Analysis reproduced via Renku")

ggsave(file.path(out_dir, "figure_diversity_ordination.png"),
       fig, width = 13, height = 5.4, dpi = 150)
message("Saved outputs/figure_diversity_ordination.png")

# provenance for the record
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
