# Optional Snakemake pipeline wiring the analysis into a reproducible DAG.
# Run the whole thing non-interactively (e.g. in a Renku session terminal or CI):
#     snakemake --cores 1
# Snakemake only re-runs a step when its inputs change, and the DAG is the
# machine-readable record of how every output was produced.

rule all:
    input:
        "outputs/figure_diversity_ordination.png",
        "outputs/alpha_diversity.csv"

rule build_phyloseq:
    input:
        counts   = "data/asv_counts.txt",
        taxonomy = "data/taxonomy.txt",
        metadata = "data/metadata.txt"
    output:
        "outputs/phyloseq.rds"
    shell:
        "Rscript analysis/01_build_phyloseq.R"

rule diversity_ordination:
    input:
        "outputs/phyloseq.rds"
    output:
        "outputs/figure_diversity_ordination.png",
        "outputs/alpha_diversity.csv"
    shell:
        "Rscript analysis/02_diversity_ordination.R"
