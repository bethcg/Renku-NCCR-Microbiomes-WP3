# Reproducible wheat rhizosphere microbiome (NCCR Microbiomes)

A small, real, end-to-end R project for a 15-minute Renku demo aimed at reproducing the core result of a
wheat-rhizosphere propagation experiment and showing how Renku makes that result
reproducible and shareable: **code + data + environment + compute in one place.**

> **The story in one line:** the biology here is *reproducible* (4 replicates
> per group cluster tightly, communities assemble the same way each cycle) — and
> Renku makes the *analysis of it* just as reproducible. One click, same figure,
> anywhere.

## What it does

1. `analysis/01_build_phyloseq.R`: builds a filtered `phyloseq` object from the
   ASV count table, taxonomy and sample metadata.
2. `analysis/02_diversity_ordination.R`: computes alpha diversity and a
   Bray-Curtis PCoA and writes the figure below to `outputs/`.

![expected output](outputs/expected_figure.png)

*Left:* Shannon diversity collapses from diverse soil, through a seed-borne
bottleneck, to a stable propagated community. *Right:* Bray-Curtis PCoA — tight
replicate clusters (reproducibility) tracing the assembly trajectory
(PCo1 ≈ 41%, PCo2 ≈ 14%).

## Data

Garrido-Sanz & Keel, *Sequential propagation of a reproducible wheat rhizosphere
microbiome* (Keel lab, UNIL).
Zenodo **10.5281/zenodo.14514438** · CC-BY-4.0. See `data/README.md`.

## Environment (conda)

`environment.yml` is the single source of truth for dependencies (R + `phyloseq`,
`vegan`, `tidyverse`, …). It is used to build the Renku session image
and to run the pipeline, so the environment is identical on a
laptop, in CI, and on RenkuLab.

## Run it

**On RenkuLab (the demo path):** open the project → start the RStudio session →
in the R console:

```r
source("analysis/01_build_phyloseq.R")
source("analysis/02_diversity_ordination.R")
```

The figure appears in `outputs/`.

**As a pipeline (optional):** `snakemake --cores 1` runs the whole DAG and only
re-runs a step when its inputs change.

**Locally:** `conda env create -f environment.yml && conda activate renku-wp3-demo`,
then the same two `source()` calls (or `snakemake`).

## Layout

```
├── analysis/
│   ├── 00_get_data.R              # provenance + re-download (not needed; data ships)
│   ├── 01_build_phyloseq.R        # data -> filtered phyloseq object
│   └── 02_diversity_ordination.R  # alpha + beta diversity -> figure
├── outputs/                       # figures & tables land here
├── environment.yml                # conda dependency spec (source of truth)
├── Dockerfile                     # custom Renku session image (RStudio + conda)
├── Snakefile                      # optional reproducible pipeline
├── DEMO_SCRIPT.md                 # the 15-minute run-of-show
└── README.md
```

## Attribution

Data © the original authors, reused under CC-BY-4.0 (Zenodo 10.5281/zenodo.14514438).
This demo scaffold is provided for teaching Renku within NCCR Microbiomes.
