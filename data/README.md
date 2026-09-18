Run `Rscript analysis/00_get_data.R` from the project root (or `snakemake --cores 1`, which
runs it automatically). It unzips the mounted Zenodo data connector
(`assessing-rhizosphere-community-...-zenodo.14514438/`, a sibling of this project) and stages
the three files the pipeline reads, renamed to match:

| Archive file (`Amplicon_sequence_analyses/`) | Staged as             |
|-----------------------------------------------|------------------------|
| `ASV_sequences.txt`                            | `data/asv_counts.txt`  |
| `Taxtable_dada2.txt`                           | `data/taxonomy.txt`    |
| `metadata.txt`                                 | `data/metadata.txt`    |
