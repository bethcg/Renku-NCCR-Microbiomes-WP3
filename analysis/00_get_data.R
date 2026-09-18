#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 00_get_data.R
# Stage the raw wheat-rhizosphere amplicon data into data/, from the
# read-only Zenodo data connector mounted as a sibling of this project.
#
# Source: Garrido-Sanz & Keel, "Sequential propagation of a reproducible
#         wheat rhizosphere microbiome" (Keel lab, UNIL).
#         Zenodo: https://doi.org/10.5281/zenodo.14514438  (CC-BY-4.0)
#
# The archive's raw file names don't match what the rest of the pipeline
# expects, so this script also renames them on the way into data/:
#   ASV_sequences.txt  -> data/asv_counts.txt
#   Taxtable_dada2.txt -> data/taxonomy.txt
#   metadata.txt       -> data/metadata.txt
# ---------------------------------------------------------------------------

data_dir <- "data"
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

rename_map <- c(
  "ASV_sequences.txt"  = "asv_counts.txt",
  "Taxtable_dada2.txt" = "taxonomy.txt",
  "metadata.txt"       = "metadata.txt"
)

if (all(file.exists(file.path(data_dir, rename_map)))) {
  message("== data/ already populated, nothing to do ==")
  quit(save = "no", status = 0)
}

# The connector is mounted next to the project root, e.g.
# ../assessing-rhizosphere-community-doi-10.5281-zenodo.14514438/
connector <- Sys.glob("../*rhizosphere*zenodo*14514438*")
if (length(connector) == 0) {
  stop("Could not find the Zenodo data connector (expected a directory next ",
       "to this project matching '*rhizosphere*zenodo*14514438*'). Mount the ",
       "'assessing-rhizosphere-community-...-zenodo.14514438' data connector ",
       "alongside this project, or place the files by hand in data/: ",
       paste(rename_map, collapse = ", "), ".")
}

zip_path <- Sys.glob(file.path(connector[1], "*", "RhizCom*.zip"))
if (length(zip_path) == 0) {
  stop(sprintf("Found the data connector at '%s' but no RhizCom*.zip inside it.",
               connector[1]))
}
zip_path <- zip_path[1]

message(sprintf("== Unzipping %s ==", zip_path))
extract_dir <- tempfile("rhizcom_")
utils::unzip(zip_path, exdir = extract_dir)

amp_dir <- Sys.glob(file.path(extract_dir, "*", "Amplicon_sequence_analyses"))
if (length(amp_dir) == 0) {
  stop("Unzipped archive did not contain an Amplicon_sequence_analyses/ folder.")
}
amp_dir <- amp_dir[1]

for (src in names(rename_map)) {
  ok <- file.copy(file.path(amp_dir, src), file.path(data_dir, rename_map[[src]]),
                   overwrite = TRUE)
  if (!ok) stop(sprintf("Failed to stage %s into %s/", src, data_dir))
}
unlink(extract_dir, recursive = TRUE)

message(sprintf("== Staged %s into %s/ ==",
                 paste(rename_map, collapse = ", "), data_dir))
