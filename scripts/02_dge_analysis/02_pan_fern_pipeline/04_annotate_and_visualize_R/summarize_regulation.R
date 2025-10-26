# ---- LOAD LIBRARIES ----
library(dplyr)
library(tidyr)
library(stringr)

# ---- CONFIGURATION ----
setwd("/mnt/ceph-hdd/projects/scc_ubmg_devries/KKlages/mrna_subproject/analysis")

# ---- PROCESS ALL SPECIES ----
subfolders <- list.dirs(path = ".", recursive = FALSE, full.names = TRUE)

for (sub in subfolders) {
  species <- basename(sub)
  result_files <- list.files(sub, pattern = "^annotated_results_.*\\.txt$", full.names = TRUE)
  if (length(result_files) == 0) next
  
  cat("Summarizing gene directions for", species, "...\n")
  
  all_gene_data <- lapply(result_files, function(file) {
    comp <- str_match(basename(file), "^annotated_results_(.+)\\.txt$")[,2]
    df <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    df %>%
      mutate(
        Comparison = comp,
        Direction = case_when(
          FDR >= 0.05 ~ "Not_DE",
          logFC > 0 ~ "Upregulated",
          logFC < 0 ~ "Downregulated"
        ),
        Significant = FDR < 0.05
      ) %>%
      select(gene_id, Comparison, Direction, Significant)
  })
  
  gene_summary <- bind_rows(all_gene_data) %>%
    pivot_wider(names_from = Comparison, values_from = c(Direction, Significant),
                names_sep = "_", values_fill = list(Direction = "NA", Significant = NA))
  
  # Save the summary
  write.table(gene_summary,
              file = file.path(sub, "gene_regulation_summary.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat("  Saved gene_regulation_summary.txt for", species, "\n")
}

cat("All gene regulation summaries created!\n")

library(dplyr)
library(tidyr)
library(stringr)

setwd("/mnt/ceph-hdd/projects/scc_ubmg_devries/KKlages/mrna_subproject/analysis")

# Get all annotated_gene_directions.txt files in subfolders
files <- list.files(path = ".", pattern = "^annotated_gene_directions.txt$", 
                   recursive = TRUE, full.names = TRUE)

for (f in files) {
  species <- basename(dirname(f))
  df <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Split Comparison into tissue1 and tissue2
  df <- df %>%
    mutate(
      tissue1 = str_extract(Comparison, "^[^_]+"),
      tissue2 = str_extract(Comparison, "(?<=_vs_).+")
    )
  
  # Upregulated matrix
  up_mat <- df %>%
    filter(Direction == "Upregulated") %>%
    group_by(tissue1, tissue2) %>%
    summarise(Upregulated = n_distinct(gene_id), .groups = "drop") %>%
    pivot_wider(names_from = tissue2, values_from = Upregulated, values_fill = 0)
  
  # Downregulated matrix
  down_mat <- df %>%
    filter(Direction == "Downregulated") %>%
    group_by(tissue1, tissue2) %>%
    summarise(Downregulated = n_distinct(gene_id), .groups = "drop") %>%
    pivot_wider(names_from = tissue2, values_from = Downregulated, values_fill = 0)
  
  # Write outputs
  write.table(up_mat, file = file.path(dirname(f), paste0(species, "_upregulated_matrix.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(down_mat, file = file.path(dirname(f), paste0(species, "_downregulated_matrix.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  cat("Matrices saved for", species, "\n")
}

