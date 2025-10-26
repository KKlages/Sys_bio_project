# ---- LOAD LIBRARIES ----
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)

# ---- CONFIGURATION ----
setwd("/mnt/ceph-hdd/projects/scc_ubmg_devries/KKlages/mrna_subproject/analysis")

# Load the Excel file with clade annotations
# Replace "your_excel_file.xlsx" with the actual path to your Excel file
excel_file <- "Genes_in_clades.xlsx"  # UPDATE THIS PATH
clade_data <- read_excel(excel_file)

# Create a lookup table for gene annotations
# Convert Excel gene names (with spaces) to match data file format (with underscores)
gene_to_clade <- data.frame()

# Extract gene IDs from each clade column
for (col_name in colnames(clade_data)) {
  if (col_name != "" && !is.na(col_name)) {
    genes_in_clade <- clade_data[[col_name]]
    genes_in_clade <- genes_in_clade[!is.na(genes_in_clade) & genes_in_clade != ""]
    
    if (length(genes_in_clade) > 0) {
      # Convert spaces to underscores to match data file format
      temp_df <- data.frame(
        gene_id = str_replace_all(genes_in_clade, "\\s+", "_"),
        clade = col_name,
        stringsAsFactors = FALSE
      )
      gene_to_clade <- rbind(gene_to_clade, temp_df)
    }
  }
}

# ---- PROCESS ALL SPECIES ----
subfolders <- list.dirs(path = ".", recursive = FALSE, full.names = TRUE)

for (sub in subfolders) {
  species <- basename(sub)
  
  # Check if gene_regulation_summary.txt exists
  summary_file <- file.path(sub, "gene_regulation_summary.txt")
  if (!file.exists(summary_file)) {
    cat("Skipping", species, "- no gene_regulation_summary.txt found\n")
    next
  }
  
  cat("Processing", species, "...\n")
  
  # Read the gene regulation summary
  gene_summary <- read.table(summary_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Get all comparison columns (those starting with "Direction_")
  direction_cols <- grep("^Direction_", colnames(gene_summary), value = TRUE)
  
  # Create long format data for each comparison
  final_results <- list()
  
  for (dir_col in direction_cols) {
    # Extract comparison name (remove "Direction_" prefix)
    comparison <- str_replace(dir_col, "^Direction_", "")
    
    # Get corresponding significance column
    sig_col <- paste0("Significant_", comparison)
    
    # Create data frame for this comparison
    comp_data <- gene_summary %>%
      select(gene_id, all_of(dir_col), all_of(sig_col)) %>%
      rename(Direction = all_of(dir_col), Significant = all_of(sig_col)) %>%
      mutate(Comparison = comparison) %>%
      # Only keep significant genes or those with actual direction calls
      filter(Direction != "Not_DE" & Direction != "NA" & !is.na(Direction))
    
    if (nrow(comp_data) > 0) {
      final_results[[comparison]] <- comp_data
    }
  }
  
  # Combine all comparisons
  if (length(final_results) > 0) {
    all_comparisons <- bind_rows(final_results)
    
    # Add clade annotations
    annotated_results <- all_comparisons %>%
      left_join(gene_to_clade, by = "gene_id") %>%
      mutate(
        Annotation = ifelse(is.na(clade), "SPL_gene", clade)
      ) %>%
      select(gene_id, Annotation, Direction, Comparison) %>%
      arrange(Annotation, gene_id, Comparison)
    
    # Save results
    output_file <- file.path(sub, "annotated_gene_directions.txt")
    write.table(annotated_results,
                file = output_file,
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    cat("  Saved annotated_gene_directions.txt for", species, "- found", nrow(annotated_results), "significant gene-comparison pairs\n")
    
    # Optional: Create a summary table showing counts by annotation and direction
    summary_table <- annotated_results %>%
      group_by(Annotation, Direction, Comparison) %>%
      summarise(Count = n(), .groups = "drop") %>%
      pivot_wider(names_from = c(Direction, Comparison), 
                  values_from = Count, 
                  values_fill = 0,
                  names_sep = "_")
    
    summary_output_file <- file.path(sub, "annotation_direction_summary.txt")
    write.table(summary_table,
                file = summary_output_file,
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    cat("  Saved annotation_direction_summary.txt for", species, "\n")
    
  } else {
    cat("  No significant differentially expressed genes found for", species, "\n")
  }
}

cat("All annotated gene direction analyses completed!\n")
