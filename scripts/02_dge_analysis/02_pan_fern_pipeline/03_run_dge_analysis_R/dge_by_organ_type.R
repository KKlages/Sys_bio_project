library(edgeR)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(readxl)

# --- Load and clean organ mapping once ---
organ_mapping <- read_excel("organ_type.xlsx")

# Clean the organ names to match column names in counts
organ_mapping$Clean_Tissue <- gsub("[^a-z0-9]+", "_", tolower(organ_mapping$`Organ name`))

manual_add <- data.frame(
  `Organ type` = c(
    "Roots", "Leaves", "Leaves", "Roots", "Leaves", "Roots",
    "Reproductive Organ", "Leaves", "Leaves", "Reproductive Organ",
    "Stem", "Stem", "Reproductive Organ",
    "Reproductive Organ", "Roots", "Leaves", "Roots",
    "Leaves", "Reproductive Organ", "Leaves"  # ← new entries for Sam
  ),
  Clean_Tissue = c(
    "root", "sterile_leaflet", "fertile_leaflet", 
    "mature_root", "mature_sterile_leaflet", "young_root",
    "strobulis", "fertile_leaflet_with_sorophores", "young_sterile_leaflet", "sporophore_with_sporangia",
    "aerial_stem", "mature_aerial_stem", "syangia",
    "fertile_frond", "mature_rhizome", "sterile_frond", "young_rhizome",
    "floating_leaves", "sporocarp", "submerged_leaves"  # ← new tissues
  ),
  stringsAsFactors = FALSE
)
manual_add <- rbind(manual_add, data.frame(
  `Organ type` = "Leaves",
  Clean_Tissue = "sterile_young_leaflet",
  stringsAsFactors = FALSE
))




# Ensure column names match
colnames(manual_add)[colnames(manual_add) == "Organ.type"] <- "Organ type"

# Combine with organ_mapping
organ_mapping <- rbind(organ_mapping[, c("Organ type", "Clean_Tissue")], manual_add)

# Create mapping dictionary
mapping_dict <- setNames(organ_mapping$`Organ type`, organ_mapping$Clean_Tissue)

# Function to analyze a single counts file
analyze_counts_file <- function(counts_file, fasta_file) {
  
  # Extract species ID from filename (first 3 letters)
  species_id <- substr(basename(counts_file), 1, 3)
  cat("Analyzing", species_id, "from file:", counts_file, "\n")
  
  # 1. Load the counts matrix
  counts <- read.csv(counts_file, header = TRUE, row.names = 1, check.names = FALSE)
  
  # 2. Extract tissue information from column names
  samples <- colnames(counts)
  tissue_names <- sub("_[0-9]+$", "", samples)
  standardized_tissues <- gsub("[^a-z0-9]+", "_", tolower(tissue_names))
  
  # Check for missing mappings before proceeding
  missing_tissues <- unique(standardized_tissues[is.na(mapping_dict[standardized_tissues])])
  if(length(missing_tissues) > 0) {
    cat("Missing mappings for tissues:\n")
    print(missing_tissues)
    stop("Some tissues could not be mapped to organ types. Please update the mapping.")
  }
  
  # Map tissue name to organ type
  organ_types <- mapping_dict[standardized_tissues]
  
  # Use organ types as the grouping variable
  tissue <- factor(organ_types)
  
  # 3. Create a DGEList object
  y <- DGEList(counts = counts, group = tissue)
  
  # 4. Filter out lowly expressed genes
  keep <- filterByExpr(y)
  y <- y[keep, , keep.lib.sizes = FALSE]
  
  # 5. Normalize library sizes
  y <- calcNormFactors(y)
  
  # 6. Create the design matrix
  design <- model.matrix(~0 + tissue)
  colnames(design) <- levels(tissue)
  
  # 7. Estimate dispersion
  y <- estimateDisp(y, design)
  
  # 8. Fit the GLM
  fit <- glmFit(y, design)
  
  # --- Extract and clean species IDs from FASTA ---
  species_genes_clean <- NULL
  if (!is.null(fasta_file) && file.exists(fasta_file)) {
    extract_fasta_ids <- function(fasta_file) {
      fasta_lines <- readLines(fasta_file)
      header_lines <- fasta_lines[grepl("^>", fasta_lines)]
      gene_ids <- gsub("^>", "", header_lines)
      gene_ids <- sapply(strsplit(gene_ids, " "), function(x) x[1])
      return(gene_ids)
    }
    
    fasta_genes <- extract_fasta_ids(fasta_file)
    species_genes <- fasta_genes[grepl(paste0("^", species_id), fasta_genes)]
    
    # Remove `.v2.1`, `.p1`, `.p`, or similar suffixes
    species_genes_clean <- gsub("\\.v2\\.1$", "", species_genes)
    species_genes_clean <- gsub("\\.p1$", "", species_genes_clean)
    species_genes_clean <- gsub("\\.p$", "", species_genes_clean)
    species_genes_clean <- unique(species_genes_clean)
    
    cat("Found", length(species_genes_clean), "SPL genes for", species_id, "in FASTA\n")
  } else {
    # If no FASTA file, use all genes starting with species_id from the counts
    species_genes_clean <- rownames(y$counts)[grepl(paste0("^", species_id), rownames(y$counts))]
    cat("FASTA file not found - using", length(species_genes_clean), "genes starting with", species_id, "from counts file\n")
  }
  
  # --- PAIRWISE COMPARISONS BETWEEN ALL ORGAN TYPES ---
  organ_types_unique <- levels(tissue)
  n_organs <- length(organ_types_unique)
  
  # Initialize storage for pairwise results
  pairwise_results <- list()
  pairwise_species_up <- list()
  pairwise_species_down <- list()
  
  # Create all pairwise contrasts
  for (i in 1:(n_organs-1)) {
    for (j in (i+1):n_organs) {
      organ1 <- organ_types_unique[i]
      organ2 <- organ_types_unique[j]
      comparison_name <- paste0(organ1, "_vs_", organ2)
      
      # Create contrast vector (organ1 vs organ2)
      contrast_vec <- rep(0, n_organs)
      contrast_vec[i] <- 1
      contrast_vec[j] <- -1
      
      # Perform likelihood ratio test
      lrt <- glmLRT(fit, contrast = contrast_vec)
      
      # Get DE genes
      de_genes <- topTags(lrt, n = Inf)$table
      up_genes <- rownames(de_genes)[de_genes$FDR < 0.05 & de_genes$logFC > 1]
      down_genes <- rownames(de_genes)[de_genes$FDR < 0.05 & de_genes$logFC < -1]
      
      # Store results
      pairwise_results[[comparison_name]] <- list(
        up = up_genes,
        down = down_genes,
        all_results = de_genes
      )
      
      # Find species-specific SPL genes
      species_up <- intersect(species_genes_clean, up_genes)
      species_down <- intersect(species_genes_clean, down_genes)
      
      pairwise_species_up[[comparison_name]] <- species_up
      pairwise_species_down[[comparison_name]] <- species_down
      
      cat(comparison_name, ": ", length(up_genes), " upregulated in ", organ1,
          " (", length(species_up), " ", species_id, " SPL); ",
          length(down_genes), " upregulated in ", organ2,
          " (", length(species_down), " ", species_id, " SPL)\n")
    }
  }
  
  # --- CREATE OVERVIEW TABLES FOR SPECIES SPL GENES ---
  create_species_overview <- function(species_lists, direction) {
    all_species <- unique(unlist(species_lists))
    
    if (length(all_species) == 0) {
      cat("No", species_id, "SPL genes found for", direction, "regulation\n")
      return(NULL)
    }
    
    overview_matrix <- matrix(0, nrow = length(all_species), ncol = length(species_lists))
    rownames(overview_matrix) <- all_species
    colnames(overview_matrix) <- names(species_lists)
    
    for (i in seq_along(species_lists)) {
      if (length(species_lists[[i]]) > 0) {
        overview_matrix[species_lists[[i]], i] <- 1
      }
    }
    
    overview_df <- as.data.frame(overview_matrix)
    overview_df$Total_Comparisons <- rowSums(overview_df)
    overview_df$Gene_ID <- rownames(overview_df)
    overview_df <- overview_df[, c("Gene_ID", names(species_lists), "Total_Comparisons")]
    
    return(overview_df)
  }
  
  # Create overview tables
  species_up_overview <- create_species_overview(pairwise_species_up, "up")
  species_down_overview <- create_species_overview(pairwise_species_down, "down")
  
  # Create output directory for this species
  output_dir <- paste0(species_id, "_organ_type_analysis_results")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # Save overview tables
  if (!is.null(species_up_overview)) {
    write.table(species_up_overview, file.path(output_dir, paste0(species_id, "_upregulated_overview.txt")), 
                sep = "\t", quote = FALSE, row.names = FALSE)
    cat(species_id, "upregulated SPL genes overview saved\n")
  }
  
  if (!is.null(species_down_overview)) {
    write.table(species_down_overview, file.path(output_dir, paste0(species_id, "_downregulated_overview.txt")), 
                sep = "\t", quote = FALSE, row.names = FALSE)
    cat(species_id, "downregulated SPL genes overview saved\n")
  }
  
  # --- CREATE SMEAR PLOTS ---
  create_smear_plot <- function(lrt_result, comparison_name, species_genes_clean, species_id) {
    results <- topTags(lrt_result, n = Inf)$table
    
    results$Gene_Type <- "Other"
    results$Gene_Type[rownames(results) %in% species_genes_clean] <- paste0(species_id, " SPL")
    
    results$Significant <- "No"
    results$Significant[results$FDR < 0.05 & abs(results$logFC) > 1] <- "Yes"
    
    p <- ggplot(results, aes(x = logCPM, y = logFC)) +
      geom_point(aes(color = Gene_Type, alpha = Significant), size = 0.8) +
      scale_color_manual(values = setNames(c("red", "black"), c(paste0(species_id, " SPL"), "Other"))) +
      scale_alpha_manual(values = c("Yes" = 0.8, "No" = 0.3)) +
      geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "blue") +
      geom_hline(yintercept = 0, linetype = "solid", color = "gray") +
      labs(title = paste("Smear Plot:", comparison_name),
           x = "Average log CPM",
           y = "log Fold Change",
           color = "Gene Type",
           alpha = "Significant") +
      theme_minimal() +
      theme(plot.title = element_text(size = 10))
    
    return(p)
  }
  
  # Create smear plots for all pairwise comparisons
  smear_plots <- list()
  for (comparison in names(pairwise_results)) {
    organs_in_comparison <- strsplit(comparison, "_vs_")[[1]]
    organ1_idx <- which(organ_types_unique == organs_in_comparison[1])
    organ2_idx <- which(organ_types_unique == organs_in_comparison[2])
    
    contrast_vec <- rep(0, n_organs)
    contrast_vec[organ1_idx] <- 1
    contrast_vec[organ2_idx] <- -1
    
    lrt <- glmLRT(fit, contrast = contrast_vec)
    smear_plots[[comparison]] <- create_smear_plot(lrt, comparison, species_genes_clean, species_id)
  }
  
  # Save individual smear plots
  for (i in seq_along(smear_plots)) {
    ggsave(file.path(output_dir, paste0("smear_plot_", names(smear_plots)[i], ".png")), 
           smear_plots[[i]], width = 8, height = 6, dpi = 300)
  }
  
  # Create combined smear plots
  if (length(smear_plots) <= 6) {
    combined_plot <- do.call(grid.arrange, c(smear_plots, ncol = 2))
    ggsave(file.path(output_dir, "all_smear_plots_combined.png"), combined_plot, 
           width = 16, height = 4 * ceiling(length(smear_plots)/2), dpi = 300)
  } else {
    plots_per_page <- 6
    n_pages <- ceiling(length(smear_plots) / plots_per_page)
    
    for (page in 1:n_pages) {
      start_idx <- (page - 1) * plots_per_page + 1
      end_idx <- min(page * plots_per_page, length(smear_plots))
      page_plots <- smear_plots[start_idx:end_idx]
      
      combined_plot <- do.call(grid.arrange, c(page_plots, ncol = 2))
      ggsave(file.path(output_dir, paste0("smear_plots_page_", page, ".png")), combined_plot, 
             width = 16, height = 4 * ceiling(length(page_plots)/2), dpi = 300)
    }
  }
  
  # --- SAVE PAIRWISE RESULTS ---
  for (comparison in names(pairwise_results)) {
    write.table(pairwise_results[[comparison]]$up, 
                file.path(output_dir, paste0("upregulated_", comparison, ".txt")),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(pairwise_results[[comparison]]$down, 
                file.path(output_dir, paste0("downregulated_", comparison, ".txt")),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    write.table(pairwise_species_up[[comparison]], 
                file.path(output_dir, paste0(species_id, "_SPL_upregulated_", comparison, ".txt")),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(pairwise_species_down[[comparison]], 
                file.path(output_dir, paste0(species_id, "_SPL_downregulated_", comparison, ".txt")),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    write.table(pairwise_results[[comparison]]$all_results, 
                file.path(output_dir, paste0("full_results_", comparison, ".txt")),
                sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  }
  
  # --- CREATE MATRICES ---
  species_up_matrix <- matrix(0, nrow = n_organs, ncol = n_organs,
                              dimnames = list(organ_types_unique, organ_types_unique))
  species_down_matrix <- matrix(0, nrow = n_organs, ncol = n_organs,
                                dimnames = list(organ_types_unique, organ_types_unique))
  
  for (comparison in names(pairwise_species_up)) {
    organs_in_comp <- strsplit(comparison, "_vs_")[[1]]
    o1 <- organs_in_comp[1]
    o2 <- organs_in_comp[2]
    
    up_count <- length(pairwise_species_up[[comparison]])
    down_count <- length(pairwise_species_down[[comparison]])
    
    species_up_matrix[o1, o2] <- up_count
    species_up_matrix[o2, o1] <- down_count
    
    species_down_matrix[o1, o2] <- down_count
    species_down_matrix[o2, o1] <- up_count
  }
  
  write.csv(species_up_matrix, file.path(output_dir, paste0(species_id, "_SPL_upregulated_by_organ_matrix.csv")))
  write.csv(species_down_matrix, file.path(output_dir, paste0(species_id, "_SPL_downregulated_by_organ_matrix.csv")))
  
  cat("\nMatrix of UPREGULATED", species_id, "SPL genes between organ types:\n")
  print(species_up_matrix)
  
  # --- CREATE SUMMARY PLOT ---
  species_up_down_summary <- data.frame(
    Comparison = names(pairwise_species_up),
    Upregulated = sapply(pairwise_species_up, length),
    Downregulated = sapply(pairwise_species_down, length)
  )
  
  species_melted <- melt(species_up_down_summary, id.vars = "Comparison", 
                         variable.name = "Direction", value.name = "SPL_Gene_Count")
  species_melted$Comparison <- factor(species_melted$Comparison,
                                      levels = species_up_down_summary$Comparison[order(rowSums(species_up_down_summary[, -1]))])
  
  species_plot <- ggplot(species_melted, aes(x = Comparison, y = SPL_Gene_Count, fill = Direction)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.5) +
    labs(
      title = paste0(species_id, " SPL genes up- and downregulated in pairwise organ comparisons"),
      x = "Organ Comparison",
      y = "Number of SPL Genes"
    ) +
    scale_fill_manual(values = c("Upregulated" = "#D73027", "Downregulated" = "#4575B4")) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 11),
      axis.ticks.length.x = unit(0.3, "cm"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.margin = margin(10, 20, 10, 20)
    )
  
  ggsave(file.path(output_dir, paste0(species_id, "_SPL_summary_by_organ_plot.png")), species_plot, 
         width = 12, height = 8, dpi = 300)
  
  cat("Analysis complete for", species_id, "- results saved in", output_dir, "\n\n")
  
  return(list(
    species_id = species_id,
    output_dir = output_dir,
    pairwise_results = pairwise_results,
    species_up_matrix = species_up_matrix,
    species_down_matrix = species_down_matrix
  ))
}

# --- MAIN EXECUTION ---

# Set the FASTA file path
fasta_file <- "F:/CoBi/sysbio/final_sequences_for_mafft.fa"

# Check if FASTA file exists
if (!file.exists(fasta_file)) {
  stop("FASTA file not found at: ", fasta_file)
}

# Find all counts files in current directory
counts_files <- list.files(pattern = "*_merged_kallisto_counts.csv$", full.names = TRUE)

cat("Found", length(counts_files), "counts files to analyze:\n")
print(counts_files)
cat("Using FASTA file:", fasta_file, "\n\n")

# Analyze each file
all_results <- list()
for (counts_file in counts_files) {
  result <- analyze_counts_file(counts_file, fasta_file)
  all_results[[result$species_id]] <- result
}

cat("All analyses complete!\n")
cat("Results saved in separate directories for each species:\n")
for (result in all_results) {
  cat("-", result$output_dir, "\n")
}