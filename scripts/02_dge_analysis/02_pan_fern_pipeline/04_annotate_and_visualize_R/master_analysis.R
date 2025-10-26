# ---- AGGREGATE RESULTS ACROSS SPECIES ----
cat("Aggregating results across all species...\n")

# Create master results directory
master_output_dir <- "master_analysis"
if (!dir.exists(master_output_dir)) dir.create(master_output_dir)

# Initialize master data structures
all_species_results <- list()
master_de_table <- data.frame()
master_logfc_table <- data.frame()
spl_de_summary <- data.frame()

# Collect results from each species
for (sub in subfolders) {
  species <- basename(sub)
  summary_file <- file.path(sub, "pairwise_dge_summary.txt")
  
  if (!file.exists(summary_file)) next
  
  cat("Collecting results from", species, "...\n")
  
  # Read the summary
  species_summary <- read.table(summary_file, header = TRUE, sep = "\t")
  
  # Read the individual results for SPL filtering
  counts_file <- file.path(sub, paste0(species, "_merged_kallisto_counts.csv"))
  if (!file.exists(counts_file)) next
  
  counts <- read.csv(counts_file, row.names = 1)
  species_fasta_ids <- intersect(all_fasta_ids, rownames(counts))
  
  # Process each contrast for this species
  for (contrast in species_summary$Contrast) {
    # Create tissue groups from contrast name
    tissues_in_contrast <- strsplit(contrast, "_vs_")[[1]]
    tissue <- factor(sub("_[0-9]+$", "", colnames(counts)))
    
    # Recreate the analysis to get logFC values for SPL genes
    y <- DGEList(counts = counts, group = tissue)
    keep <- filterByExpr(y)
    y <- y[keep, , keep.lib.sizes = FALSE]
    y <- calcNormFactors(y)
    
    design <- model.matrix(~0 + tissue)
    colnames(design) <- levels(tissue)
    y <- estimateDisp(y, design)
    fit <- glmFit(y, design)
    
    # Create contrast vector
    contrast_vec <- rep(0, length(levels(tissue)))
    contrast_vec[which(levels(tissue) == tissues_in_contrast[1])] <- 1
    contrast_vec[which(levels(tissue) == tissues_in_contrast[2])] <- -1
    
    lrt <- glmLRT(fit, contrast = contrast_vec)
    
    # Filter for SPL genes only
    spl_genes_in_contrast <- intersect(species_fasta_ids, rownames(lrt$table))
    
    if (length(spl_genes_in_contrast) > 0) {
      # Extract results for SPL genes
      spl_results <- lrt$table[spl_genes_in_contrast, ]
      spl_results$species <- species
      spl_results$contrast <- contrast
      spl_results$gene_id <- rownames(spl_results)
      spl_results$is_de <- spl_results$FDR < 0.05
      spl_results$regulation <- ifelse(spl_results$is_de, 
                                      ifelse(spl_results$logFC > 0, "up", "down"), 
                                      "not_de")
      
      # Add to master tables
      for (gene in spl_genes_in_contrast) {
        comparison_name <- paste(species, contrast, sep = "_")
        
        # DE status table
        if (nrow(master_de_table) == 0) {
          master_de_table <- data.frame(gene_id = gene, stringsAsFactors = FALSE)
        }
        if (!gene %in% master_de_table$gene_id) {
          master_de_table <- rbind(master_de_table, data.frame(gene_id = gene, stringsAsFactors = FALSE))
        }
        
        gene_row <- which(master_de_table$gene_id == gene)
        master_de_table[gene_row, comparison_name] <- spl_results[gene, "regulation"]
        
        # LogFC table
        if (nrow(master_logfc_table) == 0) {
          master_logfc_table <- data.frame(gene_id = gene, stringsAsFactors = FALSE)
        }
        if (!gene %in% master_logfc_table$gene_id) {
          master_logfc_table <- rbind(master_logfc_table, data.frame(gene_id = gene, stringsAsFactors = FALSE))
        }
        
        gene_row <- which(master_logfc_table$gene_id == gene)
        master_logfc_table[gene_row, comparison_name] <- spl_results[gene, "logFC"]
      }
    }
  }
}

# Save master tables
write.table(master_de_table, file = file.path(master_output_dir, "master_spl_de_table.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(master_logfc_table, file = file.path(master_output_dir, "master_spl_logfc_table.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)
