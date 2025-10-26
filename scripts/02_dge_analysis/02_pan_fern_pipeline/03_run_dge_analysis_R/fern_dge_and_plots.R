library(edgeR)
library(limma)
library(gridExtra)
library(ggplotify)
library(ggplot2)

# ---- CONFIGURATION ----
setwd("")
combined_meta <- read.table("combined_metadata.tsv", header = TRUE, sep = "\t")
fasta_file <- "final_sequences_for_mafft.fa" # <-- Set your fasta file path here

# ---- GENE ID EXTRACTION ----
extract_fasta_ids <- function(fasta_file) {
  lines <- readLines(fasta_file)
  header_lines <- lines[grepl("^>", lines)]
  ids <- gsub("^>", "", header_lines)
  ids <- sub("\\.v2\\.1$", "", ids)
  ids <- sub("\\.p1$", "", ids)  # Remove .p1 suffix
  ids <- sub("\\.p$", "", ids)
  unique(ids)
}

all_fasta_ids <- extract_fasta_ids(fasta_file)
cat("Total unique gene IDs in FASTA file:", length(all_fasta_ids), "\n")

# ---- PROCESS ALL SPECIES ----
subfolders <- list.dirs(path = ".", recursive = FALSE, full.names = TRUE)

for (sub in subfolders) {
  species <- basename(sub)
  # Look for merged kallisto counts file instead of summed counts
  counts_file <- file.path(sub, paste0(species, "_merged_kallisto_counts.csv"))
  if (!file.exists(counts_file)) next
  
  cat("Processing", species, "...\n")
  
  # Prepare output directory for plots
  output_dir <- file.path(sub, "pemp")
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  # Read counts - individual replicates
  counts <- read.csv(counts_file, row.names = 1)
  
  # Extract tissue groups from column names (remove replicate numbers)
  sample_names <- colnames(counts)
  tissue <- factor(sub("_[0-9]+$", "", sample_names))
  
  cat("  - Sample names:", paste(sample_names, collapse = ", "), "\n")
  cat("  - Tissue groups:", paste(levels(tissue), collapse = ", "), "\n")
  
  # Filter FASTA IDs to only those present in this species' counts
  species_fasta_ids <- intersect(all_fasta_ids, rownames(counts))
  cat("  - FASTA genes for", species, ":", length(species_fasta_ids), "\n")
  cat("  - Total genes in counts:", nrow(counts), "\n")
  cat("  - Overlap:", length(species_fasta_ids), "/", nrow(counts), 
      "(", round(length(species_fasta_ids)/nrow(counts)*100, 1), "%)\n")
  
  # DGEList and filtering with individual replicates
  y <- DGEList(counts = counts, group = tissue)
  keep <- filterByExpr(y)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  
  # Use individual samples in design matrix
  design <- model.matrix(~0 + tissue)
  colnames(design) <- levels(tissue)
  y <- estimateDisp(y, design)
  fit <- glmFit(y, design)
  
  # All pairwise contrasts
  tissues <- levels(tissue)
  n_tissues <- length(tissues)
  contrasts_list <- list()
  contrast_names <- c()
  for(i in 1:(n_tissues-1)) {
    for(j in (i+1):n_tissues) {
      contrast_name <- paste(tissues[i], "vs", tissues[j], sep = "_")
      contrast_names <- c(contrast_names, contrast_name)
      contrast_vec <- rep(0, n_tissues)
      contrast_vec[i] <- 1
      contrast_vec[j] <- -1
      contrasts_list[[contrast_name]] <- contrast_vec
    }
  }
  
  # DGE for each contrast and plot
  results_list <- list()
  de_genes_list <- list()
  plot_list <- list()
  tissue_combinations <- combn(tissues, 2, simplify = FALSE)
  for (k in seq_along(tissue_combinations)) {
    pair <- tissue_combinations[[k]]
    comparison <- paste(pair, collapse = "_vs_")
    contrast_vec <- contrasts_list[[comparison]]
    lrt <- glmLRT(fit, contrast = contrast_vec)
    genes_in_data <- rownames(lrt$table)
    
    # Store results for summary
    results_list[[comparison]] <- lrt$table
    
    # Determine DE genes (FDR < 0.05)
    de <- decideTests(lrt, p.value = 0.05)
    de_genes <- rownames(lrt$table)[as.logical(de)]
    de_genes_list[[comparison]] <- de_genes
    
    # Find annotated genes present in this comparison
    annotated_genes_in_data <- intersect(species_fasta_ids, genes_in_data)
    annotated_de_genes <- intersect(species_fasta_ids, de_genes)
    
    # Count annotated genes for this comparison
    cat("    ", comparison, ": annotated genes =", length(annotated_genes_in_data), 
        ", DE annotated =", length(annotated_de_genes), "\n")
    
    # Create enhanced smear plot with higher resolution
    png(file.path(output_dir, paste0("smear_", comparison, ".png")), 
        width = 2400, height = 1800, res = 300)  # Increased resolution
    
    # Create the base smear plot
    plotSmear(lrt, de.tags = species_fasta_ids,
              main = paste("MA Plot:", comparison, "- Annotated genes (from FASTA) as DE tags"),
              cex = 0.5, pch = 16)

    
    # Add horizontal lines at ±1 log fold-change (2-fold change)
    abline(h = c(-1, 1), col = "blue", lty = 2, lwd = 2)
    
    # Highlight annotated genes if any exist
    if(length(annotated_genes_in_data) > 0) {
      # Get indices and values for annotated genes
      gene_indices <- match(annotated_genes_in_data, rownames(lrt$table))
      logFC_vals <- lrt$table$logFC[gene_indices]
      logCPM_vals <- lrt$AveLogCPM[gene_indices]
      
      # Add annotated genes as larger red points
      points(logCPM_vals, logFC_vals, col = "red", pch = 16, cex = 1.5)
      
      # Highlight annotated DE genes with blue outline
      if(length(annotated_de_genes) > 0) {
        de_gene_indices <- match(annotated_de_genes, rownames(lrt$table))
        de_logFC_vals <- lrt$table$logFC[de_gene_indices]
        de_logCPM_vals <- lrt$AveLogCPM[de_gene_indices]
        
        # Add blue outline to DE annotated genes
        points(de_logCPM_vals, de_logFC_vals, col = "blue", pch = 1, cex = 2, lwd = 3)
      }
      
      # Add labels for DE annotated genes (limit to avoid clutter)
      if(length(annotated_de_genes) > 0 & length(annotated_de_genes) <= 400) {  
        de_gene_indices <- match(annotated_de_genes, rownames(lrt$table))
        de_logFC_vals <- lrt$table$logFC[de_gene_indices]
        de_logCPM_vals <- lrt$AveLogCPM[de_gene_indices]

        text(de_logCPM_vals, de_logFC_vals, labels = annotated_de_genes, 
             pos = 3, cex = 0.6, col = "darkred", font = 2)
      }
    }
    
    # Create comprehensive legend
    legend_items <- c("Annotated genes", "Non-annotated genes", "Annotated DE genes")
    legend_colors <- c("red", "black", "blue")
    legend_pch <- c(16, 1, 1)
    legend_cex_vals <- c(1.2, 0.8, 1.5)

    
    legend("topright", 
           legend = legend_items,
           col = legend_colors,
           pch = legend_pch,
           pt.cex = legend_cex_vals,
           cex = 0.9,
           bg = "white",
           box.lwd = 1)
    
    # Add text summary
    mtext(paste("Annotated genes:", length(annotated_genes_in_data), 
                "| Annotated DE:", length(annotated_de_genes)), 
          side = 3, line = 0.5, cex = 0.8, col = "darkblue")
    
    dev.off()
    
    # Create ggplot version for grid arrangement
    plot_data <- data.frame(
      logCPM = lrt$AveLogCPM,
      logFC = lrt$table$logFC,
      gene_id = rownames(lrt$table),
      is_de = genes_in_data %in% de_genes,
      is_annotated = genes_in_data %in% species_fasta_ids
    )

    plot_data$color_group <- ifelse(plot_data$is_annotated, "annotated", "non_annotated")
    plot_data$color_group <- factor(plot_data$color_group, levels = c("annotated", "non_annotated"))

    p <- ggplot(plot_data, aes(x = logCPM, y = logFC, color = color_group)) +
      geom_point(size = 0.5, alpha = 0.7) +
      scale_color_manual(
        values = c(
          "annotated" = "blue",
          "non_annotated" = "black"
        ),
        labels = c("Annotated genes", "Non-annotated genes")
      ) +
      geom_hline(yintercept = c(-1, 1), color = "blue", linetype = "dashed", alpha = 0.5) +
      labs(title = comparison, x = "Average log CPM", y = "log FC") +
      theme_minimal() +
      theme(plot.title = element_text(size = 10))

    plot_list[[comparison]] <- p
  }
  
  # Summary
  summary_df <- data.frame(
    Contrast = contrast_names,
    DE_genes_count = sapply(contrast_names, function(x) length(de_genes_list[[x]]))
  )
  write.table(summary_df, file = file.path(sub, "pairwise_dge_summary.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Matrix of DE gene counts
  n_de_matrix <- matrix(0, nrow = n_tissues, ncol = n_tissues)
  rownames(n_de_matrix) <- tissues
  colnames(n_de_matrix) <- tissues
  for(i in 1:length(contrast_names)) {
    tissue_pair <- strsplit(contrast_names[i], "_vs_")[[1]]
    tissue1 <- tissue_pair[1]
    tissue2 <- tissue_pair[2]
    n_de <- length(de_genes_list[[contrast_names[i]]])
    n_de_matrix[tissue1, tissue2] <- n_de
    n_de_matrix[tissue2, tissue1] <- n_de
  }
  write.table(n_de_matrix, file = file.path(sub, "pairwise_de_matrix.txt"), sep = "\t", quote = FALSE)
  
  # Arrange and save plots (only if there are plots to arrange)
  if (length(plot_list) > 0) {
    n_cols <- min(3, length(plot_list))
    n_rows <- ceiling(length(plot_list) / n_cols)

    # Fill up with null grobs if needed
    total_cells <- n_cols * n_rows
    if (length(plot_list) < total_cells) {
      for (i in seq_len(total_cells - length(plot_list))) {
        plot_list[[length(plot_list) + 1]] <- grid::nullGrob()
      }
    }

    # Arrange everything
    combined_plot <- gridExtra::arrangeGrob(
      grobs = plot_list,
      nrow = n_rows,
      ncol = n_cols
    )

    # Save the grid
    ggsave(
      filename = file.path(output_dir, "all_smear_plots_grid.png"),
      plot = combined_plot,
      width = n_cols * 4,
      height = n_rows * 3,
      dpi = 200,
      limitsize = FALSE
    )
  }
  
  cat("Completed processing", species, "\n\n")
}

cat("All species processed successfully!\n")
