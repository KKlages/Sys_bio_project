library(edgeR)
library(ggplot2)
library(gridExtra)

# 1. Load the counts matrix
counts <- read.table("counts.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# 2. Extract tissue information from column names
samples <- colnames(counts)
tissue <- factor(gsub("[0-9]+$", "", samples))

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

# --- Extract and clean Ceri IDs from FASTA ---
extract_fasta_ids <- function(fasta_file) {
  fasta_lines <- readLines(fasta_file)
  header_lines <- fasta_lines[grepl("^>", fasta_lines)]
  gene_ids <- gsub("^>", "", header_lines)
  gene_ids <- sapply(strsplit(gene_ids, " "), function(x) x[1])
  return(gene_ids)
}
#fasta_genes <- extract_fasta_ids("final_sequences_for_mafft.fa")
ceri_genes <- fasta_genes[grepl("^Ceri", fasta_genes)]
# Clean: remove .v2.1 and .p at end
ceri_genes_clean <- sub("\\.v2\\.1$", "", ceri_genes)
ceri_genes_clean <- sub("\\.p$", "", ceri_genes_clean)
ceri_genes_clean <- unique(ceri_genes_clean)

# --- PAIRWISE COMPARISONS BETWEEN ALL TISSUES ---
tissues <- levels(tissue)
n_tissues <- length(tissues)

# Initialize storage for pairwise results
pairwise_results <- list()
pairwise_ceri_up <- list()
pairwise_ceri_down <- list()

# Create all pairwise contrasts
for (i in 1:(n_tissues-1)) {
  for (j in (i+1):n_tissues) {
    tissue1 <- tissues[i]
    tissue2 <- tissues[j]
    comparison_name <- paste0(tissue1, "_vs_", tissue2)
    
    # Create contrast vector (tissue1 vs tissue2)
    contrast_vec <- rep(0, n_tissues)
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
    
    # Find Ceri genes
    ceri_up <- intersect(ceri_genes_clean, up_genes)
    ceri_down <- intersect(ceri_genes_clean, down_genes)
    
    pairwise_ceri_up[[comparison_name]] <- ceri_up
    pairwise_ceri_down[[comparison_name]] <- ceri_down
    
    cat(comparison_name, ": ", length(up_genes), " upregulated in ", tissue1,
        " (", length(ceri_up), " Ceri); ",
        length(down_genes), " upregulated in ", tissue2,
        " (", length(ceri_down), " Ceri)\n")
  }
}

# --- CREATE OVERVIEW TABLES FOR CERI GENES ---

# Function to create overview table
create_ceri_overview <- function(ceri_lists, direction) {
  # Get all unique Ceri genes
  all_ceri <- unique(unlist(ceri_lists))
  
  if (length(all_ceri) == 0) {
    cat("No Ceri genes found for", direction, "regulation\n")
    return(NULL)
  }
  
  # Create matrix
  overview_matrix <- matrix(0, nrow = length(all_ceri), ncol = length(ceri_lists))
  rownames(overview_matrix) <- all_ceri
  colnames(overview_matrix) <- names(ceri_lists)
  
  # Fill matrix
  for (i in seq_along(ceri_lists)) {
    if (length(ceri_lists[[i]]) > 0) {
      overview_matrix[ceri_lists[[i]], i] <- 1
    }
  }
  
  # Convert to data frame and add summary columns
  overview_df <- as.data.frame(overview_matrix)
  overview_df$Total_Comparisons <- rowSums(overview_df)
  overview_df$Gene_ID <- rownames(overview_df)
  
  # Reorder columns
  overview_df <- overview_df[, c("Gene_ID", names(ceri_lists), "Total_Comparisons")]
  
  return(overview_df)
}

# Create overview tables
ceri_up_overview <- create_ceri_overview(pairwise_ceri_up, "up")
ceri_down_overview <- create_ceri_overview(pairwise_ceri_down, "down")

# Save overview tables
if (!is.null(ceri_up_overview)) {
  write.table(ceri_up_overview, "ceri_upregulated_overview.txt", 
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Ceri upregulated genes overview:\n")
  print(head(ceri_up_overview))
}

if (!is.null(ceri_down_overview)) {
  write.table(ceri_down_overview, "ceri_downregulated_overview.txt", 
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Ceri downregulated genes overview:\n")
  print(head(ceri_down_overview))
}

# --- CREATE SMEAR PLOTS ---

# Function to create smear plot
create_smear_plot <- function(lrt_result, comparison_name, ceri_genes_clean) {
  # Get the results table
  results <- topTags(lrt_result, n = Inf)$table
  
  # Add gene categories
  results$Gene_Type <- "Other"
  results$Gene_Type[rownames(results) %in% ceri_genes_clean] <- "Ceri"
  
  # Add significance
  results$Significant <- "No"
  results$Significant[results$FDR < 0.05 & abs(results$logFC) > 1] <- "Yes"
  
  # Create the plot
  p <- ggplot(results, aes(x = logCPM, y = logFC)) +
    geom_point(aes(color = Gene_Type, alpha = Significant), size = 0.8) +
    scale_color_manual(values = c("Ceri" = "red", "Other" = "black")) +
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
  # Recreate the LRT for this comparison
  tissues_in_comparison <- strsplit(comparison, "_vs_")[[1]]
  tissue1_idx <- which(tissues == tissues_in_comparison[1])
  tissue2_idx <- which(tissues == tissues_in_comparison[2])
  
  contrast_vec <- rep(0, n_tissues)
  contrast_vec[tissue1_idx] <- 1
  contrast_vec[tissue2_idx] <- -1
  
  lrt <- glmLRT(fit, contrast = contrast_vec)
  smear_plots[[comparison]] <- create_smear_plot(lrt, comparison, ceri_genes_clean)
}

# Save individual smear plots
for (i in seq_along(smear_plots)) {
  ggsave(paste0("smear_plot_", names(smear_plots)[i], ".png"), 
         smear_plots[[i]], width = 8, height = 6, dpi = 300)
}

# Create a combined plot with all smear plots
if (length(smear_plots) <= 6) {
  combined_plot <- do.call(grid.arrange, c(smear_plots, ncol = 2))
  ggsave("all_smear_plots_combined.png", combined_plot, 
         width = 16, height = 4 * ceiling(length(smear_plots)/2), dpi = 300)
} else {
  # If too many plots, create multiple combined plots
  plots_per_page <- 6
  n_pages <- ceiling(length(smear_plots) / plots_per_page)
  
  for (page in 1:n_pages) {
    start_idx <- (page - 1) * plots_per_page + 1
    end_idx <- min(page * plots_per_page, length(smear_plots))
    page_plots <- smear_plots[start_idx:end_idx]
    
    combined_plot <- do.call(grid.arrange, c(page_plots, ncol = 2))
    ggsave(paste0("smear_plots_page_", page, ".png"), combined_plot, 
           width = 16, height = 4 * ceiling(length(page_plots)/2), dpi = 300)
  }
}

# --- SAVE PAIRWISE RESULTS ---
for (comparison in names(pairwise_results)) {
  # Save all DE genes
  write.table(pairwise_results[[comparison]]$up, 
              paste0("upregulated_", comparison, ".txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(pairwise_results[[comparison]]$down, 
              paste0("downregulated_", comparison, ".txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Save Ceri genes
  write.table(pairwise_ceri_up[[comparison]], 
              paste0("ceri_upregulated_", comparison, ".txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(pairwise_ceri_down[[comparison]], 
              paste0("ceri_downregulated_", comparison, ".txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Save full results table
  write.table(pairwise_results[[comparison]]$all_results, 
              paste0("full_results_", comparison, ".txt"),
              sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
}

# --- CREATE MATRIX OF CERI GENES BETWEEN TISSUE PAIRS (UP and DOWN) ---

# Initialize empty matrices
ceri_up_matrix <- matrix(0, nrow = n_tissues, ncol = n_tissues,
                         dimnames = list(tissues, tissues))
ceri_down_matrix <- matrix(0, nrow = n_tissues, ncol = n_tissues,
                           dimnames = list(tissues, tissues))

# Fill in the matrices using the pairwise_ceri_up and _down lists
for (comparison in names(pairwise_ceri_up)) {
  tissues_in_comp <- strsplit(comparison, "_vs_")[[1]]
  t1 <- tissues_in_comp[1]
  t2 <- tissues_in_comp[2]
  
  up_count <- length(pairwise_ceri_up[[comparison]])
  down_count <- length(pairwise_ceri_down[[comparison]])
  
  # Fill both directions for up and down
  ceri_up_matrix[t1, t2] <- up_count  # up in t1
  ceri_up_matrix[t2, t1] <- down_count  # up in t2
  
  ceri_down_matrix[t1, t2] <- down_count  # down in t2
  ceri_down_matrix[t2, t1] <- up_count  # down in t1
}

# Round to integers (optional)
ceri_up_matrix <- round(ceri_up_matrix)
ceri_down_matrix <- round(ceri_down_matrix)

# Save as CSV
write.csv(ceri_up_matrix, "ceri_upregulated_gene_matrix.csv")
write.csv(ceri_down_matrix, "ceri_downregulated_gene_matrix.csv")

# Also print one of them for inspection
cat("\nMatrix of UPREGULATED SPL genes between tissues:\n")
print(ceri_up_matrix)



# (Your original data prep code remains the same)
ceri_up_down_summary <- data.frame(
  Comparison = names(pairwise_ceri_up),
  Upregulated = sapply(pairwise_ceri_up, length),
  Downregulated = sapply(pairwise_ceri_down, length)
)
library(reshape2)
ceri_melted <- melt(ceri_up_down_summary, id.vars = "Comparison", 
                    variable.name = "Direction", value.name = "Ceri_Gene_Count")
ceri_melted$Comparison <- factor(ceri_melted$Comparison,
                                 levels = ceri_up_down_summary$Comparison[order(rowSums(ceri_up_down_summary[, -1]))])

# --- Store the plot in a variable ---
library(ggplot2)
ceri_plot <- ggplot(ceri_melted, aes(x = Comparison, y = Ceri_Gene_Count, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  # Remove coord_flip() to keep vertical bars
  labs(
    title = expression(italic("C. ri")~SPL~genes~up-~and~downregulated * 
                         " in pairwise tissue comparisons"),
    x = "Tissue Comparison",
    y = "Number of SPL Genes"
  ) +
  scale_fill_manual(values = c("Upregulated" = "#D73027", "Downregulated" = "#4575B4")) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),  # rotate x labels 45 degrees
    axis.text.y = element_text(size = 11),
    axis.ticks.length.x = unit(0.3, "cm"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(10, 20, 10, 20)
  )

