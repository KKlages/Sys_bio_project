library(edgeR)
library(limma)

# ---- CONFIGURATION ----
setwd("/mnt/ceph-hdd/projects/scc_ubmg_devries/KKlages/mrna_subproject/analysis")
fasta_file <- "final_sequences_for_mafft.fa"

# ---- GENE ID EXTRACTION ----
extract_fasta_ids <- function(fasta_file) {
  lines <- readLines(fasta_file)
  header_lines <- lines[grepl("^>", lines)]
  ids <- gsub("^>", "", header_lines)
  ids <- sub("\\.v2\\.1$", "", ids)
  ids <- sub("\\.p1$", "", ids)
  ids <- sub("\\.p$", "", ids)
  unique(ids)
}

all_fasta_ids <- extract_fasta_ids(fasta_file)

# ---- PROCESS ALL SPECIES ----
subfolders <- list.dirs(path = ".", recursive = FALSE, full.names = TRUE)

for (sub in subfolders) {
  species <- basename(sub)
  counts_file <- file.path(sub, paste0(species, "_merged_kallisto_counts.csv"))
  if (!file.exists(counts_file)) next
  
  cat("Processing", species, "...\n")
  
  # Read counts
  counts <- read.csv(counts_file, row.names = 1)
  sample_names <- colnames(counts)
  tissue <- factor(sub("_[0-9]+$", "", sample_names))
  
  # Filter FASTA IDs
  species_fasta_ids <- intersect(all_fasta_ids, rownames(counts))
  
  # DGE analysis
  y <- DGEList(counts = counts, group = tissue)
  keep <- filterByExpr(y)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  design <- model.matrix(~0 + tissue)
  colnames(design) <- levels(tissue)
  y <- estimateDisp(y, design)
  fit <- glmFit(y, design)
  
  # All pairwise contrasts
  tissues <- levels(tissue)
  n_tissues <- length(tissues)
  
  # Store results for annotated genes only
  all_annotated_results <- list()
  
  for(i in 1:(n_tissues-1)) {
    for(j in (i+1):n_tissues) {
      comparison <- paste(tissues[i], "vs", tissues[j], sep = "_")
      
      # Create contrast vector
      contrast_vec <- rep(0, n_tissues)
      contrast_vec[i] <- 1
      contrast_vec[j] <- -1
      
      # Perform test
      lrt <- glmLRT(fit, contrast = contrast_vec)
      
      # Get all results
      all_results <- topTags(lrt, n = Inf)$table
      
      # Filter to annotated genes only
      annotated_results <- all_results[rownames(all_results) %in% species_fasta_ids, ]
      
      # Add gene IDs as a column
      annotated_results$gene_id <- rownames(annotated_results)
      
      # Reorder columns
      annotated_results <- annotated_results[, c("gene_id", "logFC", "logCPM", "LR", "PValue", "FDR")]
      
      # Save individual comparison results
      write.table(annotated_results, 
                  file = file.path(sub, paste0("annotated_results_", comparison, ".txt")),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      
      # Store for summary
      all_annotated_results[[comparison]] <- annotated_results
      
      cat("  ", comparison, ": ", nrow(annotated_results), " annotated genes\n")
    }
  }
  
  # Create summary table
  summary_data <- data.frame(
    Species = species,
    Comparison = names(all_annotated_results),
    Total_annotated_genes = sapply(all_annotated_results, nrow),
    DE_annotated_genes = sapply(all_annotated_results, function(x) sum(x$FDR < 0.05)),
    Upregulated = sapply(all_annotated_results, function(x) sum(x$FDR < 0.05 & x$logFC > 0)),
    Downregulated = sapply(all_annotated_results, function(x) sum(x$FDR < 0.05 & x$logFC < 0)),
    stringsAsFactors = FALSE
  )
  
  # Save summary
  write.table(summary_data, 
              file = file.path(sub, "annotated_genes_summary.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat("Completed", species, "\n\n")
}

cat("All annotated gene results saved!\n")
