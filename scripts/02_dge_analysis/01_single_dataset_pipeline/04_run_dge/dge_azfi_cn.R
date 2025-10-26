library(edgeR)

# Load count matrix
x <- read.table("counts.txt",
                row.names = 1, header = TRUE, sep = "\t")

# Sample metadata
samples <- colnames(x)
# Extract experimental factors from sample names
cyanobacteria <- factor(gsub("^(c[01])n[01]_?[0-9]*$", "\\1", samples))  # c0 or c1
nitrogen      <- factor(gsub("^c[01](n[01])_?[0-9]*$", "\\1", samples))  # n0 or n1

# Rename factor levels for clarity
levels(cyanobacteria) <- c("noCyano", "cyano")
levels(nitrogen)      <- c("noNitrogen", "nitrogen")

# Create design matrix with interaction
design <- model.matrix(~ cyanobacteria * nitrogen)
colnames(design) <- c("Intercept", "Cyano", "Nitrogen", "Cyano_x_Nitrogen")

# Setup DGEList
y <- DGEList(counts = x)
y <- calcNormFactors(y)

# Filter lowly expressed genes
keep <- filterByExpr(y, design)
y <- y[keep, , keep.lib.sizes = FALSE]

# Estimate dispersion
y <- estimateDisp(y, design)

# Fit GLM
fit <- glmFit(y, design)

# Define contrasts for each comparison vs control (c0n0)
contrast1 <- c(0, 1, 0, 0)  # c1n0 vs c0n0
contrast2 <- c(0, 0, 1, 0)  # c0n1 vs c0n0
contrast3 <- c(0, 1, 1, 1)  # c1n1 vs c0n0

lrt1 <- glmLRT(fit, contrast = contrast1)
lrt2 <- glmLRT(fit, contrast = contrast2)
lrt3 <- glmLRT(fit, contrast = contrast3)

# Function to extract gene IDs from FASTA file
extract_fasta_ids <- function(fasta_file) {
  fasta_lines <- readLines(fasta_file)
  header_lines <- fasta_lines[grepl("^>", fasta_lines)]
  gene_ids <- gsub("^>", "", header_lines)
  gene_ids <- sapply(strsplit(gene_ids, " "), function(x) x[1])
  return(gene_ids)
}

# Load FASTA genes
#fasta_file <- "final_sequences_for_mafft.fa"
fasta_genes <- extract_fasta_ids(fasta_file)
azfi_genes <- fasta_genes[grepl("^Azfi", fasta_genes)]

# Find common genes with DGE results
dge_gene_names <- rownames(y)
common_genes <- intersect(azfi_genes, dge_gene_names)

print(paste("Found", length(common_genes), "common genes starting with 'Azfi'"))
print("Common genes:")
print(common_genes)

# Set plotting area to 3 rows, 1 column
# Set 3-row layout
par(mfrow = c(3, 1), mar = c(4, 4, 4, 2))  # margin bottom-left-top-right

# Helper function for a formatted MA plot
plot_custom_MA <- function(lrt_obj, title_text, azfi_de, azfi_non_de) {
  logCPM <- lrt_obj$table$logCPM
  logFC <- lrt_obj$table$logFC
  gene_names <- rownames(lrt_obj$table)
  
  # Set up empty plot
  plot(logCPM, logFC,
       type = "n",
       xlab = "Average log CPM",
       ylab = "log-fold-change",
       main = title_text,
       cex.main = 1.2)
  
  # Add threshold lines
  abline(h = c(-2, 2), col = "blue", lty = 2)
  abline(h = 0, col = "blue", lty = 1)
  
  # Plot all genes in black first
  points(logCPM, logFC, col = "black", pch = 16, cex = 0.8)
  
  # Add non-DE Azfi genes in blue
  if(length(azfi_non_de) > 0) {
    idx_non <- match(azfi_non_de, gene_names)
    idx_non <- idx_non[!is.na(idx_non)]
    points(logCPM[idx_non], logFC[idx_non], col = "blue", pch = 16, cex = 1.2)
  }
  
  # Add DE Azfi genes in red
  if(length(azfi_de) > 0) {
    idx_de <- match(azfi_de, gene_names)
    idx_de <- idx_de[!is.na(idx_de)]
    points(logCPM[idx_de], logFC[idx_de], col = "red", pch = 16, cex = 1.2)
  }
  
  # Add legend
  legend("topright", 
         legend = c("SPL DE genes", "SPL non-DE genes", "Other genes"),
         col = c("red", "blue", "black"),
         pch = 16,
         cex = 0.8)
}

# Define DE genes (FDR < 0.05)
de_genes1 <- rownames(topTags(lrt1, n = Inf, p.value = 0.05)$table)
de_genes2 <- rownames(topTags(lrt2, n = Inf, p.value = 0.05)$table)
de_genes3 <- rownames(topTags(lrt3, n = Inf, p.value = 0.05)$table)

# SPL DE and non-DE genes per contrast
azfi_de_1     <- intersect(common_genes, de_genes1)
azfi_non_de_1 <- setdiff(common_genes, de_genes1)

azfi_de_2     <- intersect(common_genes, de_genes2)
azfi_non_de_2 <- setdiff(common_genes, de_genes2)

azfi_de_3     <- intersect(common_genes, de_genes3)
azfi_non_de_3 <- setdiff(common_genes, de_genes3)

# Plot 1: Cyanobacteria effect
plot_custom_MA(
  lrt1,
  expression("MA Plot: Effect of Cyanobacteria (" * italic("c1n0 vs c0n0") * ")"),
  azfi_de_1,
  azfi_non_de_1
)

# Plot 2: Nitrogen effect
plot_custom_MA(
  lrt2,
  expression("MA Plot: Effect of Nitrogen (" * italic("c0n1 vs c0n0") * ")"),
  azfi_de_2,
  azfi_non_de_2
)

# Plot 3: Combined effect
plot_custom_MA(
  lrt3,
  expression("MA Plot: Combined Effect (" * italic("c1n1 vs c0n0") * ")"),
  azfi_de_3,
  azfi_non_de_3
)

# Reset layout
par(mfrow = c(1, 1))

# Function to classify up/down-regulated Azfi DE genes
classify_azfi_regulation <- function(lrt_obj, azfi_de_genes) {
  de_table <- lrt_obj$table[azfi_de_genes, ]
  up   <- rownames(de_table[de_table$logFC > 0, ])
  down <- rownames(de_table[de_table$logFC < 0, ])
  return(list(up = up, down = down))
}

# Classify for all 3 contrasts
azfi_reg_1 <- classify_azfi_regulation(lrt1, azfi_de_1)
azfi_reg_2 <- classify_azfi_regulation(lrt2, azfi_de_2)
azfi_reg_3 <- classify_azfi_regulation(lrt3, azfi_de_3)

# Create summary tables
up_table <- data.frame(
  Contrast = c("c1n0 vs c0n0", "c0n1 vs c0n0", "c1n1 vs c0n0"),
  Upregulated = c(length(azfi_reg_1$up), length(azfi_reg_2$up), length(azfi_reg_3$up))
)

down_table <- data.frame(
  Contrast = c("c1n0 vs c0n0", "c0n1 vs c0n0", "c1n1 vs c0n0"),
  Downregulated = c(length(azfi_reg_1$down), length(azfi_reg_2$down), length(azfi_reg_3$down))
)

# Print the tables
cat("\nUpregulated Azfi SPL genes:\n")
print(up_table)

cat("\nDownregulated Azfi SPL genes:\n")
print(down_table)

# Load if needed
if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  install.packages("VennDiagram")
}
library(VennDiagram)

# Define sets of DE Azfi SPL genes per contrast
azfi_de_sets <- list(
  Cyanobacteria = azfi_de_1,
  Nitrogen = azfi_de_2,
  Combined = azfi_de_3
)

# Print counts per set
cat("\nCounts of DE Azfi SPL genes per condition:\n")
sapply(azfi_de_sets, length)

# Compute overlaps
overlap_cyano_nitro     <- intersect(azfi_de_1, azfi_de_2)
overlap_cyano_combined  <- intersect(azfi_de_1, azfi_de_3)
overlap_nitro_combined  <- intersect(azfi_de_2, azfi_de_3)
overlap_all_three       <- Reduce(intersect, azfi_de_sets)

# Summary of overlaps
cat("\nOverlap Summary:\n")
cat("Cyanobacteria ∩ Nitrogen:", length(overlap_cyano_nitro), "\n")
cat("Cyanobacteria ∩ Combined:", length(overlap_cyano_combined), "\n")
cat("Nitrogen ∩ Combined:", length(overlap_nitro_combined), "\n")
cat("All three conditions:", length(overlap_all_three), "\n")

# Optional: Draw Venn Diagram
venn.plot <- venn.diagram(
  x = azfi_de_sets,
  filename = NULL,
  fill = c("red", "blue", "green"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  cat.col = c("red", "blue", "green"),
  margin = 0.1
)
grid.newpage()
grid.draw(venn.plot)
