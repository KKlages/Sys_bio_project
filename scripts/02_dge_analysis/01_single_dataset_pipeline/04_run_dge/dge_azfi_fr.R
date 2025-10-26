library(edgeR)


# Load count matrix
x <- read.table("counts.txt", row.names = 1, header = TRUE, sep = "\t")

# Define groups
group <- factor(c(rep("FRX", 3), rep("TLX", 3)))

# Create DGEList object
y <- DGEList(counts = x, group = group)

# Filter lowly expressed genes
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes = FALSE]

# Normalize library sizes
y <- calcNormFactors(y)

# Design matrix
design <- model.matrix(~group)

# Estimate dispersion (this is the missing step)
y <- estimateDisp(y, design)

# Fit the quasi-likelihood model
fit <- glmQLFit(y, design)

# Test for differential expression (TLX vs FRX)
qlf <- glmQLFTest(fit, coef = 2)

# View top DE genes
topTags(qlf)

# Save results
write.table(topTags(qlf, n = Inf), file = "dge_FRX_vs_TLX.txt", sep = "\t", quote = FALSE)


# Assuming edgeR and your analysis (e.g., y, qlf, detags) are already done
detags <- rownames(topTags(qlf, n = Inf, p.value = 0.05)$table)

# Extract gene IDs from FASTA
extract_fasta_ids <- function(fasta_file) {
  fasta_lines <- readLines(fasta_file)
  header_lines <- fasta_lines[grepl("^>", fasta_lines)]
  gene_ids <- gsub("^>", "", header_lines)
  gene_ids <- sapply(strsplit(gene_ids, " "), function(x) x[1])
  return(gene_ids)
}

# Extract gene IDs from your FASTA file (use the correct filename)
fasta_file <- "F:/CoBi/sysbio/final_sequences_for_mafft.fa"
fasta_genes <- extract_fasta_ids(fasta_file)

# Filter for genes starting with "Azfi"
azfi_genes <- fasta_genes[grepl("^Azfi", fasta_genes)]

# Identify DE and non-DE SPL genes
dge_gene_names <- rownames(y)
common_genes <- intersect(azfi_genes, dge_gene_names)
azfi_de_genes <- intersect(common_genes, detags)
azfi_non_de_genes <- setdiff(common_genes, detags)

# --- Plot MA Plot ---
main_title <- expression("MA Plot with "*italic("A. filiculoides")*" SPL genes highlighted")

# Create the base plot without any points
plot(qlf$table$logCPM, qlf$table$logFC, 
     type = "n",  # Don't plot points yet
     xlab = "Average log CPM", 
     ylab = "log-fold-change",
     main = main_title)

# Add horizontal reference lines
abline(h = c(-2, 2), col = "blue", lty = 2)
abline(h = 0, col = "blue", lty = 1)

# Plot all points in black first
points(qlf$table$logCPM, qlf$table$logFC, 
       col = "black", 
       pch = 16, 
       cex = 0.8)

# Add DE SPL genes (red) on top
if(length(azfi_de_genes) > 0) {
  idx_de <- match(azfi_de_genes, rownames(qlf$table))
  points(qlf$table$logCPM[idx_de], qlf$table$logFC[idx_de], col = "red", pch = 16, cex = 1)
}

# Add non-DE SPL genes (blue) on top
if(length(azfi_non_de_genes) > 0) {
  idx_non_de <- match(azfi_non_de_genes, rownames(qlf$table))
  points(qlf$table$logCPM[idx_non_de], qlf$table$logFC[idx_non_de], col = "blue", pch = 16, cex = 1)
}

# Add legend
legend("topright", 
       legend = c("SPL DE genes", "SPL non-DE genes", "Other genes"),
       col = c("red", "blue", "black"), 
       pch = c(16, 16, 16),
       cex = 0.8)
