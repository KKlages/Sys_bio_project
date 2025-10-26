# Differential Gene Expression Analysis of SPL Genes in Ferns

## Project Overview

This repository contains the complete bioinformatics pipeline for analyzing differential gene expression (DGE) of SQUAMOSA-PROMOTER BINDING PROTEIN-LIKE (SPL) transcription factors across 22+ fern species. The analysis integrates tissue-specific and environmental RNA-seq data to characterize SPL gene expression patterns and their evolutionary significance.


**Principal Investigator:** Dr. Sophie de Vries  
**Student:** Kilian Klages  
**Date:** July 3, 2025

---

## Table of Contents

1. [Project Structure](#project-structure)
2. [Prerequisites](#prerequisites)
3. [Data Sources](#data-sources)
4. [Pipeline Overview](#pipeline-overview)
5. [Detailed Workflows](#detailed-workflows)
6. [Expected Outputs](#expected-outputs)
7. [Computational Requirements](#computational-requirements)
8. [Troubleshooting](#troubleshooting)
9. [Citation](#citation)

---

## Project Structure
```
scripts/02_dge_analysis/
├── 01_single_dataset_pipeline/      # Analysis of individual experiments
│   ├── 01_download_data/            # SRA data retrieval
│   ├── 02_run_quantification/       # kallisto quantification
│   ├── 03_postprocess_and_merge/    # Count matrix generation
│   └── 04_run_dge/                  # Statistical analysis (edgeR)
│
├── 02_pan_fern_pipeline/            # Cross-species comparative analysis
│   ├── 01_prepare_and_quantify/     # Bulk processing of 22 species
│   ├── 02_collate_and_merge_results/# Data integration
│   ├── 03_run_dge_analysis_R/       # DGE by tissue and organ type
│   └── 04_annotate_and_visualize_R/ # Final figures and tables
│
└── README.md                        # This file
```

---

## Prerequisites

### Software Dependencies

#### Command-line Tools
```bash
# Sequence data tools
sra-toolkit >= 3.0.0        # prefetch, fastq-dump
fastp >= 0.23.0             # Quality control and trimming
seqkit >= 2.3.0             # FASTA/Q manipulation

# Quantification
kallisto >= 0.48.0          # Pseudoalignment and quantification

# Utilities
GNU parallel >= 20220522    # Parallel processing
```

#### R Environment (version >= 4.2.0)
```r
# Bioconductor packages
BiocManager::install(c(
  "edgeR",           # >= 3.40.0 - Differential expression
  "limma"            # >= 3.54.0 - Linear models
))

# CRAN packages
install.packages(c(
  "ggplot2",         # >= 3.4.0 - Plotting
  "pheatmap",        # >= 1.0.12 - Heatmaps
  "gridExtra",       # >= 2.3 - Plot arrangement
  "reshape2",        # >= 1.4.4 - Data reshaping
  "dplyr",           # >= 1.1.0 - Data manipulation
  "tidyr",           # >= 1.3.0 - Data tidying
  "readr",           # >= 2.1.0 - Data import
  "readxl",          # >= 1.4.0 - Excel files
  "VennDiagram",     # >= 1.7.3 - Venn diagrams
  "RColorBrewer"     # >= 1.1.3 - Color palettes
))
```

#### Python Environment (version >= 3.8)
```bash
pip install biopython pandas
```

### Installation
```bash
# Clone repository
git clone [repository_url]
cd scripts/02_dge_analysis

# Create conda environment (recommended)
conda env create -f environment.yml
conda activate fern_rnaseq_env

# Or install manually
conda install -c bioconda fastp kallisto sra-tools seqkit
```

---

## Data Sources

### RNA-seq Datasets

| Experiment | Species | Accession | Conditions | Reference |
|------------|---------|-----------|------------|-----------|
| Azolla FR light | *A. filiculoides* | PRJEB45223 | Far-red light exposure | Dijkhuizen et al. 2021 |
| Azolla symbiosis | *A. filiculoides* | PRJNA430527 | Cyanobacteria ± Nitrogen | Li et al. 2018 |
| Ceratopteris tissues | *C. richardii* | [SRA IDs] | 10 tissue types | Marchant et al. 2022 |
| Pan-fern | 22 fern species | Various | Tissue-specific | Ali et al. 2025 |

### Reference Data

- **Transcriptomes**: fernbase.org (C. richardii, A. filiculoides)
- **SPL annotations**: From phylogenetic analysis (see `../01_phylogenetics/`)
- **Organ type mapping**: `organ_type.xlsx` (manually curated)
- **Clade annotations**: `Genes_in_clades.xlsx` (from phylogeny)

---

## Pipeline Overview

### Workflow Diagram
```
Raw SRA Data
    ↓
[1. Download & Convert] → FASTQ files
    ↓
[2. Quality Control]    → fastp (trimming, filtering)
    ↓
[3. Quantification]     → kallisto (TPM & counts)
    ↓
[4. Count Matrix]       → Merge replicates
    ↓
[5. DGE Analysis]       → edgeR (pairwise comparisons)
    ↓
[6. Annotation]         → Add clade/organ info
    ↓
[7. Visualization]      → Heatmaps, MA plots, tables
```

---

## Detailed Workflows

### Workflow 1: Single Dataset Analysis

**Purpose:** Analyze individual experiments (e.g., Azolla far-red light response)

#### Step 1: Download Data
```bash
cd 01_single_dataset_pipeline/01_download_data/

# Option A: Using SRR IDs (for small datasets)
bash get_srr_ids.sh

# Option B: Using download links (for large datasets)
# 1. Edit download_links.csv with your URLs
python split.py          # Split into batches
bash loop.sh             # Download and process
```

**Expected output:** `SRR*/` directories containing `.sra` files

#### Step 2: Quantification
```bash
cd ../02_run_quantification/

# Create kallisto index (do once per species)
kallisto index -i azolla.idx path/to/transcriptome.fasta

# Run pipeline (handles: SRA→FASTQ→QC→quantification)
bash kallisto_pipeline.sh
```

**Configuration required:**
- Edit `kallisto_pipeline.sh`:
  - Set `KALLISTO_INDEX_FILENAME` to your index
  - Set `SRR_SAMPLES` array to your accessions
  - Adjust `THREADS` for your system

**Expected output:** `results/SRR*/kallisto_quant/abundance.tsv`

#### Step 3: Post-processing
```bash
cd ../03_postprocess_and_merge/

# Rename files by experimental condition
python rename_by_dict.py  # Edit dict for your experiment

# Handle biological replicates
python handle_replicates.py

# Create count matrix for edgeR
python merge_counts_to_matrix.py
```

**Expected output:** `counts.txt` (genes × samples)

#### Step 4: Differential Expression
```bash
cd ../04_run_dge/

# Choose analysis based on your experiment:

# For Azolla cyanobacteria/nitrogen experiment
Rscript dge_azfi_cn.R

# For Azolla far-red light experiment
Rscript dge_azfi_fr.R

# For Ceratopteris tissue comparison
Rscript dge_ceri_tissue.R
```

**Key parameters in R scripts:**
```r
# In each script, update:
fasta_file <- "path/to/final_sequences_for_mafft.fa"  # SPL gene IDs
FDR_threshold <- 0.05                                 # Significance cutoff
logFC_threshold <- 1                                  # 2-fold change
```

**Expected outputs:**
- MA plots with SPL genes highlighted
- DE gene lists (up/down regulated)
- Summary tables and matrices
- Venn diagrams (where applicable)

---

### Workflow 2: Pan-Fern Comparative Analysis

**Purpose:** Compare SPL expression across 22 fern species and multiple tissues

#### Step 1: Prepare Data
```bash
cd 02_pan_fern_pipeline/01_prepare_and_quantify/

# Organize downloaded FASTQ files
bash organize_raw_fastq.sh

# Build kallisto indices for all species
bash build_kallisto_indices.sh
```

**File structure required:**
```
Species_abbreviation/
├── transcriptome.fasta
└── ERR_or_SRR_ID/
    ├── *_1.fastq.gz
    └── *_2.fastq.gz
```

#### Step 2: Quantify All Species
```bash
# Submit parallel SLURM jobs (HPC)
bash slurm_submission/submit_all_jobs.sh

# Or run locally (slower)
for species in Adi Aev Ala ...; do
  cd $species
  bash run_quant.sh
  cd ..
done
```

**Expected runtime:** 2-4 hours per species (with 8 threads)

#### Step 3: Collate Results
```bash
cd ../02_collate_and_merge_results/

# Rename abundance files by tissue type
bash rename_from_metadata.sh

# Collect all TSV files into analysis directory
bash collect_renamed_files.sh

# Generate master metadata file
bash generate_final_metadata.sh
```

**Expected output:** 
- `analysis/Species/tissue_replicate.tsv` files
- `combined_metadata.tsv`

#### Step 4: Run DGE Analysis
```bash
cd ../03_run_dge_analysis_R/

# Option A: Tissue-level analysis (fine-grained)
Rscript dge_by_tissue.R

# Option B: Organ-level analysis (grouped tissues)
Rscript dge_by_organ_type.R

# Option C: Comprehensive analysis with plotting
Rscript fern_dge_and_plots.R
```

**Key configurations:**
```r
# In R scripts, verify:
setwd("path/to/analysis/")
fasta_file <- "final_sequences_for_mafft.fa"

# For organ-type analysis:
organ_mapping <- read_excel("organ_type.xlsx")  # Must be present
```

**Expected outputs per species:**
- `Species_analysis_results/` directory
  - MA plots for all pairwise comparisons
  - Gene lists (up/down regulated)
  - Count matrices
  - Summary tables

#### Step 5: Annotate and Visualize
```bash
cd ../04_annotate_and_visualize_R/

# Create tables of annotated genes
Rscript create_table.R

# Summarize regulation patterns
Rscript summarize_regulation.R

# Add clade annotations
Rscript annotate_dge_results.R

# Generate final heatmaps (Figure 6)
Rscript create_final_heatmap.R
```

**Required input files:**
- `Genes_in_clades.xlsx` - SPL gene to superclade mapping
- `organ_type.xlsx` - Tissue to organ type mapping
- All DGE results from Step 4

**Expected outputs:**
- `heatmap_by_clade.png` - Main visualization
- `heatmap_by_order.png`
- `heatmap_by_family.png`
- Summary tables with statistics

---

## Expected Outputs

### Primary Results Files

1. **Count Matrices**
   - `*_merged_kallisto_counts.csv` - Raw counts per species
   - `counts.txt` - edgeR input format

2. **DGE Results**
   - `annotated_results_*.txt` - Full results per comparison
   - `*_upregulated_overview.txt` - Summary of upregulated genes
   - `*_downregulated_overview.txt` - Summary of downregulated genes

3. **Figures**
   - `smear_plot_*.png` - MA plots showing SPL genes
   - `heatmap_*.png` - Expression patterns across species
   - `*_summary_plot.png` - Bar plots of DE gene counts

4. **Tables**
   - `gene_regulation_summary.txt` - Per-gene DE status
   - `pairwise_dge_summary.txt` - Comparison statistics
   - `*_matrix.csv` - Count matrices for visualization

### File Naming Conventions

- Species codes: 3-letter abbreviations (Adi, Aev, Ala, etc.)
- Tissues: snake_case (sterile_leaflet, young_root)
- Comparisons: `tissue1_vs_tissue2`
- Directions: `upregulated`, `downregulated`

---

## Computational Requirements

### Minimum Requirements

- **CPU:** 8 cores
- **RAM:** 16 GB
- **Storage:** 100 GB (for single species analysis)
- **Storage:** 1 TB (for full pan-fern analysis)

### Recommended (HPC)

- **CPU:** 16+ cores per job
- **RAM:** 32 GB per job
- **Scheduler:** SLURM (scripts provided)
- **Partitions:** Medium (24-48 hours runtime)

### Runtime Estimates

| Step | Single Species | All Species (22) |
|------|---------------|------------------|
| Download | 1-2 hours | 10-20 hours |
| Quantification | 2-4 hours | 48-96 hours (parallel) |
| DGE Analysis | 30 min | 10-12 hours |
| Visualization | 10 min | 1-2 hours |

**Tip:** Use SLURM array jobs for parallel species processing:
```bash
sbatch --array=1-22 process_species.slurm
```

---

## Troubleshooting

### Common Issues

#### 1. kallisto Index Not Found

**Error:** `ERROR: Kallisto index not found`

**Solution:**
```bash
# Create index for your transcriptome
kallisto index -i species.idx transcriptome.fasta

# Verify index
ls -lh species.idx  # Should be 10-100 MB
```

#### 2. Missing Gene IDs in DGE Results

**Error:** `No SPL genes found for species XXX`

**Cause:** Gene ID format mismatch between FASTA and count matrix

**Solution:**
```r
# Check formats
head(rownames(counts))        # e.g., "Adi_g104042"
head(species_fasta_ids)       # Should match exactly

# If mismatched, clean IDs:
species_fasta_ids <- gsub("\\.v2\\.1$", "", species_fasta_ids)
species_fasta_ids <- gsub("\\.p$", "", species_fasta_ids)
```

#### 3. Tissue to Organ Mapping Errors

**Error:** `Missing mappings for tissues: [list]`

**Solution:**
```r
# Add missing tissues to manual_add in dge_by_organ_type.R
manual_add <- data.frame(
  `Organ type` = c("Roots", "Leaves"),
  Clean_Tissue = c("your_missing_tissue", "another_tissue")
)
```

#### 4. Memory Errors in R

**Error:** `cannot allocate vector of size X Gb`

**Solutions:**
```r
# Increase memory limit (Unix)
options(java.parameters = "-Xmx16g")

# Process species individually instead of all at once
# Use data.table for large datasets
library(data.table)
counts <- fread("large_counts.csv")
```

#### 5. SLURM Job Failures

**Error:** `Job exceeded memory/time limit`

**Solution:**
```bash
# Increase resources in .slurm file
#SBATCH --mem=64G          # Was 32G
#SBATCH --time=72:00:00    # Was 48:00:00

# Check actual usage
sacct -j JOBID --format=JobID,MaxRSS,Elapsed
```

---

## Quality Control Checklist

Before finalizing results, verify:

- [ ] All expected output files generated
- [ ] No FASTQ files with anomalously low read counts
- [ ] kallisto quantification success rate > 95%
- [ ] Number of mapped reads per sample documented
- [ ] PCA/MDS plots show expected clustering by condition
- [ ] Biological replicates correlate (Pearson R > 0.95)
- [ ] Number of DE genes is reasonable (not 0, not 50%+ of genome)
- [ ] MA plots show expected symmetry around logFC = 0
- [ ] SPL genes are present in final results
- [ ] Heatmaps reproduce patterns in report figures

---

## Key Analysis Decisions

### Filtering Thresholds
```r
# Gene filtering (edgeR::filterByExpr defaults)
min_count <- 10              # Minimum count across all samples
min_total_count <- 15        # Minimum total count per gene
min_prop <- 0.7              # Minimum proportion of samples

# Significance thresholds
FDR_cutoff <- 0.05           # False discovery rate
logFC_cutoff <- 1            # 2-fold change (log2)
```

### Statistical Methods

- **Normalization:** TMM (edgeR default)
- **Dispersion:** GLM with empirical Bayes shrinkage
- **Testing:** Likelihood ratio test (LRT) for multi-factor, QLF for simple designs
- **Multiple testing:** Benjamini-Hochberg FDR correction

---

## File Formats

### Count Matrix Format
```
target_id    sample1    sample2    sample3
gene_001     150        200        180
gene_002     50         45         52
```

### Metadata Format
```
species    sample       tissue          replicate
Adi        leaves_1     leaves          1
Adi        leaves_2     leaves          2
Adi        root_1       root            1
```

### Annotation Format (Excel)
```
Superclade SPL1/12    Superclade SPL7    Superclade SPL8
Adi_g12345           Adi_g67890         Adi_g11111
Aev_g22222           Aev_g33333         Aev_g44444
```

---

## Citation

If you use this pipeline, please cite:
```
Klages, K. (2025). A phylogenomic and differential gene expression 
analysis of the SPL gene family in ferns. Internship report for 
M.Bio.310 Bioinformatik der Systembiologie. University of Göttingen.
```

**Associated Publications:**
- Ali et al. (2025). Nature Plants 11, 1028–1048
- Li et al. (2018). Nature Plants 4, 460–472
- Dijkhuizen et al. (2021). Front. Plant Sci. 12, 693039

---

## Contact

**Analyst:** Kilian Klages  
**Email:** kilianklages@icloud.com  
**Supervisor:** Dr. Sophie de Vries  
**Institution:** University of Göttingen

---

## License

[Specify license - e.g., MIT, GPL-3.0, CC-BY-4.0]

---

## Acknowledgments

- Dr. Sophie de Vries for supervision and streptophyte/alga sequence data
- de Vries group HPC cluster for computational resources
- fernbase.org for reference transcriptomes
- Authors of original RNA-seq datasets (see Data Sources)

---

## Changelog

### Version 1.0 (2025-07-03)
- Initial release with complete pipeline
- Phylogenetic and DGE analysis workflows
- Pan-fern comparative analysis
- Automated visualization scripts

---

**Last Updated:** October 26, 2025
