library(readxl)
library(tidyr)
library(dplyr)

# Load the Excel file (adjust path accordingly)
superclade_df <- read_excel("F:/CoBi/sysbio/data/Genes_in_clades.xlsx")

# Convert from wide to long: one gene per row with superclade info
superclade_long <- superclade_df %>%
  pivot_longer(
    cols = everything(),      # all columns
    names_to = "Superclade",  # new column for superclade names (former columns)
    values_to = "Gene_ID"     # new column for gene IDs (values)
  ) %>%
  filter(!is.na(Gene_ID), Gene_ID != "") %>%  # remove empty rows
  mutate(Gene_ID = trimws(Gene_ID))            # trim whitespace

# View the result
print(superclade_long)

library(dplyr)
library(tidyr)
library(stringr)
library(readr)

species_list <- c("Adi", "Aev", "Ala", "Als", "Aob", "Aop", "Cba", "Dac", "Dcu",
                  "Dde", "Ehy", "Lfl", "Len", "Msp", "Nbi", "Ore", "Pir", "Pnu",
                  "Ppi", "Sam", "Spa", "Tin")

overview_data_list <- list()

for (sp in species_list) {
  
  # Define file paths
  up_file <- file.path(paste0(sp, "_organ_type_analysis_results"), paste0(sp, "_upregulated_overview.txt"))
  down_file <- file.path(paste0(sp, "_organ_type_analysis_results"), paste0(sp, "_downregulated_overview.txt"))
  
  # Function to read and reshape overview files
  read_and_pivot <- function(file, species, direction) {
    if (!file.exists(file)) {
      warning(paste("File not found:", file))
      return(NULL)
    }
    df <- read_tsv(file, col_types = cols())
    
    if (!"Gene_ID" %in% colnames(df)) {
      warning(paste("No Gene_ID column in", file))
      return(NULL)
    }
    
    # Pivot longer: one row per gene per organ comparison
    df_long <- df %>%
      pivot_longer(
        cols = -c(Gene_ID, Total_Comparisons),
        names_to = "Comparison",
        values_to = "Present"
      ) %>%
      filter(Present == 1) %>%
      mutate(
        Species = species,
        Direction = direction
      )
    
    return(df_long)
  }
  
  # Read upregulated overview
  up_long <- read_and_pivot(up_file, sp, "Upregulated")
  # Read downregulated overview
  down_long <- read_and_pivot(down_file, sp, "Downregulated")
  
  # Combine if both exist
  combined <- bind_rows(up_long, down_long)
  
  overview_data_list[[sp]] <- combined
}

# Combine all species overview data
all_overview <- bind_rows(overview_data_list)

# Preview results
head(all_overview)


all_overview <- all_overview %>%
  mutate(Gene_ID = gsub("_", " ", Gene_ID))


# Now filter superclade genes that appear in all_overview
present_superclades <- superclade_long %>%
  filter(Gene_ID %in% all_overview$Gene_ID) %>%
  distinct(Superclade)

cat("Sample from all_overview:\n")
print(head(all_overview$Gene_ID, 10))

cat("\nSample from superclade_long:\n")
print(head(superclade_long$Gene_ID, 10))
# 2. Join overview with superclade info by Gene_ID
overview_with_sc <- all_overview %>%
  inner_join(superclade_long, by = "Gene_ID") %>%
  filter(Present == 1)  # keep only genes present in the comparison

# 3. Extract Organ from Comparison column, e.g. "Leaves_vs_Roots" → "Leaves"
overview_with_sc <- overview_with_sc %>%
  mutate(Organ = str_extract(Comparison, "^[^_]+"))

# 4. Summarize: count of genes per Superclade × Organ × Direction
summary_table <- overview_with_sc %>%
  group_by(Superclade, Organ, Direction) %>%
  summarise(Gene_Count = n(), .groups = "drop")

# 5. View summary
print(summary_table)

library(dplyr)
library(stringr)

# Clean gene IDs in overview (if not done)
all_overview <- all_overview %>%
  mutate(Gene_ID = gsub("_", " ", Gene_ID))

# Join overview with superclade info by Gene_ID
overview_with_sc <- all_overview %>%
  inner_join(superclade_long, by = "Gene_ID") %>%
  filter(Present == 1)  # only genes present in comparison

# Extract Organ from Comparison (first tissue before underscore)
overview_with_sc <- overview_with_sc %>%
  mutate(Organ = str_extract(Comparison, "^[^_]+"))

# Select relevant columns: Gene_ID, Organ, Direction, Superclade
gene_organ_direction <- overview_with_sc %>%
  select(Gene_ID, Organ, Direction, Superclade)

# Optionally arrange or remove duplicates
gene_organ_direction <- gene_organ_direction %>%
  distinct() %>%
  arrange(Superclade, Organ, Direction, Gene_ID)

# View result
print(gene_organ_direction)

# Save to CSV if you want
write.csv(gene_organ_direction, "Genes_per_Organ_and_Direction.csv", row.names = FALSE)

# Load required libraries
library(dplyr)
library(knitr)

# --- Step 1: Ensure the annotated symmetrical dataframe is loaded ---
# This dataframe has 650 rows, representing 650 regulation perspectives.
if (!exists("annotated_symmetrical_df")) {
  cat("Loading annotated symmetrical details from 'all_species_spl_annotated_symmetrical.csv'...\n")
  annotated_symmetrical_df <- readr::read_csv("all_species_spl_annotated_symmetrical.csv", col_types = readr::cols())
}

# --- Step 2: Create the Final Summary by Counting Rows ---

# This is the simplest and most direct method.
# We are NOT using distinct(). We are counting every row in each group.
# This guarantees the sum of the final table will equal the number of rows
# in the input dataframe (650).

final_summary_650 <- annotated_symmetrical_df %>%
  # Group by the categories for the final table. Note we use `Organ = Organ`
  # which was already renamed from Organ1 in your script.
  group_by(Superclade, Organ, Direction) %>%
  # Use n() to COUNT THE ROWS in each group.
  summarise(
    Count = n(),
    .groups = "drop" # Ungroup after summarising
  ) %>%
  # Arrange for readability
  arrange(Superclade, Organ, desc(Direction))


# --- Step 3: Print the Final Table (Sum = 650) ---

cat("\n--- Final Summary Table (Total Count = 650) ---\n\n")
cat(paste("The sum of counts in this table is:", sum(final_summary_650$Count), "\n\n"))

kable(
  final_summary_650,
  caption = "Count of Gene Regulation Events by Superclade, Organ, and Direction.",
  col.names = c("Superclade", "Organ", "Direction", "Count")
)

# --- Verification Step ---
# Let's verify that the total Up and Down counts now balance perfectly across the whole dataset.
verification_summary <- final_summary_650 %>%
  group_by(Direction) %>%
  summarise(Total_Count = sum(Count))

cat("\n\n--- Verification of Balanced Counts (Total Up vs. Total Down) ---\n")
print(verification_summary)

# Load required libraries
library(dplyr)
library(tidyr)
library(readr)
library(stringr)

# --- 1. SETUP ---

# Define the list of all species IDs to process
species_list <- c("Adi", "Aev", "Ala", "Als", "Aob", "Aop", "Cba", "Dac", "Dcu",
                  "Dde", "Ehy", "Lfl", "Len", "Msp", "Nbi", "Ore", "Pir", "Pnu",
                  "Ppi", "Sam", "Spa", "Tin")

# This list will hold the detailed dataframe from each species
all_species_details_list <- list()

# --- 2. FUNCTION TO EXTRACT GENE IDs FOR EACH COMPARISON ---

# This function is identical to the one from the previous script.
# It processes one species at a time.
extract_gene_details <- function(matrix_data, direction, results_path, species) {
  
  all_comparisons <- list()
  organ_names <- rownames(matrix_data)
  
  for (organ1 in organ_names) {
    for (organ2 in organ_names) {
      gene_count <- matrix_data[organ1, organ2]
      
      if (gene_count == 0 || organ1 == organ2) {
        next
      }
      
      comparison_name <- paste0(organ1, "_vs_", organ2)
      
      direction_for_file <- if (direction == "Upregulated") "upregulated" else "downregulated"
      
      gene_list_file <- file.path(
        results_path,
        paste0(species, "_SPL_", direction_for_file, "_", comparison_name, ".txt")
      )
      
      if (!file.exists(gene_list_file)) {
        warning(paste("File not found for", species, ", skipping:", basename(gene_list_file)))
        next
      }
      
      gene_ids <- readLines(gene_list_file, warn = FALSE)
      
      if (length(gene_ids) > 0) {
        comparison_df <- tibble(
          Species = species,
          Organ1 = organ1,
          Organ2 = organ2,
          Direction = direction,
          Gene_ID = gene_ids
        )
        all_comparisons[[length(all_comparisons) + 1]] <- comparison_df
      }
    }
  }
  
  return(bind_rows(all_comparisons))
}

# --- 3. LOOP THROUGH ALL SPECIES AND PROCESS THEIR DATA ---

cat("Starting data extraction for all species...\n")

for (sp in species_list) {
  
  cat(paste0("--- Processing ", sp, " ---\n"))
  
  # Define paths for the current species
  results_dir <- paste0(sp, "_organ_type_analysis_results")
  up_matrix_file <- file.path(results_dir, paste0(sp, "_SPL_upregulated_by_organ_matrix.csv"))
  down_matrix_file <- file.path(results_dir, paste0(sp, "_SPL_downregulated_by_organ_matrix.csv"))
  
  # Check if required files exist before proceeding
  if (!dir.exists(results_dir)) {
    warning(paste("Results directory not found for", sp, ". Skipping."))
    next
  }
  if (!file.exists(up_matrix_file) || !file.exists(down_matrix_file)) {
    warning(paste("Matrix files not found for", sp, ". Skipping."))
    next
  }
  
  # Read the matrices
  up_matrix <- read.csv(up_matrix_file, row.names = 1, check.names = FALSE)
  down_matrix <- read.csv(down_matrix_file, row.names = 1, check.names = FALSE)
  
  # Extract detailed gene information for the current species
  upregulated_details_df <- extract_gene_details(up_matrix, "Upregulated", results_dir, sp)
  downregulated_details_df <- extract_gene_details(down_matrix, "Downregulated", results_dir, sp)
  
  # Combine the up and down results for the current species
  species_gene_details <- bind_rows(upregulated_details_df, downregulated_details_df)
  
  # Add the resulting dataframe to our master list
  if (nrow(species_gene_details) > 0) {
    all_species_details_list[[sp]] <- species_gene_details
  }
  
  cat(paste0("Finished processing ", sp, ". Found ", nrow(species_gene_details), " gene regulation events.\n"))
}

# --- 4. COMBINE ALL SPECIES DATA INTO A SINGLE MASTER DATAFRAME ---

cat("\nCombining results from all species into a single master dataframe...\n")
master_gene_details_df <- bind_rows(all_species_details_list)

# --- 5. VIEW THE FINAL DATAFRAME ---

cat("\n--- Final Master Detailed Gene Dataframe ---\n")
cat(paste("Total rows (gene regulation events) found across all species:", nrow(master_gene_details_df), "\n"))
cat(paste("Total unique genes found:", n_distinct(master_gene_details_df$Gene_ID), "\n\n"))

# Print the first few rows of the final result
print(head(master_gene_details_df))

# Print the last few rows to see data from other species
print(tail(master_gene_details_df))

# Save the final dataframe to a CSV file for future use (e.g., joining with annotations)
write.csv(master_gene_details_df, "all_species_spl_regulation_details.csv", row.names = FALSE)
cat("\nMaster dataframe saved to 'all_species_spl_regulation_details.csv'\n")

# Load required libraries
library(dplyr)
library(readr) # In case you are loading the file from disk

# --- Step 1: Ensure the master dataframe is loaded ---

# If you have the 'master_gene_details_df' in your environment from the last script, you can skip this.
# If you are starting a new session, load it from the CSV file you saved.
if (!exists("master_gene_details_df")) {
  cat("Loading master details from 'all_species_spl_regulation_details.csv'...\n")
  master_gene_details_df <- read_csv("all_species_spl_regulation_details.csv", col_types = cols())
}

# --- Step 2: Create the mirrored dataframe ---

cat("Creating the mirrored version of the dataframe...\n")

# Create a copy of the original dataframe to modify
mirrored_df <- master_gene_details_df

# Swap the Organ1 and Organ2 columns
# We use a temporary variable 'temp_organ' to hold the original Organ1 values
mirrored_df <- mirrored_df %>%
  mutate(
    temp_organ = Organ1,
    Organ1 = Organ2,
    Organ2 = temp_organ
  ) %>%
  select(-temp_organ) # Remove the temporary column

# Flip the 'Direction'
mirrored_df <- mirrored_df %>%
  mutate(
    Direction = if_else(Direction == "Upregulated", "Downregulated", "Upregulated")
  )

# --- Step 3: Combine the original and mirrored dataframes ---

cat("Combining original and mirrored dataframes into a final symmetrical table...\n")

# Use bind_rows() to stack the original and the new mirrored data on top of each other
symmetrical_gene_details_df <- bind_rows(master_gene_details_df, mirrored_df) %>%
  # Arrange the final dataframe for better readability
  arrange(Species, Gene_ID, Organ1)

# --- Step 4: Verify and View the Results ---

cat("\n--- Verification of the Symmetrical Dataframe ---\n")
cat("Original dataframe rows:", nrow(master_gene_details_df), "\n")
cat("Symmetrical dataframe rows:", nrow(symmetrical_gene_details_df), "(Should be exactly double)\n\n")

# Let's look at our example gene 'Adi_g104042' to see the result
cat("--- Example: Adi_g104042 ---\n")
verification_example <- symmetrical_gene_details_df %>%
  filter(Gene_ID == "Adi_g104042")

print(verification_example)

# Display the first few rows of the new, complete dataframe
cat("\n--- Head of the Symmetrical Dataframe ---\n")
print(head(symmetrical_gene_details_df, 10))

# Save the symmetrical dataframe to a new file for future use
write.csv(symmetrical_gene_details_df, "all_species_spl_symmetrical_regulation.csv", row.names = FALSE)
cat("\nSymmetrical dataframe saved to 'all_species_spl_symmetrical_regulation.csv'\n")

# Load required libraries
library(dplyr)
library(readr)
library(readxl)
library(tidyr)
library(stringr)

# --- 1. Load the Symmetrical Gene Details Data ---

# This dataframe should be in your environment from the previous step.
# If not, load it from the saved file.
if (!exists("symmetrical_gene_details_df")) {
  cat("Loading symmetrical details from 'all_species_spl_symmetrical_regulation.csv'...\n")
  symmetrical_gene_details_df <- read_csv("all_species_spl_symmetrical_regulation.csv", col_types = cols())
}

# --- 2. Load and Prepare the Superclade Annotation Data ---

# Define the path to your Excel file
superclade_file_path <- "F:/CoBi/sysbio/data/Genes_in_clades.xlsx"

if (!file.exists(superclade_file_path)) {
  stop("Annotation file not found at:", superclade_file_path)
}

cat("Loading and processing superclade annotations...\n")
superclade_df <- read_excel(superclade_file_path)

# Convert from wide to long format AND FIX THE GENE ID FORMAT
superclade_long <- superclade_df %>%
  pivot_longer(
    cols = everything(),
    names_to = "Superclade",
    values_to = "Gene_ID"
  ) %>%
  filter(!is.na(Gene_ID) & Gene_ID != "") %>%
  mutate(Gene_ID = trimws(Gene_ID)) %>%
  # ------------------------------------------------------------------ #
  # <<< THIS IS THE CRITICAL FIX >>>
  # Replace the first space with an underscore to match the other dataframe
  # e.g., 'Ehy g05026' becomes 'Ehy_g05026'
  mutate(Gene_ID = str_replace(Gene_ID, " ", "_"))
# ------------------------------------------------------------------ #

# --- Verification Step: Check if the formats now match ---
cat("\n--- Verifying Gene ID Formats ---\n")
cat("Sample from symmetrical_gene_details_df:\n")
print(head(symmetrical_gene_details_df$Gene_ID, 3))
cat("\nSample from superclade_long (after fixing):\n")
print(head(superclade_long$Gene_ID, 3))
cat("\n")


# --- 3. Annotate the Symmetrical Dataframe ---

cat("Annotating the symmetrical dataframe with superclade information...\n")

# Use a LEFT JOIN to add Superclade info to every gene regulation event.
# This will now work correctly because the IDs match.
annotated_symmetrical_df <- symmetrical_gene_details_df %>%
  left_join(superclade_long, by = "Gene_ID")

# --- 4. Handle Unannotated Genes ---

# After the join, genes not found in 'superclade_long' will have NA for Superclade.
# Replace these NAs with the label "Unannotated".
annotated_symmetrical_df <- annotated_symmetrical_df %>%
  mutate(
    Superclade = if_else(is.na(Superclade), "Unannotated", Superclade)
  )

# --- 5. Final Touches and Verification ---

# Reorder columns for better readability
annotated_symmetrical_df <- annotated_symmetrical_df %>%
  select(Species, Gene_ID, Superclade, Organ = Organ1, Direction, Compared_To = Organ2) %>%
  arrange(Species, Gene_ID, Organ)


cat("\n--- Final Annotated Symmetrical Dataframe ---\n")
cat(paste("Total rows:", nrow(annotated_symmetrical_df), "\n"))
cat(paste("Total unique genes:", n_distinct(annotated_symmetrical_df$Gene_ID), "\n\n"))

# View the head of the final, fully annotated table
print(head(annotated_symmetrical_df, 10))

# Verify that both annotated and unannotated genes are present
cat("\n--- Superclade Summary (Post-Fix) ---\n")
print(count(annotated_symmetrical_df, Superclade))

# Save the final, most useful dataframe to a file
write.csv(annotated_symmetrical_df, "all_species_spl_annotated_symmetrical.csv", row.names = FALSE)
cat("\nFully annotated symmetrical dataframe saved to 'all_species_spl_annotated_symmetrical.csv'\n")
