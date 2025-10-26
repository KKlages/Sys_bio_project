# Load required libraries
library(pheatmap)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(tibble)
library(readr)

# Define species list
species_list <- c("Adi", "Aev", "Ala", "Als", "Aob", "Aop", "Cba", "Dac", "Dcu",
                  "Dde", "Ehy", "Lfl", "Len", "Msp", "Nbi", "Ore", "Pir", "Pnu",
                  "Ppi", "Sam", "Spa", "Tin")

# Initialize lists to store data
up_data_list <- list()
down_data_list <- list()

# Read data files with error handling
for (sp in species_list) {
  up_file <- file.path(paste0(sp, "_organ_type_analysis_results"), 
                       paste0(sp, "_SPL_upregulated_by_organ_matrix.csv"))
  down_file <- file.path(paste0(sp, "_organ_type_analysis_results"), 
                         paste0(sp, "_SPL_downregulated_by_organ_matrix.csv"))
  
  # Read upregulated data
  if (file.exists(up_file)) {
    tryCatch({
      up_df <- read.csv(up_file, row.names = 1, check.names = FALSE)
      if (nrow(up_df) > 0 && ncol(up_df) > 0) {
        up_long <- up_df %>%
          rownames_to_column("Organ1") %>%
          pivot_longer(cols = -Organ1, names_to = "Organ2", values_to = "Count") %>%
          mutate(Species = sp, Count = as.numeric(Count))
        up_data_list[[sp]] <- up_long
      }
    }, error = function(e) {
      warning(paste("Error reading", up_file, ":", e$message))
    })
  } else {
    warning(paste("File not found:", up_file))
  }
  
  # Read downregulated data
  if (file.exists(down_file)) {
    tryCatch({
      down_df <- read.csv(down_file, row.names = 1, check.names = FALSE)
      if (nrow(down_df) > 0 && ncol(down_df) > 0) {
        down_long <- down_df %>%
          rownames_to_column("Organ1") %>%
          pivot_longer(cols = -Organ1, names_to = "Organ2", values_to = "Count") %>%
          mutate(Species = sp, Count = as.numeric(Count))
        down_data_list[[sp]] <- down_long
      }
    }, error = function(e) {
      warning(paste("Error reading", down_file, ":", e$message))
    })
  } else {
    warning(paste("File not found:", down_file))
  }
}

# Combine all data
up_all <- bind_rows(up_data_list)
down_all <- bind_rows(down_data_list)

# Check if we have data
if (nrow(up_all) == 0 && nrow(down_all) == 0) {
  stop("No data was successfully loaded. Please check file paths and formats.")
}

# Create species to scientific name mapping
species_mapping <- c(
  "Ehy" = "Equisetum hyemale", "Pnu" = "Psilotum nudum", "Ore" = "Ophioglossum reticulatum",
  "Aev" = "Angiopteris evecta", "Aob" = "Abrodictyum obscurum", "Dcu" = "Dicranopteris curranii",
  "Lfl" = "Lygodium flexuosum", "Ala" = "Adiantum latifolium", "Sam" = "Salvinia molesta",
  "Cba" = "Cibotium barometz", "Als" = "Alsophila latebrosa", "Len" = "Lindsaea ensifolia",
  "Msp" = "Microlepia speluncae", "Pir" = "Pleocnemia irregularis",
  "Nbi" = "Nephrolepis biserrata", "Tin" = "Tectaria incisa", "Dde" = "Diplazium proliferum",
  "Ppi" = "Pyrrosia piloselloides", "Aop" = "Amblovenatum opulentum", "Dac" = "Davallia denticulata",
  "Spa" = "Stenochlaena palustris", "Adi" = "Actinostachys digitata"
)

# Create taxonomy information
taxonomy_info <- tribble(
  ~Species, ~Clade, ~Order, ~Family,
  "Ehy", "Eusporangiate", "Equisetales", "Equisetaceae",
  "Pnu", "Eusporangiate", "Psilotales", "Psilotaceae",
  "Aev", "Eusporangiate", "Marattiales", "Marattiaceae",
  "Ore", "Eusporangiate", "Ophioglossales", "Ophioglossaceae",
  "Aob", "Leptosporangiate", "Hymenophyllales", "Hymenophyllaceae",
  "Dcu", "Leptosporangiate", "Gleicheniales", "Gleicheniaceae",
  "Lfl", "Leptosporangiate", "Schizaeales", "Lygodiaceae",
  "Adi", "Leptosporangiate", "Schizaeales", "Schizaeaceae",
  "Sam", "Core leptosporangiate", "Salviniales", "Salviniaceae",
  "Cba", "Core leptosporangiate", "Cyatheales", "Cibotiaceae",
  "Als", "Core leptosporangiate", "Cyatheales", "Cyatheaceae",
  "Len", "Core leptosporangiate", "Polypodiales", "Lindsaeaceae",
  "Ala", "Core leptosporangiate", "Polypodiales", "Pteridaceae",
  "Msp", "Core leptosporangiate", "Polypodiales", "Dennstaedtiaceae",
  "Pir", "Core leptosporangiate", "Polypodiales", "Dryopteridaceae",
  "Nbi", "Core leptosporangiate", "Polypodiales", "Nephrolepidaceae",
  "Tin", "Core leptosporangiate", "Polypodiales", "Tectariaceae",
  "Dde", "Core leptosporangiate", "Polypodiales", "Athyriaceae",
  "Aop", "Core leptosporangiate", "Polypodiales", "Thelypteridaceae",
  "Ppi", "Core leptosporangiate", "Polypodiales", "Polypodiaceae",
  "Dac", "Core leptosporangiate", "Polypodiales", "Davalliaceae",
  "Spa", "Core leptosporangiate", "Polypodiales", "Blechnaceae"
)

# Calculate summaries (excluding self-comparisons)
up_summary <- up_all %>%
  filter(Organ1 != Organ2, !is.na(Count)) %>%
  group_by(Species, Organ = Organ1) %>%
  summarise(SPL_Up = sum(Count, na.rm = TRUE), .groups = "drop")

down_summary <- down_all %>%
  filter(Organ1 != Organ2, !is.na(Count)) %>%
  group_by(Species, Organ = Organ1) %>%
  summarise(SPL_Down = sum(Count, na.rm = TRUE), .groups = "drop")

# Combine summaries
combined_summary <- full_join(up_summary, down_summary, by = c("Species", "Organ")) %>%
  replace_na(list(SPL_Up = 0, SPL_Down = 0)) %>%
  mutate(Net_Change = SPL_Up - SPL_Down) %>%
  left_join(taxonomy_info, by = "Species") %>%
  mutate(Scientific_Name = species_mapping[Species])

# Remove rows with missing taxonomy or scientific names
combined_summary <- combined_summary %>%
  filter(!is.na(Scientific_Name), !is.na(Clade))

# Check if we have data for heatmap
if (nrow(combined_summary) == 0) {
  stop("No valid data available for heatmap generation after processing.")
}

# Prepare heatmap matrix
heatmap_data <- combined_summary %>%
  select(Scientific_Name, Organ, Net_Change) %>%
  pivot_wider(names_from = Organ, values_from = Net_Change, values_fill = 0)

# Convert to matrix
heatmap_matrix <- as.matrix(heatmap_data[, -1])
rownames(heatmap_matrix) <- heatmap_data$Scientific_Name

# Function to create and save heatmap with extended left margin
create_heatmap <- function(grouping_var, filename, title, display_only = FALSE) {
  
  # Create annotation dataframe
  annotation_data <- combined_summary %>%
    select(Scientific_Name, !!sym(grouping_var)) %>%
    distinct() %>%
    filter(Scientific_Name %in% rownames(heatmap_matrix)) %>%
    arrange(!!sym(grouping_var)) %>%
    column_to_rownames("Scientific_Name")
  
  colnames(annotation_data) <- grouping_var
  
  # Reorder matrix to match annotation order
  heatmap_matrix_ordered <- heatmap_matrix[rownames(annotation_data), , drop = FALSE]
  
  # Create gaps for visual separation
  annotation_groups <- annotation_data[[grouping_var]]
  gaps_positions <- which(diff(as.numeric(factor(annotation_groups))) != 0)
  
  # Set up color palette
  max_val <- max(abs(heatmap_matrix_ordered), na.rm = TRUE)
  if (max_val == 0) max_val <- 1  # Avoid division by zero
  
  colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
  breaks <- seq(-max_val, max_val, length.out = 101)
  
  # Calculate width with extra 2 cm (approximately 0.79 inches) for left margin
  base_width <- 7
  extended_width <- base_width + 0  # Adding ~2 cm
  
  # Create heatmap with extended width and adjusted cell spacing
  p <- pheatmap(
    heatmap_matrix_ordered,
    color = colors,
    breaks = breaks,
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    annotation_row = annotation_data,
    main = title,
    fontsize = 10,
    fontsize_row = 8,
    fontsize_col = 8,
    angle_col = 45,
    border_color = "grey90",
    annotation_names_row = TRUE,
    annotation_legend = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    gaps_row = gaps_positions,
    treeheight_row = 0,
    treeheight_col = 50,
    # Extended width to accommodate longer species names
    cellwidth = if(!display_only) 15 else NA,  # Slightly wider cells
    cellheight = if(!display_only) 12 else NA, # Consistent cell height
    filename = if(display_only) NA else filename,
    width = if(!display_only) extended_width else NULL,
    height = if(!display_only) 6 else NULL,
    dpi = if(!display_only) 300 else NULL
  )
  
  if (!display_only) {
    cat(paste("Saved:", filename, "with extended width:", extended_width, "inches\n"))
  }
  
  return(p)
}
# Function to create ggplot version
create_ggplot_heatmap <- function(grouping_var, filename, title) {
  
  plot_data <- combined_summary %>%
    select(Scientific_Name, Organ, Net_Change, !!sym(grouping_var)) %>%
    mutate(
      Scientific_Name = factor(Scientific_Name),
      Group = !!sym(grouping_var)
    ) %>%
    arrange(Group, Scientific_Name)
  
  # Reorder factor levels
  plot_data$Scientific_Name <- factor(plot_data$Scientific_Name, 
                                      levels = unique(plot_data$Scientific_Name))
  
  p <- ggplot(plot_data, aes(x = Organ, y = Scientific_Name, fill = Net_Change)) +
    geom_tile(color = "white", size = 0.1) +
    scale_fill_gradient2(
      low = "#2166AC", 
      mid = "white", 
      high = "#B2182B",
      midpoint = 0, 
      name = "Net Change"
    ) +
    facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 8),
      strip.text.y = element_text(angle = 0, size = 10),
      plot.title = element_text(hjust = 0.5, size = 12),
      legend.position = "right",
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey90")
    ) +
    labs(title = title, x = "Organ", y = "Species")
  
  ggsave(filename, plot = p, width = 12, height = 10, dpi = 300)
  cat(paste("Saved ggplot version:", filename, "\n"))
  
  return(p)
}

# Generate heatmaps
cat("Creating heatmaps...\n")

# Display in RStudio (optional)
if (interactive()) {
  cat("Displaying Clade heatmap...\n")
  create_heatmap("Clade", "", "Net SPL Gene Expression by Clade", display_only = TRUE)
}

# Save heatmaps
create_heatmap("Clade", "heatmap_by_clade.png", "Net SPL Gene Expression by Clade")
create_heatmap("Order", "heatmap_by_order.png", "Net SPL Gene Expression by Order")
create_heatmap("Family", "heatmap_by_family.png", "Net SPL Gene Expression by Family")

# Save ggplot versions
create_ggplot_heatmap("Clade", "heatmap_by_clade_ggplot.png", "Net SPL Gene Expression by Clade")
create_ggplot_heatmap("Order", "heatmap_by_order_ggplot.png", "Net SPL Gene Expression by Order")
create_ggplot_heatmap("Family", "heatmap_by_family_ggplot.png", "Net SPL Gene Expression by Family")

cat("All heatmap generation completed successfully!\n")

# Print summary statistics
cat("\nSummary Statistics:\n")
cat("Species processed:", length(unique(combined_summary$Species)), "\n")
cat("Organs found:", length(unique(combined_summary$Organ)), "\n")
cat("Total observations:", nrow(combined_summary), "\n")

create_clustered_heatmap_by_profile <- function() {
  # Prepare annotation (Order for visual grouping)
  annotation_data <- combined_summary %>%
    select(Scientific_Name, Order) %>%
    distinct() %>%
    filter(Scientific_Name %in% rownames(heatmap_matrix)) %>%
    column_to_rownames("Scientific_Name")
  
  colnames(annotation_data) <- "Order"
  
  # Match order to matrix
  heatmap_matrix_ordered <- heatmap_matrix[rownames(annotation_data), , drop = FALSE]
  
  # Set up color palette
  max_val <- max(abs(heatmap_matrix_ordered), na.rm = TRUE)
  if (max_val == 0) max_val <- 1
  
  colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
  breaks <- seq(-max_val, max_val, length.out = 101)
  
  # Create clustered heatmap
  pheatmap(
    heatmap_matrix_ordered,
    color = colors,
    breaks = breaks,
    cluster_rows = TRUE,       # Cluster by expression profile
    cluster_cols = TRUE,
    annotation_row = annotation_data,
    show_rownames = TRUE,
    show_colnames = TRUE,
    main = "SPL Net Expression (Clustered by Profile, Annotated by Order)",
    fontsize = 10,
    fontsize_row = 8,
    fontsize_col = 8,
    angle_col = 45,
    border_color = "grey90",
    filename = "heatmap_clustered_by_expression_profile.png",
    width = 10,
    height = 8,
    dpi = 300
  )
}
create_clustered_heatmap_by_profile()