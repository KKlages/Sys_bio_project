# --- Load and clean organ mapping once ---
library(readxl)
organ_mapping <- read_excel("organ_type.xlsx")

# Clean the organ names to match column names in counts
organ_mapping$Clean_Tissue <- gsub("[^a-z0-9]+", "_", tolower(organ_mapping$`Organ name`))

manual_add <- data.frame(
  `Organ type` = c("Roots", "Leaves", "Leaves"),  # add more if needed
  Clean_Tissue = c("root", "sterile_leaflet", "fertile_leaflet")
)

# Assuming organ_mapping has Clean_Tissue and Organ type columns
organ_mapping <- rbind(
  organ_mapping[, c("Organ type", "Clean_Tissue")],
  manual_add
)

# Then recreate mapping_dict
mapping_dict <- setNames(organ_mapping$`Organ type`, organ_mapping$Clean_Tissue)

missing_tissues <- unique(standardized_tissues[is.na(mapping_dict[standardized_tissues])])
if(length(missing_tissues) > 0) {
  cat("Missing mappings for tissues:\n")
  print(missing_tissues)
}
