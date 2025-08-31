#Assignment_2

# Create common subfolders for project organization
dir.create("Raw_data")     
dir.create("Script")   
dir.create("Results") 

# Define input and output folders
input_dir <- "Raw_data" 
output_dir <- "Results"

# Create output folder if it doesn't exist
if(!dir.exists(output_dir)){
  dir.create(output_dir)
}

# Import Data from CSV

data <- read.csv(file.choose())

# List of files to process
files_to_process <- c("DEGs_Data_1.csv", "DEGs_Data_2.csv")

# Prepare empty list to store results
result_list <- list()

# Define gene classification function
classify_gene <- function(logFC, padj) {
  if (logFC > 1 & padj < 0.05) {
    return("Upregulated")
  } else if (logFC < -1 & padj < 0.05) {
    return("Downregulated")
  } else {
    return("Not_Significant")
  }
}

# Loop through each file
for (file_name in files_to_process) {
  
  cat("\nProcessing:", file_name, "\n")
  
  input_file_path <- file.path(input_dir, file_name)
  
  # Import dataset
  data <- read.csv(input_file_path, stringsAsFactors = FALSE)
  cat("File imported. Checking for missing values...\n")
  
  # Replace missing padj values with 1
  if ("padj" %in% names(data)) {
    na_count <- sum(is.na(data$padj))
    cat("Missing padj values:", na_count, "\n")
    data$padj[is.na(data$padj)] <- 1
  }
  
  # Rename if needed
  if ("log2FoldChange" %in% names(data)) {
    data$logFC <- data$log2FoldChange
  }
  
  # Ensure numeric
  data$logFC <- as.numeric(data$logFC)
  data$padj <- as.numeric(data$padj)
  
  # Create empty status column
  data$status <- character(nrow(data))
  
  # Classify each gene using for loop
  for (i in 1:nrow(data)) {
    data$status[i] <- classify_gene(data$logFC[i], data$padj[i])
  }
  
  cat("Classification completed.\n")
  
  # Save results to list in R
  result_list[[file_name]] <- data
  
  # Save processed file to Results folder
  output_file_path <- file.path(output_dir, paste0("Processed_", file_name))
  write.csv(data, output_file_path, row.names = FALSE)
  cat("Results saved to:", output_file_path, "\n")
  
  # Print summary
  cat("Summary for", file_name, ":\n")
  print(table(data$status))
  cat("---------------------------\n")
}

# Access individual result sets
results_1 <- result_list[["DEGs_Data_1.csv"]]
results_2 <- result_list[["DEGs_Data_2.csv"]]

save.image(file = "Swati_Sudipta_Sahoo_Class_2_Assignment.RData")

