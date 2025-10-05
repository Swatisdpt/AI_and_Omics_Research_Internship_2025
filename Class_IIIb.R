# Install and load required packages (only run installation once)
if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c("ArrayExpress","affy","arrayQualityMetrics"))

# Install CRAN packages for data manipulation
install.packages("dplyr")

# Load Required Libraries
library(ArrayExpress)             # Download GEO datasets (series matrix, raw CEL files)
library(affy)                 # Pre-processing of Affymetrix microarray data (RMA normalization)
library(arrayQualityMetrics)  # QC reports for microarray data
library(dplyr)                # Data manipulation
library(oligo)

# Download from ArrayExpress
#ae_data <- ArrayExpress("E-MTAB-6102")
#Pre-processed data not available for this dataset, so showing error.

# List all CEL files
cel_files <- list.celfiles( "C:/Users/swati/OneDrive/Documents/AI&Omics/AI_Omics_Internship_2025/E-MTAB-6102", full.names = TRUE)

# Read the raw data into an oligo FeatureSet
raw_data <- read.celfiles(cel_files)

# View a summary
raw_data

getwd()
setwd("C:/Users/swati/OneDrive/Documents/AI&Omics/AI_Omics_Internship_2025")
arrayQualityMetrics(expressionset = raw_data,
                    outdir = "Results/QC_Raw_Data",
                    force = TRUE,
                    do.logtransform = TRUE)
dir.create("Module_II")
file.rename("E-MTAB-6102", "Module_II/E-MTAB-6102")  # move CEL folder
file.rename("Results", "Module_II/Results")          # move QC results
file.rename("3B.R", "Module_II/3B.R")  # move R script

# RMA (Robust Multi-array Average) Normalization
normalized_data <- rma(raw_data)

# QC after data normalization 
arrayQualityMetrics(expressionset = normalized_data,
                    outdir = "Results/QC_Normalized_Data",
                    force = TRUE)

# Extract normalized expression values into a data frame
processed_data <- as.data.frame(exprs(normalized_data))

dim(processed_data)   # Dimensions: number of probes × number of samples

# Filter Low-Variance Transcripts (“soft” intensity based filtering)
# Calculate median intensity per probe across samples
row_median <- rowMedians(as.matrix(processed_data))

# 1. Perform quality control before and after normalization and  note down how many you found before and after normalization
# Ans. 16 out of 31 samples are found as outliers before, and 16 after normalization, 
# but outliers before were present in multiple detection methods whereas, after 
# normalization all the outliers were found only by MA plots.

# Visualize distribution of probe median intensities
hist(row_median,
     breaks = 100,
     freq = FALSE,
     main = "Median Intensity Distribution")

# Set a threshold to remove low variance probes (dataset-specific, adjust accordingly)
threshold <- 0.8
abline(v = threshold, col = "black", lwd = 2) 

# Select probes above threshold
indx <- row_median > threshold 
filtered_data <- processed_data[indx, ] 

nrow(processed_data)      # Before filtering
nrow(filtered_data)       # After filtering

# 2. Normalize the data and then apply filtering to remove low-intensity probes 
# and note how many transcripts remain. 
# 15734 before, and 9834 after filtering

# Overwrite processed data with filtered dataset
processed_data <- filtered_data 

#phenotype data preparation
sdrf_path <- "C:/Users/swati/OneDrive/Documents/AI&Omics/AI_Omics_Internship_2025/Module_II/E-MTAB-6102/E-MTAB-6102.sdrf.txt"

# Load the SDRF file into R as a data frame
sdrf_data <- read.delim(sdrf_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 3. Use the phenotype information to define your target groups and re-label them (e.g normal vs cancer)

# Replace processed_data column names with SDRF first column
colnames(processed_data) <- sdrf_data[, 1]

# Verify new column names
head(colnames(processed_data))

# Phenotype Data Preparation
class(sdrf_data$Factor.Value.disease.)


# Define experimental groups (normal vs cancer)
groups <- factor(sdrf_data$Factor.Value.disease.,
                 levels = c("normal", "bladder carcinoma"),
                 label = c("normal", "cancer"))

class(groups)
levels(groups)
