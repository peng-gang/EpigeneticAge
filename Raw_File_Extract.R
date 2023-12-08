#Use this script to extract methylation profiling data from raw data file as some times the series matrix file might not have the beta values in it

#Load the required libraries
#BiocManager::install("minfi")
#BiocManager::install("IlluminaHumanMethylationEPICmanifest")
library(minfi)
library(IlluminaHumanMethylationEPICmanifest)


# Set your path to the .idat files of the dataset
idatPath <- "/geode2/home/u030/harsurya/Carbonate/Methylation_Datasets/DNAMethylation/GSE179849_RAW/"

# Retrieve all the compressed idat files in the set idatpath
gz_files <- list.files(path = idatPath, pattern = "\\.gz$", full.names = TRUE)

# Decompress files in the same path and delete the compressed file
for (gz_file in gz_files) {
  system(paste("gzip -d", shQuote(gz_file)))
  # Delete compressed file
  unlink(gz_file)
}

# Load IDAT files
rgSet <- read.metharray.exp(base = idatPath, recursive = TRUE)


# Preprocess the data
preprocessIllumina(rgSet)


# Extract beta values
beta_values <- getBeta(rgSet)


#Transpose the betavalues matrix and convert it into a data frame
beta_values_df <- data.frame(t(beta_values))


#Use the data frame for further analysis
