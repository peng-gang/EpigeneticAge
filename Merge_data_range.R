library(readxl)
library(httr)
library(stringr)
library(data.table)
library(dplyr)
library(readxl)



between <- function(data,column,min,max){
  data <- data[data[[column]] >= min]
  data <- data[data[[column]] <= max]
  return(data)
}



Merge_samples_in_range <- function(min,max){
  
  print("Collecting dataset information")
  datasets <- read_excel("/geode2/home/u030/harsurya/Carbonate/Methylation Datasets/DNAMethylation/Human_Methylation_Datasets.xlsx",sheet = "All datasets")
  GSE <- datasets$`Accession Number`[datasets$Remarks == "Useful"]
  
  cumulative_samples <- data.table()
  cumulative_matrix <- data.table()
  
  print("Initializing merging of datasets")
  for(i in 1:length(GSE)){
    
    print(paste0("Merging dataset: ",GSE[i]))
    matrix <- fread(paste0("/geode2/home/u030/harsurya/Carbonate/Methylation Datasets/Processed_Data/",GSE[i],"/methylation_data.csv"))
    sample_data <-fread(paste0("/geode2/home/u030/harsurya/Carbonate/Methylation Datasets/Processed_Data/",GSE[i],"/cleaned_metadata.csv"))
    row_names <- matrix$V1
    matrix<- matrix[, -1]  # Remove the first column from the data frame
    if(identical(colnames(matrix),sample_data$sample_id)) {print("Samples are in order")} else{print("Samples are not in order")}
    matrix <- as.matrix(matrix)
    rownames(matrix) <- row_names
    matrix <- t(matrix)
    matrix <- as.data.table(matrix)
    rownames(matrix) <- sample_data$sample_id
    filtered_samples <- between(sample_data,"age",min,max)
    matrix <- matrix[rownames(matrix) %in% filtered_samples$sample_id]
    
    if(i==1){
      
      cumulative_samples <- filtered_samples
      cumulative_matrix <- matrix
      rownames(cumulative_matrix) <- cumulative_samples$sample_id
    }
    
    cumulative_samples <- bind_rows(cumulative_samples,filtered_samples)
    cumulative_matrix <- merge(cumulative_matrix,matrix, all = TRUE, sort = FALSE)
    rownames(cumulative_matrix) <- cumulative_samples$sample_id
    
  }
  
  print("Finished merging")
  
  return_tables <- list(table1 = cumulative_matrix, table2 = cumulative_samples)
  return(return_tables)
  
}

dataset_information <- Merge_samples_in_range(1,10)

methylation_data <- dataset_information$table1
sample_data <- dataset_information$table2


