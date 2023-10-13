library(readxl)
library(httr)
library(stringr)
library(data.table)
library(dplyr)
library(readxl)

datasets <- read_excel("/geode2/home/u030/harsurya/Carbonate/Methylation Datasets/DNAMethylation/Human_Methylation_Datasets.xlsx",sheet = "All datasets")

GSE <- datasets$`Accession Number`[datasets$Remarks == "Useful"]

# GSE <- c("GSE101764","GSE106648","GSE114134","GSE124076","GSE124366")

#GSE114135
#GSE146377

for (i in 1:length(GSE)){
  print(paste0("Merging dataset: ",GSE[i]))
  if(i == 1 ){
    matrix_cum <- fread(paste0("/geode2/home/u030/harsurya/Carbonate/Methylation Datasets/Processed_Data/",GSE[i],"/methylation_data.csv"))
    sample_data_cum <-fread(paste0("/geode2/home/u030/harsurya/Carbonate/Methylation Datasets/Processed_Data/",GSE[i],"/cleaned_metadata.csv"))
    row_names <- matrix_cum$V1
    matrix_cum<- matrix_cum[, -1]  # Remove the first column from the data frame
    if(identical(colnames(matrix_cum),sample_data_cum$sample_id)) {print("Samples are in order")} else{print("Samples are not in order")}
    matrix_cum <- as.matrix(matrix_cum)
    rownames(matrix_cum) <- row_names
    matrix_cum <- t(matrix_cum)
    matrix_cum <- as.data.table(matrix_cum)
    rownames(matrix_cum) <- sample_data_cum$sample_id
  }
  
  else{
    
    matrix <- fread(paste0("/geode2/home/u030/harsurya/Carbonate/Methylation Datasets/Processed_Data/",GSE[i],"/methylation_data.csv"))
    sample_data <-fread(paste0("/geode2/home/u030/harsurya/Carbonate/Methylation Datasets/Processed_Data/",GSE[i],"/cleaned_metadata.csv"))
    row_names <- matrix$V1
    matrix<- matrix[, -1]  # Remove the first column from the data frame
    if(identical(colnames(matrix),sample_data$sample_id)) {print("Samples are in order")} else{print("Samples are not in order")}
    matrix <- as.matrix(matrix)
    rownames(matrix) <- row_names
    matrix<- t(matrix)
    matrix <- as.data.table(matrix)
    rownames(matrix) <- sample_data$sample_id
    
    matrix_cum <- merge(matrix_cum,matrix,all = TRUE, sort = FALSE)
    rownames(matrix_cum) <- c(sample_data_cum$sample_id,sample_data$sample_id)
    sample_data_cum <- bind_rows(sample_data_cum,sample_data)
    
  }
}

# matrix_cum <- fread(paste0("/geode2/home/u030/harsurya/Carbonate/Methylation Datasets/Processed_Data/",GSE,"/methylation_data.csv"))
# row_names <- matrix_cum$V1
# matrix_cum<- matrix_cum[, -1]  # Remove the first column from the data frame
# matrix_cum <- as.matrix(matrix_cum)
# rownames(matrix_cum) <- row_names
# matrix_cum <- t(matrix_cum)
# matrix_cum <- as.data.table(matrix_cum)
# sample_data_cum <-fread(paste0("/geode2/home/u030/harsurya/Carbonate/Methylation Datasets/Processed_Data/",GSE,"/cleaned_metadata.csv"))
# rownames(matrix_cum) <- sample_data_cum$sample_id
# 
# 
# matrix_cum$V1[1]
# matrix_cum[1,]
# 
# matrix_cum<- matrix_cum[, -1]  # Remove the first column from the data frame
# head(matrix)
# 
# file <- paste0("/geode2/home/u030/harsurya/Carbonate/Methylation Datasets/Processed_Data/",GSE,"/matrix.txt.gz")
# gse <- getGEO(filename = file,GSEMatrix=TRUE)
# sample_data <- read.csv(paste0("/geode2/home/u030/harsurya/Carbonate/Methylation Datasets/Processed_Data/",GSE,"/cleaned_metadata.csv"))
# matrix <- gse@assayData[["exprs"]]
# matrix <- t(as.matrix(matrix))
# matrix<- as.data.table(matrix)
# 
# 
# 
# 
# 
# GSE <- "GSE106648"
# 
# file <- paste0("/geode2/home/u030/harsurya/Carbonate/Methylation Datasets/Processed_Data/",GSE,"/matrix.txt.gz")
# gse <- getGEO(filename = file,GSEMatrix=TRUE)
# sample_data2 <- read.csv(paste0("/geode2/home/u030/harsurya/Carbonate/Methylation Datasets/Processed_Data/",GSE,"/cleaned_metadata.csv"))
# matrix2 <- gse@assayData[["exprs"]]
# matrix2 <- t(as.matrix(matrix2))
# matrix2<- as.data.table(matrix2)
# 
# 
# 
# 
# 
# 
# 
# GSE <- "GSE114134"
# 
# file <- paste0("/geode2/home/u030/harsurya/Carbonate/Methylation Datasets/Processed_Data/",GSE,"/matrix.txt.gz")
# gse <- getGEO(filename = file,GSEMatrix=TRUE)
# sample_data2 <- read.csv(paste0("/geode2/home/u030/harsurya/Carbonate/Methylation Datasets/Processed_Data/",GSE,"/cleaned_metadata.csv"))
# matrix2 <- gse@assayData[["exprs"]]
# matrix2 <- t(as.matrix(matrix2))
# matrix2<- as.data.table(matrix2)
# 
# 
# matrix <- merge(matrix,matrix2,all = TRUE, sort = FALSE)
# rownames(matrix) <- c(sample_data$sample_id,sample_data2$sample_id)
# sample_data <- bind_rows(sample_data,sample_data2)
# 
# 
# GSE <- "GSE114135"
# 
# file <- paste0("/geode2/home/u030/harsurya/Carbonate/Methylation Datasets/Processed_Data/",GSE,"/matrix.txt.gz")
# gse <- getGEO(filename = file,GSEMatrix=TRUE)
# sample_data2 <- read.csv(paste0("/geode2/home/u030/harsurya/Carbonate/Methylation Datasets/Processed_Data/",GSE,"/cleaned_metadata.csv"))
# matrix2 <- gse@assayData[["exprs"]]
# matrix2 <- t(as.matrix(matrix2))
# matrix2<- as.data.table(matrix2)
# 
# 
# matrix <- merge(matrix,matrix2,all = TRUE, sort = FALSE)
# rownames(matrix) <- c(sample_data$sample_id,sample_data2$sample_id)
# sample_data <- bind_rows(sample_data,sample_data2)

print(identical(rownames(matrix_cum),sample_data_cum$sample_id))

library(ggplot2)

plot_df <- data.frame(matrix_cum$cg00002473,sample_data_cum$age)
colnames(plot_df) <- c("Beta_values","age")



ggplot(plot_df,aes(age,Beta_values)) +
  geom_smooth(se = FALSE) +
  theme_classic()
