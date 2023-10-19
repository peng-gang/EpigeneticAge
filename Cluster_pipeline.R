library(readxl)
library(httr)
library(stringr)
library(data.table)
library(dplyr)
library(readxl)


Collect_CpG_sites <- function(){
  
  print("Collecting dataset information")
  datasets <- read_excel("/geode2/home/u030/harsurya/Carbonate/Methylation_Datasets/DNAMethylation/Human_Methylation_Datasets.xlsx",sheet = "All datasets")
  GSE <- datasets$`Accession Number`[datasets$Remarks == "Useful"]
  
  cumulative_samples <- data.table()
  cumulative_matrix <- data.table()
  

  for(i in 1:length(GSE)){
    
    print(paste0("Collecting CpG sites from dataset: ",GSE[i]))
    
    if(i==1){
      
      CpG_sites_all <- fread(paste0("/geode2/home/u030/harsurya/Carbonate/Methylation_Datasets/Processed_Data/",GSE[i],"/methylation_data.csv"), select = 1)
      
    }
    else{
      
      CpG_sites <- fread(paste0("/geode2/home/u030/harsurya/Carbonate/Methylation_Datasets/Processed_Data/",GSE[i],"/methylation_data.csv"), select = 1)
      CpG_sites_all <- rbind(CpG_sites_all,CpG_sites)
      
    }
    
    
  }
  
  colnames(CpG_sites_all) <- c("CpG_sites")
  
  print("All CpG sites info: ")
  
  print(head(CpG_sites_all))
  print(nrow(CpG_sites_all))
  
  CpG_sites_unique <- data.table()
  
  CpG_sites_unique$CpG_sites <- unique(CpG_sites_all$CpG_sites)
  
  print("Unique CpG sites info: ")
  
  print(head(CpG_sites_unique))
  print(nrow(CpG_sites_unique))
  
  return(CpG_sites_unique)
  
}


Cluster_CpG_sites <- function(CpG_sites){
  
  n = 100000
  CpG_sites_clusters <- data.table(replicate(0, numeric(0), n = n))
  start = 1
  end = 100000
  
  cluster_names <- c()
  for (i in 1:(as.integer(nrow(all_CpG_sites)/100000)+1)) {
    
    cluster_name <- paste0("cluster_",i)
    
    cluster_names <- c(cluster_names,cluster_name)
  }
  
  
  for(i in 1:(as.integer(nrow(all_CpG_sites)/100000)+1))
    
    if(i==(as.integer(nrow(all_CpG_sites)/100000)+1)){
      
      cluster_sites <- CpG_sites[start:end, ]
      CpG_sites_clusters[[cluster_names[i]]] <- cluster_sites$CpG_sites
      
    }else{
      
      cluster_sites <- CpG_sites[start:end, ]
      CpG_sites_clusters[[cluster_names[i]]] <- cluster_sites$CpG_sites
      start <- start + 100000
      end <- end + 100000
      
    }
  
  
  return(CpG_sites_clusters[,-1])
  
}


Find_correlation <- function(Clustered_CpG_sites){
  
  print("Collecting dataset information")
  datasets <- read_excel("/geode2/home/u030/harsurya/Carbonate/Methylation_Datasets/DNAMethylation/Human_Methylation_Datasets.xlsx",sheet = "All datasets")
  GSE <- datasets$`Accession Number`[datasets$Remarks == "Useful"]
  
for (j in 1:ncol(Clustered_CpG_sites)){
  
  print(paste0("Working with cluster: ",j))
  
  current_cluster <- Clustered_CpG_sites[, ..j]
  cluster_col_name <- paste0("cluster_",j)
  current_cluster <- current_cluster[[cluster_col_name]]
  
  if(j==ncol(Clustered_CpG_sites)){
    current_cluster <- na.omit(current_cluster)
  }
  
  print("Initializing merging of datasets")
  Initiated <- "No"
  for(i in 1:length(GSE)){
    
    print(paste0("Merging dataset: ", GSE[i]))
    matrix <- fread(paste0("/geode2/home/u030/harsurya/Carbonate/Methylation_Datasets/Processed_Data/",GSE[i],"/methylation_data.csv"))
    sample_data <-fread(paste0("/geode2/home/u030/harsurya/Carbonate/Methylation_Datasets/Processed_Data/",GSE[i],"/cleaned_metadata.csv"))
    row_names <- matrix$V1
    matrix<- matrix[, -1]  # Remove the first column from the data frame
    if(identical(colnames(matrix),sample_data$sample_id)) {print("Samples are in order")} else{print("Samples are not in order")}
    matrix <- as.matrix(matrix)
    rownames(matrix) <- row_names
    matrix <- t(matrix)
    matrix <- as.data.table(matrix)
    rownames(matrix) <- sample_data$sample_id
    
    missing_columns <- setdiff(current_cluster,colnames(matrix))
    
    print(paste0("No of Missing CpG sites in the dataset: ", length(missing_columns)))
    
    if(length(missing_columns) == 0){
      filtered_table <- matrix[,..current_cluster, with = FALSE]
      rownames(filtered_table) <- sample_data$sample_id
    }else if(length(missing_columns == length(current_cluster))){
      
      rm(matrix)
      rm(sample_data)
      next
      
    }else{
      current_cluster_temp <- current_cluster[!current_cluster %in% missing_columns]
      filtered_table <- matrix[,..current_cluster_temp, with = FALSE]
      for(k in 1:length(missing_columns)){
        filtered_table[[missing_columns[k]]] = NA
      }
      rownames(filtered_table) <- sample_data$sample_id
    }
    
    setcolorder(filtered_table, neworder = current_cluster)
    
    if(i==1 || Initiated == "No"){
      cumulative_filtered_table <- filtered_table
      cumulative_sample_data <- sample_data
      Initiated <- "Yes"
    }else{
      cumulative_filtered_table <- rbind(cumulative_filtered_table,filtered_table) 
      cumulative_sample_data <- bind_rows(cumulative_sample_data,sample_data)
    }
    
    rm(filtered_table)
    rm(matrix)
    rm(sample_data)
    
  }
  
  
  print("Finished merging")
  print("Removing sites with less than 80% data")
  
  na_counts <- colSums(is.na(cumulative_filtered_table))
  threshold <- 0.2 * nrow(cumulative_filtered_table)
  columns_to_remove <- which(na_counts > threshold)
  cumulative_filtered_table <- cumulative_filtered_table[, -columns_to_remove, with = FALSE]
  
  print("Calculating correlation for all sites with adequate data")
  
  for(m in 1:ncol(cumulative_filtered_table)){
    
    cor <- cor.test(cumulative_filtered_table[[m]],cumulative_sample_data$age)
    cor_df <- data.table(colnames(cumulative_filtered_table)[m],cor$statistic,cor$estimate,cor$p.value,cor$method)
    colnames(cor_df) <- c("CpG_site","Statistic","Estimate","pValue","Method")
    
    if (m == 1){
      cum_cor_df <- cor_df
    }else{
      cum_cor_df <- rbind(cum_cor_df,cor_df)
    }
  }
  
  if(j==1){
    
    total_cum_cor_df <- cum_cor_df
    
  }else{
    
    total_cum_cor_df <- rbind(total_cum_cor_df,cum_cor_df)
    
  }
  
}
  
  
  
  return(total_cum_cor_df)
  
}



all_CpG_sites <- Collect_CpG_sites()


Clustered_CpG_sites <- Cluster_CpG_sites(all_CpG_sites)


Correlation_data <- Find_correlation(Clustered_CpG_sites)


write.csv(Correlation_data, file = "/geode2/home/u030/harsurya/Carbonate/Methylation_Datasets/Correlation_data.csv",row.names = FALSE)

