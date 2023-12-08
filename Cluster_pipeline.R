#The code in this script will be useful to gather all unique CpG sites from all the available datasets and cluster the sites with desired size

#Load required libraries
library(readxl)
library(httr)
library(stringr)
library(data.table)
library(dplyr)
library(readxl)


#Function to collect all the unique CpG sites from all the available datasets
Collect_CpG_sites <- function(){
  
  print("Collecting dataset information")
  
  #Read the information of datasets from the Human_Methylation_Datasets.xlsx and sheet All datasets. Each dataset has the following variables,
  #Accession Number - GEO accession number of the dataset
  #Database	- Database from where the dataset has been downloaded
  #Remarks - Each dataset has one of the following remarks
    ##Invalid cpg sites - These datasets do not have valid CpG sites
    ##Useful - These datasets are downloaded and processed and can be used for analysis
    ##Redundant data - Datasets with redundant samples
    ##File read error - Not able to read the datasets using getGEO() function
    ##Empty matrix - Series matrix file does not include the methylation data - have to download raw file and process it using Raw_File_Extract.R
    ##27K Beadchip - Datasets are generated using 27K bead chip platform
  #Link	- Link to the GEO page of the dataset
  #Citations - Citations included in the dataset
  #Tissue Type	- Sample tissue type used in the dataset
  #Sample Size - Number of samples in the dataset
  #ftp link - Link to download the series matrix file of the dataset	
  #Min age - Minimum age of the samples
  #Max age - Maximum age of the samples
  #Cancer dataset - Whether the dataset is cancer dataset or not
  #Bead Chip - Platform used for methylation profiling
  datasets <- read_excel("/geode2/home/u030/harsurya/Carbonate/Methylation_Datasets/DNAMethylation/Human_Methylation_Datasets.xlsx",sheet = "All datasets")
  
  #Get the accession numbers of the datasets which are marked as useful into the variable GSE
  GSE <- datasets$`Accession Number`[datasets$Remarks == "Useful"]
  
  #Create empty datasets to store the data
  cumulative_samples <- data.table()
  cumulative_matrix <- data.table()
  
  #For loop to iterate through GSE to get the unique CpG sites from all the datasets one at a time
  for(i in 1:length(GSE)){
    
    #Print the accession number of the current dataset from which the CpG sites are being collected
    print(paste0("Collecting CpG sites from dataset: ",GSE[i]))
    
    #If else block to decide whether to append the data to the cumulative data - 1st iteration does not need to append the data
    if(i==1){
      
      #Read only the first column of the methylation dataset as it has all the IDs of CpG sites
      CpG_sites_all <- fread(paste0("/geode2/home/u030/harsurya/Carbonate/Methylation_Datasets/Processed_Data/",GSE[i],"/methylation_data.csv"), select = 1)
      
    }
    else{
      
      #Read and append the CpG sites to the already read data
      CpG_sites <- fread(paste0("/geode2/home/u030/harsurya/Carbonate/Methylation_Datasets/Processed_Data/",GSE[i],"/methylation_data.csv"), select = 1)
      CpG_sites_all <- rbind(CpG_sites_all,CpG_sites)
      
    }
    
    
  }
  
  #Assign column name to the CpG sites data table
  colnames(CpG_sites_all) <- c("CpG_sites")
  
  #Printing Head and Tail of the CpG sites data table
  print("All CpG sites info: ")
  
  print(head(CpG_sites_all))
  print(nrow(CpG_sites_all))
  
  #Create a data table and collect unique CpG sites in it
  CpG_sites_unique <- data.table()
  CpG_sites_unique$CpG_sites <- unique(CpG_sites_all$CpG_sites)
  
  #Printing Head and Tail of the CpG sites unique data table to compare with previous data table which may contain redundant CpG sites
  print("Unique CpG sites info: ")
  
  print(head(CpG_sites_unique))
  print(nrow(CpG_sites_unique))
  
  #Return the Unique CpG sites
  return(CpG_sites_unique)
  
}


#Function to group CpG sites into different clusters
#Takes the input of all Unique CpG sites that is obtained from the function Collect_CpG_sites in this script
Cluster_CpG_sites <- function(CpG_sites){
  
  #Number of CpG sites to assign to each cluster. n is set to 10,000 CpG sites per cluster. Change it for your convenience by keeping in mind of the computing resources
  n = 10000
  
  #Create an empty data table with the size of n rows and 1 column with default value of zero for all rows
  #Use this datatable to store the CpG sites that are clustered into one cluster
  CpG_sites_clusters <- data.table(replicate(0, numeric(0), n = n))
  
  #Define the range of CpG sites to cluster in current iteration
  #In first iteration it will be 1 to 10,000 (default n) CpG sites of all the unique CpG sites
  start = 1
  end = n
  
  #Cluster names vector to store the names of each cluster 
  cluster_names <- c()
  
  #For loop to generate the cluster names for those many clusters as can be obtained by dividing the total unique CpG sites by the size of n
  for (i in 1:(as.integer(nrow(all_CpG_sites)/n)+1)) {
    
    #Generate each cluster name and append it to the vector
    cluster_name <- paste0("cluster_",i)
    cluster_names <- c(cluster_names,cluster_name)
    
  }
  
  #For loop to divide all the unique CpG sites into clusters
  for(i in 1:(as.integer(nrow(all_CpG_sites)/n)+1))
    
    #If else block to decide if the current iteration is the last one or not
    #If block is executed only for the last iteration
    if(i==(as.integer(nrow(all_CpG_sites)/n)+1)){
      
      #In the last iteration the remaining CpG sites are clustered together and added as a new column to CpG_sites_clusters data table
      cluster_sites <- CpG_sites[start:end, ]
      CpG_sites_clusters[[cluster_names[i]]] <- cluster_sites$CpG_sites
      
    }else{
      
      #In all other iterations the Cpg sites of the current iteration are added as a new column to the CpG_sites_clusters data table and the next clusters range is updated
      cluster_sites <- CpG_sites[start:end, ]
      CpG_sites_clusters[[cluster_names[i]]] <- cluster_sites$CpG_sites
      start <- start + n
      end <- end + n
      
    }
  
  
  #Return the data table that contains the clusters of CpG sites with each column containing CpG sites of one cluster
  return(CpG_sites_clusters[,-1])
  
  #Note: The data table only contains the CpG sites IDs present in the each cluster. This clustered data is used to actually cluster and save the methylation data in the function Find_correlation in this script
}



#Function to calculate correlation between Beta values of CpG sites and age
#The function also combines methylation data cluster wise from all the available datasets 
#Takes the input of clusteres CpG sites as returned by the function Cluster_CpG_sites in this script
Find_correlation <- function(Clustered_CpG_sites){
  
  
  print("Collecting dataset information")
  
  #Read the information of datasets from the Human_Methylation_Datasets.xlsx and sheet All datasets. Each dataset has the following variables,
  #Accession Number - GEO accession number of the dataset
  #Database	- Database from where the dataset has been downloaded
  #Remarks - Each dataset has one of the following remarks
  ##Invalid cpg sites - These datasets do not have valid CpG sites
  ##Useful - These datasets are downloaded and processed and can be used for analysis
  ##Redundant data - Datasets with redundant samples
  ##File read error - Not able to read the datasets using getGEO() function
  ##Empty matrix - Series matrix file does not include the methylation data - have to download raw file and process it using Raw_File_Extract.R
  ##27K Beadchip - Datasets are generated using 27K bead chip platform
  #Link	- Link to the GEO page of the dataset
  #Citations - Citations included in the dataset
  #Tissue Type	- Sample tissue type used in the dataset
  #Sample Size - Number of samples in the dataset
  #ftp link - Link to download the series matrix file of the dataset	
  #Min age - Minimum age of the samples
  #Max age - Maximum age of the samples
  #Cancer dataset - Whether the dataset is cancer dataset or not
  #Bead Chip - Platform used for methylation profiling
  datasets <- read_excel("/geode2/home/u030/harsurya/Carbonate/Methylation_Datasets/DNAMethylation/Human_Methylation_Datasets.xlsx",sheet = "All datasets")
  #Get the accession numbers of the datasets which are marked as useful into the variable GSE
  GSE <- datasets$`Accession Number`[datasets$Remarks == "Useful"]

#For loop to iterate through CpG sites of each cluster and gather methylation data for CpG sites in each cluster from all datasets and then save it to the local computer  
for (j in 1:ncol(Clustered_CpG_sites)){
  
  #Print the cluster number for whose CpG sites the data is being collected in this iteration
  print(paste0("Working with cluster: ",j))
  
  #Extract the CpG sites from the current cluster
  current_cluster <- Clustered_CpG_sites[, ..j]
  cluster_col_name <- paste0("cluster_",j)
  current_cluster <- current_cluster[[cluster_col_name]]
  
  #If the iteration is processing the last cluster omit NAs as the last cluster may not contain all 10000 CpG sites
  if(j==ncol(Clustered_CpG_sites)){
    current_cluster <- na.omit(current_cluster)
  }
  
  #Collect methylation data for all CpG sites in the current cluster from all the available datasets
  print("Initializing merging of datasets")
  #A boolean variable Initiated is used to keep track of the data merge status - whether it is strated or not
  Initiated <- "No"
  #For loop to iterate through each dataset and collect the methylation data
  for(i in 1:length(GSE)){
    
    #Print the current dataset being used to collect methylation data
    print(paste0("Merging dataset: ", GSE[i]))
    
    #Read the methylation data of the current dataset from its folder
    #The code assumes the following data organisation structure
    ##|--Processed_Data
    #   |
    #   |
    #   |---All the datasets with each dataset information in a separete folder that is named afer its GEO accession number
    #       |
    #       |
    #       |---Each dataset folder has the following four files
    #           matrix.txt.gz - The series matrix file of the particular dataset - downloaded using Data_Cleaning.R
    #           methylation_data.csv - a csv file containing the beta values for the dataset - Extracted using Data_Cleaning.R
    #           metadata.csv - a csv file containing the sample information - Extracted using Data_Cleaning.R
    #           cleaned_metadata.csv - a csv file containing formatted sample information mainly containing sample id and age  - Extracted from metadata.csv using Data_Cleaning.R
    
    #Read the methylation data from the current dataset from methylation_data.csv
    matrix <- fread(paste0("/geode2/home/u030/harsurya/Carbonate/Methylation_Datasets/Processed_Data/",GSE[i],"/methylation_data.csv"))
    
    #Read the sample data for the current dataset using cleaned_metadata.csv
    sample_data <-fread(paste0("/geode2/home/u030/harsurya/Carbonate/Methylation_Datasets/Processed_Data/",GSE[i],"/cleaned_metadata.csv"))
    
    #Save the row names of the matrix i.e., the names of the CpG sites
    row_names <- matrix$V1
    
    #The first column contains the CpG site Ids, it can be removed as they are saved in the row_names
    matrix<- matrix[, -1]  # Remove the first column from the data frame
    
    #Check if the samples in the matrix and in the sample data are in same order or not
    if(identical(colnames(matrix),sample_data$sample_id)) {print("Samples are in order")} else{print("Samples are not in order")}
    
    #The beta values matrix is read as data table - convert it to the matrix
    matrix <- as.matrix(matrix)
    
    #Assign the rownames as they are lost after conversion to matrix
    rownames(matrix) <- row_names
    
    #Transpose the matrix to make CpG sites as columns and samples as rows  - will be easy for correlation
    matrix <- t(matrix)
    
    #Convert the matrix into data table and assign the sample ids as row names
    matrix <- as.data.table(matrix)
    rownames(matrix) <- sample_data$sample_id
    
    #The current dataset may or may not contain all the CpG sites present in the current cluster - getting the  missing CpG sites into missing_columns
    missing_columns <- setdiff(current_cluster,colnames(matrix))
    print(paste0("No of Missing CpG sites in the dataset: ", length(missing_columns)))
    
    #If else block to gather methylation data based on number of missing columns
    #If the number of missing columns is 0 it means all the CpG sites in the current cluster has data in the current dataset
    if(length(missing_columns) == 0){
      filtered_table <- matrix[,..current_cluster, with = FALSE]
      rownames(filtered_table) <- sample_data$sample_id
      
    #If all the CpG sites have no data in the current dataset - ignore the dataset for this cluster
    }else if(length(missing_columns == length(current_cluster))){
      
      rm(matrix)
      rm(sample_data)
      next
      
    #If there are certain CpG sites that do not have data in the current dataset assign NA for beta values for all the samples in this dataset for these CpG sites
    }else{
      
      #Retrieve those CpG sites which have data in the current dataset
      current_cluster_temp <- current_cluster[!current_cluster %in% missing_columns]
      filtered_table <- matrix[,..current_cluster_temp, with = FALSE]
      
      #For all other CpG sites which have no data in this dataset assign NA values
      for(k in 1:length(missing_columns)){
        filtered_table[[missing_columns[k]]] = NA
      }
      rownames(filtered_table) <- sample_data$sample_id
    }
    
    #Order the columns based on the current_cluster vector
    setcolorder(filtered_table, neworder = current_cluster)
    
    #If else block to merge the data based on status of the data merge
    #If the dataset being processed in this iteration is the first one, Append the data to the cumulative table and change Initiated to Yes
    if(i==1 || Initiated == "No"){
      #Assign the filtered data in this cluster to the cumulative data table
      cumulative_filtered_table <- filtered_table
      
      #Add an extra column GEO accesion to the cumulative table to keep track of the dataset to which a particular sample belongs to
      GEO_Accesion <- rep(c(GSE[i]),times = nrow(sample_data))
      GEO_Accesion <- data.frame(GEO_Accesion)
      sample_data$GEO_Accesion <- GEO_Accesion$GEO_Accesion
      
      #Assign the sample data to the cumulative sample data
      cumulative_sample_data <- sample_data
      
      #Change the initiated status to Yes
      Initiated <- "Yes"
      
    #If the current iteration is not the first dataset that is being processed  
    }else{
      
      #Append the current dataset data to the cumulative table
      cumulative_filtered_table <- rbind(cumulative_filtered_table,filtered_table)
      
      #Add an extra column GEO accesion to the cumulative table to keep track of the dataset to which a particular sample belongs to
      GEO_Accesion <- rep(c(GSE[i]),times = nrow(sample_data))
      GEO_Accesion <- data.frame(GEO_Accesion)
      sample_data$GEO_Accesion <- GEO_Accesion$GEO_Accesion
      
      #Append the current sample data to the cumulative sample data
      cumulative_sample_data <- bind_rows(cumulative_sample_data,sample_data)
    }
    
    #Remove the current dataset tables as the data is already appended to the cumulative tables
    rm(filtered_table)
    rm(matrix)
    rm(sample_data)
    
  }
  
  
  print("Finished merging")
  print("Removing sites with less than 80% data")
  
  #We only want to keep those CpG sites that have beta values for atleast 80% of the total samples
  #Since for some CpG sites NA values are added if they do not have data for a particular dataset,
  #They might not have beta values for more than 80% of the data. Such CpG sites must be discarded
  #Get the number of NA counts for each CpG site
  na_counts <- colSums(is.na(cumulative_filtered_table))
  
  #Set the threshold that the number of NAs for any CpG site shall not be over 20%
  threshold <- 0.2 * nrow(cumulative_filtered_table)
  
  #Retrieve the CpG sites that have less than 80% of the data
  columns_to_remove <- which(na_counts > threshold)
  
  #Remove these columns from the cumulative dataset
  cumulative_filtered_table <- cumulative_filtered_table[, -columns_to_remove, with = FALSE]
  
  #If certain samples have age as NA, then we do not need such samples information as we cannot correlate the sample information with age
  print("Removing samples with age = NA")
  
  #Retrieve the samples that age as NA
  samples_to_remove <- is.na(cumulative_sample_data$age)
  
  #Remove these samples from cumulative data tables
  cumulative_filtered_table <- cumulative_filtered_table[!samples_to_remove, ]
  cumulative_sample_data <- cumulative_sample_data[!samples_to_remove, ]
  
  #Update the sample names
  rownames(cumulative_filtered_table) <- cumulative_sample_data$sample_id
  
  #Removing those samples with age > 100 to avoid collecting samples where age is in days
  print("Removing samples with age > 100")
  
  #Retrieve the samples that has age > 100
  samples_to_remove <- cumulative_sample_data$age > 100
  
  #Remove these samples from cumulative data tables
  cumulative_filtered_table <- cumulative_filtered_table[!samples_to_remove, ]
  cumulative_sample_data <- cumulative_sample_data[!samples_to_remove, ]
  
  #Update the sample names
  rownames(cumulative_filtered_table) <- cumulative_sample_data$sample_id
  
  #Save the Methylation data and sample data of the current cluster
  write.csv(cumulative_filtered_table, file = paste0("/geode2/home/u030/harsurya/Carbonate/Clusters/cluster_",j,"_beta_values.csv"))
  write.csv(cumulative_sample_data, file = paste0("/geode2/home/u030/harsurya/Carbonate/Clusters/cluster_",j,"_sample_data.csv"))
  
  #Calculate the correlation coefficient for all CpG sites between their beta values and  sample age
  print("Calculating correlation for all sites with adequate data")
  
  #For loop to calculate correlation coefficient for all CpG sites between their beta values and  sample age
  for(m in 1:ncol(cumulative_filtered_table)){
    
    #Use cor.test function to calculate the correlation and save the results for each CpG site in each iteration
    cor <- cor.test(cumulative_filtered_table[[m]],cumulative_sample_data$age)
    
    #Save the required statistics in a data table and assign the appropriate column names
    cor_df <- data.table(colnames(cumulative_filtered_table)[m],cor$statistic,cor$estimate,cor$p.value,cor$method)
    colnames(cor_df) <- c("CpG_site","Statistic","Estimate","pValue","Method")
    
    #If this is the first iteration no data binding is required
    if (m == 1){
      cum_cor_df <- cor_df
      
    #If this is not first iteration bind the new row to the existing data table  
    }else{
      cum_cor_df <- rbind(cum_cor_df,cor_df)
    }
  }
  
  #If else block to bind correlation results of CpG sites from each cluster into a single data table
  #If the current iteration is processing the first cluster
  if(j==1){
    
    #Assign the correlation information to the cumulative table
    total_cum_cor_df <- cum_cor_df
    
  #If the current iteration is not processing the first cluster
  }else{
    
    #Append the correlation information of current cluster to the cumulative data table
    total_cum_cor_df <- rbind(total_cum_cor_df,cum_cor_df)
    
  }
  
}
  
  
  #Return the correlation data for all CpG sites from all the clusters
  return(total_cum_cor_df)
  
}


##Main

#Call the function to get all Unique CpG sites from the available datasets
all_CpG_sites <- Collect_CpG_sites()


#Call the function by passing all unique CpG sited to cluster the CpG sites
Clustered_CpG_sites <- Cluster_CpG_sites(all_CpG_sites)

#Call this function by passing the CpG sites cluster information to save the methylation data for each cluster and to obtain correlation data for all CpG sites
Correlation_data <- Find_correlation(Clustered_CpG_sites)

