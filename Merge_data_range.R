#Use this script to merge the data from all the datasets from all the samples between a specific age range
#Can be used only if the number of datasets avaiable is smaller
#Too many datasets will consume the RAM memory and the execution is halted
#Load required libraries
library(readxl)
library(httr)
library(stringr)
library(data.table)
library(dplyr)
library(readxl)


#Function to filter the given data and return that data that have sample age in between the given range
#Takes input of data and the column in which the data has to be filtered and the min and max age
between <- function(data,column,min,max){
  #Retrieve all the data with corresponding sample age greater than min age
  data <- data[data[[column]] >= min]
  
  #Filter out the data which have corresponding sample age greater than max age
  data <- data[data[[column]] <= max]
  
  #return the data
  return(data)
}


#Function to merge all the data from all the datasets for the samples with in the given range
Merge_samples_in_range <- function(min,max){
  
  
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
  datasets <- read_excel("/geode2/home/u030/harsurya/Carbonate/Methylation Datasets/DNAMethylation/Human_Methylation_Datasets.xlsx",sheet = "All datasets")
  #Get the accession numbers of the datasets which are marked as useful into the variable GSE
  GSE <- datasets$`Accession Number`[datasets$Remarks == "Useful"]
  
  #Create empty datasets to store the data
  cumulative_samples <- data.table()
  cumulative_matrix <- data.table()
  
  print("Initializing merging of datasets")
  
  #For loop to iterate through CpG sites of each dataset and gather methylation data for CpG sites from all datasets and then save it to the local computer  
  for(i in 1:length(GSE)){
    
    #Print the current dataset being used to collect methylation data
    print(paste0("Merging dataset: ",GSE[i]))
    
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
    matrix <- fread(paste0("/geode2/home/u030/harsurya/Carbonate/Methylation Datasets/Processed_Data/",GSE[i],"/methylation_data.csv"))
    
    #Read the sample data for the current dataset using cleaned_metadata.csv
    sample_data <-fread(paste0("/geode2/home/u030/harsurya/Carbonate/Methylation Datasets/Processed_Data/",GSE[i],"/cleaned_metadata.csv"))
    
    #Save the row names of the matrix i.e., the names of the CpG sites
    row_names <- matrix$V1
    
    #The first column contains the CpG site Ids, it can be removed as they are saved in the row_names
    matrix<- matrix[, -1]
    
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
    
    #Filter the samples that are between the given age range using between function in this script and retain only them
    filtered_samples <- between(sample_data,"age",min,max)
    matrix <- matrix[rownames(matrix) %in% filtered_samples$sample_id]
    
    #If the current iteration is the first iteration just add the data to the cumulative data tables if not append the current data to them.
    if(i==1){
      
      #Assign the current dataset information to the cumulative table
      cumulative_samples <- filtered_samples
      cumulative_matrix <- matrix
      rownames(cumulative_matrix) <- cumulative_samples$sample_id
      
    }else{
      
      #Append the current dataset information to the cumulative table
      #merge function takes time to execute as it has map each CpG site from each dataset to the CpG sites in the cumulative tables and then append the data
      cumulative_samples <- bind_rows(cumulative_samples,filtered_samples)
      cumulative_matrix <- merge(cumulative_matrix,matrix, all = TRUE, sort = FALSE)
      rownames(cumulative_matrix) <- cumulative_samples$sample_id
      
    }
    
    
    
  }
  
  #Return the data tables as a list
  print("Finished merging")
  
  return_tables <- list(table1 = cumulative_matrix, table2 = cumulative_samples)
  return(return_tables)
  
}


#Main

#Read the age values from the execution command if the script is executed from terminal
args <- commandArgs(trailingOnly = TRUE)

#Assign the desired age values
if(length(args) == 0){
  min = 1
  max = 70
}else{
  min = as.numeric(args[1])
  max = as.numeric(args[2])
}

#Call the function Merge_samples_in_range with the desired age range 
print(paste0("Merging samples between ages ",min," and ",max))
dataset_information <- Merge_samples_in_range(min,max)

#Collect the methylation and sample information from the returned list
methylation_data <- dataset_information$table1
sample_data <- dataset_information$table2


