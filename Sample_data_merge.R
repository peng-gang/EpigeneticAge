#Use this script to merge all samples from all the available datasets 
#Load required libraries
library(readxl)
library(httr)
library(stringr)
library(data.table)
library(dplyr)
library(readxl)

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
print("Collecting dataset information")
datasets <- read_excel("/geode2/home/u030/harsurya/Carbonate/Methylation_Datasets/DNAMethylation/Human_Methylation_Datasets.xlsx",sheet = "All datasets")

#Get the accession numbers of the datasets which are marked as useful into the variable GSE
GSE <- datasets$`Accession Number`[datasets$Remarks == "Useful"]


#For loop to iterate through each dataset and collect and collate sample information
for(i in 1:length(GSE)){
  
  #Print the accession number of the data set being processed in the current iteration
  print(paste0("Collecting Samples from dataset: ",GSE[i]))
  
  #For the first iteration there is no need to append the data. Hence, read the data directly to the cumulative sample data table
  if(i==1){
    
    cum_sample_data <- fread(paste0("/geode2/home/u030/harsurya/Carbonate/Methylation_Datasets/Processed_Data/",GSE[i],"/cleaned_metadata.csv"))
    
  }#For all other iterations read and append the sample information to the cumulative sample data table
  else{
    
    sample_data <- fread(paste0("/geode2/home/u030/harsurya/Carbonate/Methylation_Datasets/Processed_Data/",GSE[i],"/cleaned_metadata.csv"))
    cum_sample_data <- bind_rows(cum_sample_data,sample_data)
    
    
    
    
  }
  
  
}


#Save the cumulative sample data table to local computer
write.csv(cum_sample_data, file = "/geode2/home/u030/harsurya/Carbonate/Methylation_Datasets/cum_sample_data.csv",row.names = FALSE)
