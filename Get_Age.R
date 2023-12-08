#Use this script to identify the minimum and maximum ages for samples of a particular dataset
#Load the required libraries
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
datasets <- read_excel("/geode2/home/u030/harsurya/Carbonate/Methylation_Datasets/DNAMethylation/Human_Methylation_Datasets.xlsx",sheet = "All datasets")

#Assign the accession number of the dataset for which we want to find the min and max ages
Accession_Number <- "GSE179849"

#Create a folder with the name of Accession_Number to save all the dataset files in there
dir.create(Accession_Number)

#Retrieve the url of the series matrix file for the dataset using Accession_Number 
url <- datasets$`ftp link`[datasets$`Accession Number`==Accession_Number]

#Set a destination to download the series matrix file
dest <- paste0(getwd(),"/",Accession_Number,"/matrix.txt.gz")

#Download the file
download.file(url,dest)

#Read the file using getGEO()
file <- paste0(Accession_Number,"/matrix.txt.gz")
gse <- getGEO(filename = file,GSEMatrix=TRUE)

#Get the age data using appropriate column name. It varies from dataset to dataset
age <- gse@phenoData@data$characteristics_ch1.8


#For loop to convert all the age values into numericals by removing strings if any
for(i in 1:length(age)){
  value <- age[i]
  value <- str_extract(value, "\\d+")
  age[i] <- as.numeric(floor(as.numeric(value)))
}

#Convert the type of the vector to numerics
age <- as.numeric(age)

#Sort the age values
age <- sort(age)

#Print head and tail of the ages
head(age,n = 10)
tail(age,n = 10)

#Print minimum and maximum ages of the samples and manually update this information for the dataset in the Human_Methylation_Datasets.xlsx
min(age,na.rm = TRUE)
max(age,na.rm = TRUE)

