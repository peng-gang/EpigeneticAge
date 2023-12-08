#This script is useful to download a particular dataset and to process metadata to obtain cleaned metadata
#The script has to be executed step by step by looking at the meta data and changing values as desired
#install.packages("BiocManager")
#BiocManager::install("GEOquery")
#BiocManager::install("readxl")
#BiocManager::install("httr")
#BiocManager::install("stringr")


#Load the required packages
library("GEOquery")
library(readxl)
library(httr)
library(stringr)

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
#Read the datasets information
methylation_data<- read_excel("Human_Methylation_Datasets.xlsx")


#For loop to fetch the GEO accesion number and download the raw data file and extract the metadata
for (i in 1:nrow(methylation_data)){
  
  #Create a folder with the accession number of the dataset processing in the current iteration
  dir.create(methylation_data$`Accession Number`[i])
  
  #Obtain the URL of the series matrix file
  url <- methylation_data$`ftp link`[i]
  
  #Set the destination for the matrix file to be downloaded
  dest <- paste0(getwd(),"/",methylation_data$`Accession Number`[i],"/matrix.txt.gz")
  
  #Download the dataset's matrix file
  download.file(url,dest)
  
  #Read the downloaded file using getGEO() function
  file <- paste0(methylation_data$`Accession Number`[i],"/matrix.txt.gz")
  gse <- getGEO(filename = file,GSEMatrix=TRUE)
  
  #Retrieve the methylation data and metadata and save them in the respective files
  write.csv(gse@phenoData@data,file = paste0(methylation_data$`Accession Number`[i],"/","metadata.csv"))
  write.csv(gse@assayData[["exprs"]],file = paste0(methylation_data$`Accession Number`[i],"/","methylation_data.csv"))
}

#It is not always possible to download and process all the datasets successfully using a loop
#Some times the download may fail or sometimes the getGEO() might throw an error
#So it is recommended to manually download the dataset by executing one step at once
#Code to download a particular dataset manually
#Retrieve the GEO accession of the dataset to be downloaded
GSE <- "GSE42861"

#Create a folder with the accession number of the dataset
dir.create(GSE)

#Obtain the URL of the series matrix file
url <- methylation_data$`ftp link`[methylation_data$`Accession Number` == GSE][1]

#Set the destination for the matrix file to be downloaded
dest <- paste0(getwd(),"/",GSE,"/matrix.txt.gz")

#Download the dataset's matrix file
download.file(url,dest)

#Read the downloaded file using getGEO() function
file <- paste0(GSE,"/matrix.txt.gz")
gse <- getGEO(filename = file,GSEMatrix=TRUE)

#Retrieve the methylation data and metadata and save them in the respective files
write.csv(gse@phenoData@data,file = paste0(GSE,"/","metadata.csv"))
write.csv(gse@assayData[["exprs"]],file = paste0(GSE,"/","methylation_data.csv"))

#After the data is downloaded we have to clean the sample data and save only those variables that are important
# Cleaning each dataset manually

#Retrieve the GEO accession number for which the data has to be cleaned
GSE <- paste0("GSE","42861")

#Read the metadata
data <- read.csv(paste0(GSE,'/metadata.csv'))


#Extract the sample_id and age and other variables such as tissue type and disease state
#One has to change the column names to retrieve the data as they might not be same across all datasets
cleaned_data <- data.frame(data$geo_accession,data$characteristics_ch1.14,data$type,data$tissue.ch1)

#Assign column names
colnames(cleaned_data) <- c("sample_id","age","type","tissue")

#View sample data
head(cleaned_data)

#Extract only numerical values for age as some times the age data might contain strings
for(i in 1:nrow(cleaned_data)){
  value <- cleaned_data$age[i]
  value <- str_extract(value, "\\d+")
  cleaned_data$age[i] <- floor(as.numeric(value))
}


#Extract gender data and append it to cleaned data
cleaned_data$gender <- data$characteristics_ch1.1

#Clean the gender data for uniformity, it shall only contain either male for males and female for females and not in any other format
for(i in 1:nrow(cleaned_data)){
  value <- cleaned_data$gender[i]
  value <- sub("gender: ", "", value)
  cleaned_data$gender[i] <- tolower(value)
  # if(tolower(value) == "m"){
  #   cleaned_data$gender[i] <- "male"
  # }else if (tolower(value) == "f"){
  #   cleaned_data$gender[i] <- "female"
  # }
}


#Change unknown values to NA in any column
for(i in 1:nrow(cleaned_data)){
  if(cleaned_data$gender[i] == "na"){
    cleaned_data$gender[i] <- NA
  }
}


#Extract the type of data append it to cleaned data. Ignore if already done
cleaned_data$type <- data$type

#Extract the disease_state data append it to cleaned data. Ignore if already done
cleaned_data$disease_state <- data$characteristics_ch1.4

#Clean the disease_state data
for(i in 1:nrow(cleaned_data)){
  value <- cleaned_data$disease_state[i]
  #value <- sub("disease state: ", "", value)
  value <- gsub(".*?;","",value)
  cleaned_data$disease_state[i] <- value
}


#Extract the treatment_condition data append it to cleaned data.
cleaned_data$treatment_condition <- data$characteristics_ch1.2

#Add treatment_condition to cleaned data
for(i in 1:nrow(cleaned_data)){
  value <- cleaned_data$treatment_condition[i]
  value <- sub("treatment: ", "", value)
  if(tolower(value) == "control"){
    cleaned_data$treatment_condition[i] <- "no"
  }else {
    cleaned_data$treatment_condition[i] <- "yes"
  }
}

#Extract the treatment_details data append it to cleaned data.
cleaned_data$treatment_details <- data$characteristics_ch1.2

#Add treatment details to cleaned data
for(i in 1:nrow(cleaned_data)){
  value <- cleaned_data$treatment_details[i]
  value <- sub("treatment: ", "", value)
  cleaned_data$treatment_details[i] <- tolower(value)
  # if(tolower(value) == "no"){
  #   cleaned_data$treatment_details[i] <- "control"
  # }else {
  #   cleaned_data$treatment_details[i] <- "t2d medication"
  # }
}


#Extract the tissue data append it to cleaned data.
cleaned_data$tissue <- data$characteristics_ch1.3

#Add tissue details
for(i in 1:nrow(cleaned_data)){
  value <- cleaned_data$tissue[i]
  value <- sub("tissue: ", "", value)
  cleaned_data$tissue[i] <- tolower(value)
  # if(tolower(value) == "no"){
  #   cleaned_data$treatment_details[i] <- "control"
  # }else {
  #   cleaned_data$treatment_details[i] <- "t2d medication"
  # }
}

#Remove unwanted columns if desired
cleaned_data$disease_state <- NULL


#Write the cleaned data to a csv file
write.csv(cleaned_data,file = paste0(GSE,"/","cleaned_metadata.csv"),row.names = FALSE)
