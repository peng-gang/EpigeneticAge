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

#Read the datasets information
methylation_data<- read_excel("Human_Methylation_Datasets.xlsx")

methylation_data$`Accession Number`[27]

#for loop to fetch the GEO accesion number and download the raw data file and extract the metadata
for (i in 27:nrow(methylation_data)){
  
  dir.create(methylation_data$`Accession Number`[i])
  url <- methylation_data$`ftp link`[i]
  dest <- paste0(getwd(),"/",methylation_data$`Accession Number`[i],"/matrix.txt.gz")
  download.file(url,dest)
  file <- paste0(methylation_data$`Accession Number`[i],"/matrix.txt.gz")
  gse <- getGEO(filename = file,GSEMatrix=TRUE)
  write.csv(gse@phenoData@data,file = paste0(methylation_data$`Accession Number`[i],"/","metadata.csv"))
  write.csv(gse@assayData[["exprs"]],file = paste0(methylation_data$`Accession Number`[i],"/","methylation_data.csv"))
}



#Removing the datasets that have problem in downloading
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE42861",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE43414",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE51032",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE198904",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE210843",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE213478",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE224363",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE224364",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE224365",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE224624",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE227815",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE104210",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE55763",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE68838",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE77696",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE104210",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE112611",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE46376",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE152026",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE154566",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE157131",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE67202",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE72365",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE72368",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE174422",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE179847",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE179849",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE180474",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE185445",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE185920",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE190931",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE201754",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE203332",]


methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE167202",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE127368",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE172368",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE172365",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE146376",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE115278",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE124413",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE128235",]
methylation_data <- methylation_data[!methylation_data$`Accession Number` == "GSE131989",]

write.csv(methylation_data,file = paste0("Filtered_Methylation_Datasets.csv"))


#Datasets with some problematic or wrong data
#GSE185445 - mother age and father age - ignored
#GSE213478 - Sex in numericals and age in groups - lower age saved - deidentified data

#Code to download a particular dataset manually
GSE <- "GSE42861"
dir.create(GSE)
url <- methylation_data$`ftp link`[methylation_data$`Accession Number` == GSE][1]
dest <- paste0(getwd(),"/",GSE,"/matrix.txt.gz")
download.file(url,dest)
file <- paste0(GSE,"/matrix.txt.gz")
gse <- getGEO(filename = file,GSEMatrix=TRUE)
write.csv(gse@phenoData@data,file = paste0(GSE,"/","metadata.csv"))
write.csv(gse@assayData[["exprs"]],file = paste0(GSE,"/","methylation_data.csv"))

# Cleaning each dataset manually
#Saving the row number of a particular dataset to read the information for cleaning
# k <- 40
GSE <- paste0("GSE","42861")
#Read the data
# data <- read.csv(paste0(methylation_data$`Accession Number`[k],'/metadata.csv'))
data <- read.csv(paste0(GSE,'/metadata.csv'))
cleaned_data <- read.csv(paste0(GSE,'/cleaned_metadata.csv'))

#Extract the sample_id and age
cleaned_data <- data.frame(data$geo_accession,data$characteristics_ch1.14,data$type,data$tissue.ch1)

#Assign column names
colnames(cleaned_data) <- c("sample_id","age","type","tissue")

#View sample data
head(cleaned_data)

#Extract only numerical values for age
for(i in 1:nrow(cleaned_data)){
  value <- cleaned_data$age[i]
  value <- str_extract(value, "\\d+")
  cleaned_data$age[i] <- floor(as.numeric(value))
}

#Extract gender data
cleaned_data$gender <- data$characteristics_ch1.1

#Clean the gender data for uniformity
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

#Extract the type of data
cleaned_data$type <- data$type

#Extract the disease_state data
cleaned_data$disease_state <- data$characteristics_ch1.4

#Clean the disease_state data
for(i in 1:nrow(cleaned_data)){
  value <- cleaned_data$disease_state[i]
  #value <- sub("disease state: ", "", value)
  value <- gsub(".*?;","",value)
  cleaned_data$disease_state[i] <- value
}

# #Write the cleaned data to a csv file
# write.csv(cleaned_data,file = paste0(methylation_data$`Accession Number`[k],"/","cleaned_metadata.csv"),row.names = FALSE)
# 
# 
# #Add additional variables to the cleaned data
# #Read the data
# cleaned_data <- read.csv(paste0(methylation_data$`Accession Number`[k],'/cleaned_metadata.csv'))

cleaned_data$treatment_condition <- data$characteristics_ch1.2

#Add treatment details
for(i in 1:nrow(cleaned_data)){
  value <- cleaned_data$treatment_condition[i]
  value <- sub("treatment: ", "", value)
  if(tolower(value) == "control"){
    cleaned_data$treatment_condition[i] <- "no"
  }else {
    cleaned_data$treatment_condition[i] <- "yes"
  }
}


cleaned_data$treatment_details <- data$characteristics_ch1.2

#Add treatment details
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

#Remove unwanted columns
cleaned_data$disease_state <- NULL



cleaned_data$tissue <- data$characteristics_ch1.3

#Add treatment details
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



#Write the cleaned data to a csv file
# write.csv(cleaned_data,file = paste0(methylation_data$`Accession Number`[k],"/","cleaned_metadata.csv"),row.names = FALSE)
#Write the cleaned data to a csv file
write.csv(cleaned_data,file = paste0(GSE,"/","cleaned_metadata.csv"),row.names = FALSE)

g <- getGEOSuppFiles("GSE185445")

print(methylation_data,n=100)
