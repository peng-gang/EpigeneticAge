#Use this script to read the clustered CpG sites data and sample data and to obtain correlation plots for certain CpG sites by selectively excluding certain datasets

#Load the required libraries
library(factoextra)
library(readxl)
library(httr)
library(stringr)
library(data.table)
library(dplyr)
library(readxl)
library(ggplot2)
library(stats)

#Read the desired clustered CpG sites file - These files are generated after the CpG sites are clustered using the script Cluster_pipeline.R
beta_values <- read.csv("cluster_1_beta_values.csv")

#Retrieve CpG IDs from the first column of the data frame and assign them as  rownames for the data frame and remove the first column.
rownames(beta_values) <- beta_values$X
beta_values <- beta_values[,-1]

#Read the desired sample data file - These files are generated after the CpG sites are clustered using the script Cluster_pipeline.R
sample_data <- read.csv("cluster_1_sample_data.csv")

#Retrieve sample IDs from the first column of the data frame and assign them as  rownames for the data frame and remove the first column.
rownames(sample_data) <- sample_data$sample_id
sample_data <- sample_data[,-1]


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

#The inner for loop below generates the correlation plots by excluding one data set for each iteration
#The outer for loop also makes a function call to KO_plots that takes in the vector of datasets that needs to be excluded and generates a single correlation plot by excluding all those datasets at once
#Select the desired condition to include or exclude certain datasets. Here a dataset with GSE accession number GSE124076 is desired. Hence, the rest of the accession numbers are passed on to KO_plots function
#Whatever GEO accession is stored in the GSE2, it is not used for correlation plot in the KO_plots function
condition <- GSE!="GSE124076"
GSE2 <- GSE[condition]



#The outer loop is executed for the first 5 CpG sites - change the condition as desired
for(j in 1:length(head(colnames(beta_values)))){
  
  #Retrieve the jth CpG site from the first 5 CpG sites
  CpG_site <- head(colnames(beta_values))[j]
  
  #For loop to generate correlation plots by excluding one dataset at a time
  #This will help to know the affect that this dataset is causing to the correlation plot
  for(i in 1:length(GSE)){
    
    #Generate a boolean column that has value TRUE for all the samples except the samples belonging to the current dataset in GSE
    idx <- sample_data$GEO_Accesion != GSE[i]
    
    #Retrieve the desired beta values and sample data using the boolean column
    beta_values_tmp <- beta_values[idx,]
    sample_data_tmp <- sample_data[idx,]
    
    #Combine the beta values and sample age into a data frame to generate plots
    plot_df <- data.frame(beta_values_tmp[[CpG_site]],sample_data_tmp$age)
    
    #Assign the column names
    colnames(plot_df) <- c(CpG_site,"age")
    
    #Generate and save the correlation plot
    gg <- ggplot(plot_df,aes(age,get(CpG_site))) +
      geom_smooth(method = "gam") +
      geom_point() +
      theme_light() +
      labs(title = paste0(CpG_site," : ",GSE[i]," Knock out"), 
           x = "Age", 
           y = paste0(CpG_site," Beta values") ) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    #Save the correlation plot in a desired location
    ggsave(paste0("Cluster_results/Knockouts/",CpG_site,"/",GSE[i],"_KO.png"),plot =gg)
  }
  
  #Call the function KO_plots to remove all the desired datasets in GSE2 from the plot at once
  KO_plots(GSE2,CpG_site)
}

#function to generate a correlation plot by ignoring all the desired datasets at once
#It takes the vector of GEO accession numbers GSE2 that are to be ignored and the CpG site
KO_plots <- function(KO,CpG_site){
  
  #Generate a string to include the GEO accession numbers of all excluding datasets to use as name for the plot
  KO_str <- KO[1]
  for(i in 2:length(KO)){
    KO_str <- paste0(KO_str,"_",KO[i])
  }
  
  #Create a boolean variable to filter the samples that do not belong to the datasets in KO
  idx <- !(sample_data$GEO_Accesion %in% KO)
  
  #Use the boolean variable to extract the desired beta values and sample data
  beta_values_tmp <- beta_values[idx,]
  sample_data_tmp <- sample_data[idx,]
  
  #Combine the beta values and sample age into a data frame to generate plots
  plot_df <- data.frame(beta_values_tmp[[CpG_site]],sample_data_tmp$age)
  
  #Assign the column names
  colnames(plot_df) <- c(CpG_site,"age")
  
  #Generate and save the correlation plot
  gg <- ggplot(plot_df,aes(age,get(CpG_site))) +
    geom_smooth(method = "gam") +
    geom_point() +
    theme_light() +
    labs(title = paste0(CpG_site," : ",KO_str," Knock out"), 
         x = "Age", 
         y = paste0(CpG_site," Beta values") ) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  #Save the correlation plot in a desired location
  ggsave(paste0("Cluster_results/Knockouts/",CpG_site,"/All_KOs/",KO_str,"_KO.png"),plot =gg)
}



