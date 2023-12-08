#Use this script to generate plots for correlation between beta values of CpG sites and age from a single dataset
#Load required libraries
library(readxl)
library(httr)
library(stringr)
library(data.table)
library(dplyr)
library(readxl)
library(ggplot2)



#The function to generate the correlation plots for beta values of CpG sites and age
#The function generates plots only for first 10 CpG sites in the dataset by default - Adjust the range in for loop to get the plots for desired CpG sites
#The function takes the accession number of the dataset
find_methylation_change <- function(Dataset){
  
    #Read the methylation data of the dataset
    matrix <- fread(paste0("/geode2/home/u030/harsurya/Carbonate/Methylation_Datasets/Processed_Data/",Dataset,"/methylation_data.csv"))
    
    #Read sample data of the dataset
    sample_data <-fread(paste0("/geode2/home/u030/harsurya/Carbonate/Methylation_Datasets/Processed_Data/",Dataset,"/cleaned_metadata.csv"))
    
    #Extract row names of the matrix which are the CpG site IDs
    row_names <- matrix$V1
    
    #Remove column one in the matrix as it only contains the CpG IDs and they are retrieved in the previous step
    matrix<- matrix[, -1]
    
    #Check if the samples in the matrix and in the sample data are in same order or not
    if(identical(colnames(matrix),sample_data$sample_id)) {print("Samples are in order")} else{print("Samples are not in order")}
    
    #The beta values matrix is read as data table - convert it to the matrix
    matrix<- as.matrix(matrix)
    
    #Assign the rownames as they are lost after conversion to matrix
    rownames(matrix) <- row_names
    
    
    #Transpose the matrix to make CpG sites as columns and samples as rows  - will be easy for correlationmatrix <- t(matrix)
    matrix <- as.data.table(matrix)
    
    #Assign the sample ids as row names
    rownames(matrix) <- sample_data$sample_id

    #For loop to iterate through 1st 10 columns to generate correlation plot for first 10 CpG sites
    #Change the range in for loop to get plots for desired CpG sites
    for(i in 1:10){
      
      #Retrieve the current CpG site
      site <- colnames(matrix)[i]
      
      #Prepare a temporary dataframe containing the beta values of the site and the sample age
      plot_df <- data.frame(matrix[[site]],sample_data$age)
      
      #Assign column names to the plot dataframe
      colnames(plot_df) <- c(site,"age")
    

    #Generate the plot and save in gg
      gg <- ggplot(plot_df,aes(age,get(site))) +
        geom_point()+
        geom_smooth(method = "gam") +
        theme_light() +
        labs(title = paste0("Change in methylation at cpg site ",site), 
            x = "Age", 
            y = paste0(site," Beta values") ) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
      #Save the plot to the desired location
      ggsave(paste0("Cluster_results/Single_Dataset/",Dataset,"/",site,".png"),plot =gg)
    
  }
  
    
    
  
}


#Call the function to generate correlation plots for the CpG sites by passing the GEO accession of the dataset
find_methylation_change("GSE53740")

