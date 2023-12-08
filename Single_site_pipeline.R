#Use this script to generate correlation plot between beta values from all datasets for a single CpG site and sample age
#Load the required libraries
library(readxl)
library(httr)
library(stringr)
library(data.table)
library(dplyr)
library(readxl)
library(ggplot2)



#The function takes CpG ID as the input and collates the beta values for this CpG site and generates correlation plot for different age groups
find_methylation_change <- function(cpg_site){

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
  
  
  #For loop for iterating through each dataset using variable GSE and combining beta values for the CpG site in cpg_site variable
  for (i in 1:length(GSE)){
    
    #Print the Dataset that is being processed in the current iteration
    print(paste0("Merging dataset: ",GSE[i]))
    
    #If this is the first iteration then we need not append anything to the cumulative data frame
    if(i == 1){
      
      #Read the methylation data from the current dataset 
      matrix <- fread(paste0("/geode2/home/u030/harsurya/Carbonate/Methylation_Datasets/Processed_Data/",GSE[i],"/methylation_data.csv"))
      
      #Read the sample data for the current dataset
      sample_data_cum <-fread(paste0("/geode2/home/u030/harsurya/Carbonate/Methylation_Datasets/Processed_Data/",GSE[i],"/cleaned_metadata.csv"))
      
      #Retrieve the CpG site names into the variable row_names
      row_names <- matrix$V1
      
      #Remove the first column from the data frame as the CpG site names were already retrieved in the previous step
      matrix<- matrix[, -1]
      
      #Check if the samples in the matrix and in the sample data are in same order or not
      if(identical(colnames(matrix),sample_data_cum$sample_id)) {print("Samples are in order")} else{print("Samples are not in order")}
      
      #The beta values matrix is read as data table - convert it to the matrix
      matrix<- as.matrix(matrix)
      
      #Assign the rownames as they are lost after conversion to matrix
      rownames(matrix) <- row_names
      
      #Transpose the matrix to make CpG sites as columns and samples as rows  - will be easy for correlation
      matrix <- t(matrix)
      
      #Convert the matrix into data table and assign the sample ids as row names
      matrix <- as.data.table(matrix)
      rownames(matrix) <- sample_data_cum$sample_id
      
      #Check if the CpG site has data in the current dataset as it is not necessary that a particular CpG site will have data in all the datasets
      #If the CpG site is present in the dataset
      if(cpg_site %in% colnames(matrix)){
        
        #Retrieve the data for the current CpG site from the current dataset and assign that to cumulative data table and assign column names
        cpg_dt_cum <- data.table(rownames(matrix),matrix[[cpg_site]])
        colnames(cpg_dt_cum) <- c("sample_id_cpg",paste0(cpg_site," Beta values"))
        
        #If the CpG site is not present in the dataset
      } else{
        
        #Generate NA values for the CpG site for the current dataset and assign the same to the cumulative data table and assign column names
        na_values <- rep(NA,nrow(matrix))
        cpg_dt_cum <- data.table(rownames(matrix),na_values)
        colnames(cpg_dt_cum) <- c("sample_id_cpg",paste0(cpg_site," Beta values"))
      }
    
    }
    #If the current iteration is not the first iteration
    else{
      
      
      #Read the methylation data from the current dataset 
      matrix <- fread(paste0("/geode2/home/u030/harsurya/Carbonate/Methylation_Datasets/Processed_Data/",GSE[i],"/methylation_data.csv"))
      
      #Read the sample data for the current dataset
      sample_data <-fread(paste0("/geode2/home/u030/harsurya/Carbonate/Methylation_Datasets/Processed_Data/",GSE[i],"/cleaned_metadata.csv"))
      
      #Retrieve the CpG site names into the variable row_names
      row_names <- matrix$V1
      
      #Remove the first column from the data frame as the CpG site names were already retrieved in the previous step
      matrix<- matrix[, -1]  
      
      #Check if the samples in the matrix and in the sample data are in same order or not
      if(identical(colnames(matrix),sample_data$sample_id)) {print("Samples are in order")} else{print("Samples are not in order")}
      
      #The beta values matrix is read as data table - convert it to the matrix
      matrix <- as.matrix(matrix)
      
      #Assign the rownames as they are lost after conversion to matrix
      rownames(matrix) <- row_names
      
      #Transpose the matrix to make CpG sites as columns and samples as rows  - will be easy for correlation
      matrix<- t(matrix)
      
      #Convert the matrix into data table and assign the sample ids as row names
      matrix <- as.data.table(matrix)
      rownames(matrix) <- sample_data$sample_id
    
      
      #Check if the CpG site has data in the current dataset as it is not necessary that a particular CpG site will have data in all the datasets
      #If the CpG site is present in the dataset
      if(cpg_site %in% colnames(matrix)){
        
        #Retrieve the data for the current CpG site from the current dataset and assign column names
        cpg_dt <- data.table(rownames(matrix),matrix[[cpg_site]])
        colnames(cpg_dt) <- c("sample_id_cpg",paste0(cpg_site," Beta values"))
        
        #If the CpG site is not present in the dataset
      } else{
        
        #Generate NA values for the CpG site for the current dataset  and assign column names
        na_values <- rep(NA,nrow(matrix))
        cpg_dt <- data.table(rownames(matrix),na_values)
        colnames(cpg_dt) <- c("sample_id_cpg",paste0(cpg_site," Beta values"))
        
        
      }
    
      #Append the data obtained from the current dataset to the cumulative data table
      cpg_dt_cum <- rbind(cpg_dt_cum,cpg_dt)
      sample_data_cum <- bind_rows(sample_data_cum,sample_data)
    
    }

  }

    #Column bind the beta values for the given CpG site and the corresponding sample data  
    cpg_plot_dt <- cbind(sample_data_cum,cpg_dt_cum)

    #If the total number of beta values are not more than 80% of the total number of samples then we ignore thr CpG site as it does not have sufficient values
    #Calculate the total percent of beta values for the given CpG site
    na_values_count <- sum(is.na(cpg_plot_dt[[paste0(cpg_site," Beta values")]]))
    total_per <- ((nrow(cpg_plot_dt) - na_values_count)/nrow(cpg_plot_dt)) * 100

  #If the total percent is greater than 80%
  if(total_per > 80){
  
    #Remove the rows where sample age is NA
    cpg_plot_dt <- cpg_plot_dt[complete.cases(cpg_plot_dt$age), ]
    
    #Remove the rows which have Beta value as NA
    cpg_plot_dt <- cpg_plot_dt[complete.cases(cpg_plot_dt[[paste0(cpg_site," Beta values")]]), ]
  
    #Retain only those rows with sample age less than or equal to 100 to avoid samples that have age in days
    cpg_plot_dt <- cpg_plot_dt[cpg_plot_dt$age <= 100]
  
    #First generate a plot that shows correlation between the beta values and age in all age groups
    plot_all <- ggplot(cpg_plot_dt,aes(age,get(paste0(cpg_site," Beta values")))) +
      geom_smooth() +
      theme_light() +
      labs(title = paste0("Change in methylation at cpg site ",cpg_site," in all age groups"), 
          x = "Age", 
          y = paste0(cpg_site," Beta values") ) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
    #Create an empty plot list to save plots of different age groups
    plot_list <- list()
    
    
    #Use variable k to iterate through plot list
    k=0
    
    #Use variable i to iterate through different age ranges
    i = 1
    
    #While loop to generate plots for each age group
    while (i <= 100){
      
      #Retrieve the Beta values and sample data of the given CpG site into a temporary data table
      cpg_plot_dt_temp <- cpg_plot_dt
      
      #Filter the rows where age is between i and i+9 which will be 0 to 10, 11 to 20, 21 to 30, etc, ....upto 91 to 100 for each iteration
      cpg_plot_dt_temp <- cpg_plot_dt_temp[cpg_plot_dt_temp$age > (i-1)]
      cpg_plot_dt_temp <- cpg_plot_dt_temp[cpg_plot_dt_temp$age <= (i+9)]
    
      #Generate the plot for the current age group
      plot <-  ggplot(cpg_plot_dt_temp,aes(age,get(paste0(cpg_site," Beta values")))) +
        geom_smooth() +
        theme_light() +
        labs(title = paste0("Change in methylation at cpg site ",cpg_site," in age group ",i-1," to ",i+9), 
            x = "Age", 
            y = paste0(cpg_site," Beta values") ) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
      #Increment the k and in the place of it in the plot list save the current plot
      k=k+1
      plot_list[[k]] <- plot
    
      #Increment i by 10 to get to the next age group
      i <- i+10
    }
  
    #Append the plot for age groups to the plot list
    plot_list[[k+1]] <- plot_all
  
    #Arrange the plots in the plot list in a pattern to print to pdf
    pl <- ggpubr::ggarrange(plotlist = plot_list,ncol = 2,nrow = 6)
    pdf(paste0(cpg_site,"_methylation_change.pdf"),height = 30, width = 15)
    print(pl)
  
    dev.off()
  
  #If CpG values number is less than the 80% of the total samples display the following message
  }else{
    print("Cpg values count is less than 80% of total samples!")
  }

}





#Main
#Call the function by passing a CpG ID to get the correlation plots for that site
find_methylation_change("cg00415665")



