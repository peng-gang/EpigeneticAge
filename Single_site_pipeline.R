library(readxl)
library(httr)
library(stringr)
library(data.table)
library(dplyr)
library(readxl)
library(ggplot2)




find_methylation_change <- function(cpg_site){

  datasets <- read_excel("/geode2/home/u030/harsurya/Carbonate/Methylation Datasets/DNAMethylation/Human_Methylation_Datasets.xlsx",sheet = "All datasets")

  GSE <- datasets$`Accession Number`[datasets$Remarks == "Useful"]

  # GSE <- c("GSE101764","GSE106648","GSE114134","GSE124076","GSE124366")

  #GSE114135
  #GSE146377



  for (i in 1:length(GSE)){
    
    print(paste0("Merging dataset: ",GSE[i]))
    if(i == 1 ){
      matrix <- fread(paste0("/geode2/home/u030/harsurya/Carbonate/Methylation Datasets/Processed_Data/",GSE[i],"/methylation_data.csv"))
      sample_data_cum <-fread(paste0("/geode2/home/u030/harsurya/Carbonate/Methylation Datasets/Processed_Data/",GSE[i],"/cleaned_metadata.csv"))
      row_names <- matrix$V1
      matrix<- matrix[, -1]  # Remove the first column from the data frame
      if(identical(colnames(matrix),sample_data_cum$sample_id)) {print("Samples are in order")} else{print("Samples are not in order")}
      matrix<- as.matrix(matrix)
      rownames(matrix) <- row_names
      matrix <- t(matrix)
      matrix <- as.data.table(matrix)
      rownames(matrix) <- sample_data_cum$sample_id
      if(cpg_site %in% colnames(matrix)){
        cpg_dt_cum <- data.table(rownames(matrix),matrix[[cpg_site]])
        colnames(cpg_dt_cum) <- c("sample_id_cpg",paste0(cpg_site," Beta values"))
      } else{
        na_values <- rep(NA,nrow(matrix))
        cpg_dt_cum <- data.table(rownames(matrix),na_values)
        colnames(cpg_dt_cum) <- c("sample_id_cpg",paste0(cpg_site," Beta values"))
      }
    
    }
  
    else{
    
      matrix <- fread(paste0("/geode2/home/u030/harsurya/Carbonate/Methylation Datasets/Processed_Data/",GSE[i],"/methylation_data.csv"))
      sample_data <-fread(paste0("/geode2/home/u030/harsurya/Carbonate/Methylation Datasets/Processed_Data/",GSE[i],"/cleaned_metadata.csv"))
      row_names <- matrix$V1
      matrix<- matrix[, -1]  # Remove the first column from the data frame
      if(identical(colnames(matrix),sample_data$sample_id)) {print("Samples are in order")} else{print("Samples are not in order")}
      matrix <- as.matrix(matrix)
      rownames(matrix) <- row_names
      matrix<- t(matrix)
      matrix <- as.data.table(matrix)
      rownames(matrix) <- sample_data$sample_id
    
      if(cpg_site %in% colnames(matrix)){
        cpg_dt <- data.table(rownames(matrix),matrix[[cpg_site]])
        colnames(cpg_dt) <- c("sample_id_cpg",paste0(cpg_site," Beta values"))
      } else{
        na_values <- rep(NA,nrow(matrix))
        cpg_dt <- data.table(rownames(matrix),na_values)
        colnames(cpg_dt) <- c("sample_id_cpg",paste0(cpg_site," Beta values"))
      }
    
      cpg_dt_cum <- rbind(cpg_dt_cum,cpg_dt)
      sample_data_cum <- bind_rows(sample_data_cum,sample_data)
    
    }

  }

    cpg_plot_dt <- cbind(sample_data_cum,cpg_dt_cum)

    na_values_count <- sum(is.na(cpg_plot_dt[[paste0(cpg_site," Beta values")]]))

    total_per <- ((nrow(cpg_plot_dt) - na_values_count)/nrow(cpg_plot_dt)) * 100

  if(total_per > 80){
  
    cpg_plot_dt <- cpg_plot_dt[complete.cases(cpg_plot_dt$age), ]
    cpg_plot_dt <- cpg_plot_dt[complete.cases(cpg_plot_dt[[paste0(cpg_site," Beta values")]]), ]
  
    cpg_plot_dt <- cpg_plot_dt[cpg_plot_dt$age <= 100]
  
  
    plot_all <- ggplot(cpg_plot_dt,aes(age,get(paste0(cpg_site," Beta values")))) +
      geom_smooth() +
      theme_light() +
      labs(title = paste0("Change in methylation at cpg site ",cpg_site," in all age groups"), 
          x = "Age", 
          y = paste0(cpg_site," Beta values") ) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
    plot_list <- list()
  
    k=0
    i = 1
    while (i <= 100){
      cpg_plot_dt_temp <- cpg_plot_dt
      cpg_plot_dt_temp <- cpg_plot_dt_temp[cpg_plot_dt_temp$age > (i-1)]
      cpg_plot_dt_temp <- cpg_plot_dt_temp[cpg_plot_dt_temp$age <= (i+9)]
    
      plot <-  ggplot(cpg_plot_dt_temp,aes(age,get(paste0(cpg_site," Beta values")))) +
        geom_smooth() +
        theme_light() +
        labs(title = paste0("Change in methylation at cpg site ",cpg_site," in age group ",i-1," to ",i+9), 
            x = "Age", 
            y = paste0(cpg_site," Beta values") ) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
      k=k+1
      plot_list[[k]] <- plot
    
      i <- i+10
    }
  
    plot_list[[k+1]] <- plot_all
  
    pl <- ggpubr::ggarrange(plotlist = plot_list,ncol = 2,nrow = 6)
    pdf(paste0(cpg_site,"_methylation_change.pdf"),height = 30, width = 15)
    print(pl)
  
    dev.off()
  
  }else{
    print("Cpg values count is less than 80% of total samples!")
  }

}

find_methylation_change("cg00000734")



