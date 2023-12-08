The scripts in the project are useful to download methylation profiling data and to see the trends in methylation change with age.

The following scripts are useful to download the data:
1. Data_Cleaning.R - Has the code to download and process a particular dataset using getGEO() function.
2. Raw_File_Extract.R - Has the code to download the raw files and process them to obtain the methylation data. This is useful when methylation data is not in the series matrix file.


The next step is to collate the data from all the datasets as high number of samples will give us better trends. The following scripts will help to achieve the same:
1. Merge_data_range.R - Has the code to merge methylation data from all the available datasets with in a desired age range.
2. Sample_data_merge.R - Has the code to merge sample data from all the datasets into one file.
3. Cluster_pipeline.R - Merging data from all datasets is huge task and hence it is better to cluster the CpG sites into different clusters with desired size and combine methylation data of CpG sites present in a particular cluster. This will help to process some CpG sites at once and will reduce the computing time and resources. The code in this script will help to achieve this objective.


Lastly we want to analyse the trends in methylation change. The following scripts will help to achieve the same:
1. Get_Age.R - Has the code to tell the minimum age and maximum age of the samples in a particular dataset. This will help to decide if that dataset is to be included for analysis or not.
2. Single_Dataset.R - Has the code to generate correlation plot between beta values of desired CpG site and sample ages from a particular dataset.
3. Single_site_pipeline.R - Has the code to combine methylation data from all the datasets for a given CpG site and to generate correlation plots for all and different age groups.
4. KO_plots.R - Cluster_pipeline.R will give the clustered CpG sites data. Using this data the KO_plots.R script will generate the correlation plots for the desired CpG sites by excluding one dataset at a time. It also has code to exclude many undesired datasets at once from the correlation plot. While the former procedure will help us to know the effect each dataset is having on the correlation plot, the latter one will help to understand the combined effect of certain datasets on the correlation plot.
