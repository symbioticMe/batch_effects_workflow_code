#!/usr/bin/env Rscript --vanilla


#load the libraries
library(readr)
library(proBatch)

# load the data
proteome_df_AgingMice = read_csv("data_AgingMice/2_interim_data/raw_proteome_AgingMice.csv")

# quantile normalize of log2 transformed data matrix 
normalized_df_AgingMice = normalize_data_df(proteome_df_AgingMice, 
                                            normalize_func = "quantile", 
                                            log_base = NULL,
                                            keep_all = 'default')

# save new data frame
write_csv(normalized_df_AgingMice, 
          path = "data_AgingMice/3_data_for_plots/normalized_proteome_AgingMice.csv")

#normalize the proteome with requant values
proteome_df_withRequants = read_csv('data_AgingMice/2_interim_data/raw_proteome_AgingMice_with_requants.csv')
normalized_df_withRequants = normalize_data_df(proteome_df_withRequants, 
                                            normalize_func = "quantile", 
                                            log_base = NULL,
                                            keep_all = 'default')

# save new data frame
write_csv(normalized_df_withRequants, 
          path = "data_AgingMice/3_data_for_plots/normalized_proteome_AgingMice_with_requants.csv")

