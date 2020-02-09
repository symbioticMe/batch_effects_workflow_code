#!/usr/bin/env Rscript --vanilla

#load the libraries
library(readr)
library(proBatch)

# load the data
sample_annotation_AgingMice = read_csv("data_AgingMice/3_data_for_plots/sample_annotation_AgingMice.csv")
normalized_df_AgingMice = read_csv(file = "data_AgingMice/3_data_for_plots/normalized_proteome_AgingMice.csv")

# Batch effect correction 
loess_fit_75 <- adjust_batch_trend_df(normalized_df_AgingMice, sample_annotation_AgingMice,
                                      span = 0.75, no_fit_imputed = FALSE)

write_csv(loess_fit_75, 
          path = "data_AgingMice/3_data_for_plots/adjusted_fit_df_agingMice.csv")


loess_fit_25 <- adjust_batch_trend_df(normalized_df_AgingMice, sample_annotation_AgingMice,
                                      span = 0.25, no_fit_imputed = FALSE)
write_csv(loess_fit_25, 
          path = "data_AgingMice/3_data_for_plots/adjusted_fit_25_agingMice_no_requants.csv")



loess_fit_150 <- adjust_batch_trend_df(normalized_df_AgingMice, sample_annotation_AgingMice,
                                      span = 1.5, no_fit_imputed = FALSE)

write_csv(loess_fit_150, 
          path = "data_AgingMice/3_data_for_plots/adjusted_fit_150_agingMice_no_requants.csv")
