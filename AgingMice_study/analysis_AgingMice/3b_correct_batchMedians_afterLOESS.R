#load the libraries
library(readr)
library(proBatch)

# load the data
sample_annotation_AgingMice = read_csv("data_AgingMice/3_data_for_plots/sample_annotation_AgingMice.csv")
loess_fit_75_df = read_csv("data_AgingMice/3_data_for_plots/adjusted_fit_df_agingMice.csv")

batchCorrected_df_AgingMice = center_feature_batch_medians_df(loess_fit_75_df, 
                                                              sample_annotation_AgingMice, 
                                                              batch_col = 'MS_batch',  
                                                              feature_id_col = 'peptide_group_label', 
                                                              sample_id_col = 'FullRunName',
                                                              measure_col = 'Intensity', 
                                                              keep_all = T, 
                                                              no_fit_imputed = F, qual_col = NULL)

# save new data frame
write_csv(batchCorrected_df_AgingMice, 
          path = "data_AgingMice/3_data_for_plots/batchCorrected_proteome_AgingMice.csv")
