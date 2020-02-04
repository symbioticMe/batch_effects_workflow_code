library(readr)
library(proBatch)

normalized_transitome_PanCancer = read_csv('data_PanCancer/2_interim_data/normalized_transitome_PanCancer.csv')
sample_annotation_PanCancer = read_csv('data_PanCancer/3_data_for_plots/sample_annotation_PanCancer.csv')

corrected_transitome_proBatch = center_feature_batch_means_df(normalized_transitome_PanCancer, 
                                                              sample_annotation =  sample_annotation_PanCancer, 
                                                              feature_id_col = 'transition', 
                                                              sample_id_col = 'sample_name', 
                                                              no_fit_imputed = F, 
                                                              batch_col = 'Batch')

write_csv(corrected_transitome_proBatch, 
          'data_PanCancer/2_interim_data/batchCorrected_transitome_PanCancer.csv')