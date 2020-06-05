library(readr)
library(proBatch)


#load the data
raw_proteome_PanCancer = read_csv('data_PanCancer/2_interim_data/raw_transitome_PanCancer.csv')

raw_proteome_PanCancer = log_transform_df(raw_proteome_PanCancer,
                                          log_base = 2, offset = 1)
normalized_transitome_PanCancer = quantile_normalize_df(raw_proteome_PanCancer, 
                                               feature_id_col = 'transition', 
                                               sample_id_col = 'sample_name', 
                                               no_fit_imputed = F)

write_csv(normalized_transitome_PanCancer, 'data_PanCancer/2_interim_data/normalized_transitome_PanCancer.csv')
