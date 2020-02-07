#!/usr/bin/env Rscript --vanilla


#load the libraries
library(readr)
library(ggplot2)
library(dplyr)
library(proBatch)

# load data
sample_annotation_AgingMice = read_csv("data_AgingMice/3_data_for_plots/sample_annotation_AgingMice.csv")
normalized_df_AgingMice = read_csv(file = "data_AgingMice/3_data_for_plots/normalized_proteome_AgingMice.csv")

Sys.time()
normalized_matrix_AgingMice = long_to_matrix(normalized_df_AgingMice, 
                                             feature_id_col = 'peptide_group_label',
                                             measure_col = 'Intensity', 
                                             sample_id_col = 'FullRunName')
Sys.time()

pvca_naOmit <- prepare_PVCA_df(normalized_matrix_AgingMice, 
                         sample_annotation_AgingMice, 
                         fill_the_missing = NULL,
                         sample_id_col = 'FullRunName',
                         feature_id_col = 'peptide_group_label',
                         technical_factors = c('MS_batch', "digestion_batch"),
                         biological_factors = c('EarTag', 'Strain', "Diet", "Sex", "Age_Days"))
write_csv(pvca_naOmit, 
          path = "plots_AgingMice/interim_data_for_plots/pvca_naOmit_normalized.csv")
Sys.time()

peptide_count = normalized_df_AgingMice %>% 
  group_by(peptide_group_label) %>%
  summarise(n = sum(!is.na(Intensity))) 

peptides_complete_90 = peptide_count %>%
  filter(n/nrow(sample_annotation_AgingMice) > 0.9) %>%
  pull(peptide_group_label)

normalized_matrix_AgingMice_90 = long_to_matrix(normalized_df_AgingMice %>%
                                                  filter(peptide_group_label %in% peptides_complete_90), 
                                                feature_id_col = 'peptide_group_label',
                                                measure_col = 'Intensity', 
                                                sample_id_col = 'FullRunName')
pvca_90 <- prepare_PVCA_df(normalized_matrix_AgingMice_90, sample_annotation_AgingMice, 
                     sample_id_col = 'FullRunName',
                     feature_id_col = 'peptide_group_label',
                     technical_factors = c('MS_batch', "digestion_batch"),
                     biological_factors= c('EarTag', 'Strain', "Diet", "Sex", "Age_Days"))
# save correlation gg data 
write_csv(pvca_90, 
          path = "plots_AgingMice/interim_data_for_plots/pvca_normalized_90.csv")
dim(normalized_matrix_AgingMice_90)
Sys.time()


peptides_complete_70 = peptide_count %>%
  filter(n/nrow(sample_annotation_AgingMice) > 0.7) %>%
  pull(peptide_group_label)

normalized_matrix_AgingMice_70 = long_to_matrix(normalized_df_AgingMice %>%
                                                  filter(peptide_group_label %in% peptides_complete_70), 
                                                feature_id_col = 'peptide_group_label',
                                                measure_col = 'Intensity', 
                                                sample_id_col = 'FullRunName')
pvca_70 <- prepare_PVCA_df(normalized_matrix_AgingMice_70, sample_annotation_AgingMice, 
                     sample_id_col = 'FullRunName',
                     feature_id_col = 'peptide_group_label',
                     technical_factors = c('MS_batch', "digestion_batch"),
                     biological_factors= c('EarTag', 'Strain', "Diet", "Sex", "Age_Days"))
# save correlation gg data 
write_csv(pvca_70, 
          path = "plots_AgingMice/interim_data_for_plots/pvca_normalized_70.csv")

peptides_complete_50 = peptide_count %>%
  filter(n/nrow(sample_annotation_AgingMice) > 0.5) %>%
  pull(peptide_group_label)
dim(normalized_matrix_AgingMice_70)
Sys.time()


normalized_matrix_AgingMice_50 = long_to_matrix(normalized_df_AgingMice %>%
                                                  filter(peptide_group_label %in% peptides_complete_50), 
                                                feature_id_col = 'peptide_group_label',
                                                measure_col = 'Intensity', 
                                                sample_id_col = 'FullRunName')
pvca_50 <- prepare_PVCA_df(normalized_matrix_AgingMice_50, sample_annotation_AgingMice, 
                     sample_id_col = 'FullRunName',
                     feature_id_col = 'peptide_group_label',
                     technical_factors = c('MS_batch', "digestion_batch"),
                     biological_factors= c('EarTag', 'Strain', "Diet", "Sex", "Age_Days"))
# save correlation gg data 
write_csv(pvca_50, 
          path = "plots_AgingMice/interim_data_for_plots/pvca_normalized_50.csv")
dim(normalized_matrix_AgingMice_50)
Sys.time()