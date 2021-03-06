#!/usr/bin/env Rscript --vanilla


#load the libraries
library(proBatch)
library(ggplot2)
library(dplyr)
library(readr)

proteome_df_AgingMice = read_csv("data_AgingMice/2_interim_data/raw_proteome_AgingMice.csv")
peptide_df_AgingMice = read_csv("data_AgingMice/1_original_data/peptide_annotation.csv")

proteome_df_AgingMice_50 = proteome_df_AgingMice %>%
  group_by(peptide_group_label) %>%
  mutate(n_na = sum(is.na(Intensity))) %>%
  filter(n_na < 0.5 * 375)
  
proteome_log_AgingMice = long_to_matrix(proteome_df_AgingMice_50, 
                                        feature_id_col = 'peptide_group_label',
                                        measure_col = 'Intensity', 
                                        sample_id_col = 'FullRunName')

# Peptide correlation of raw matrix
peptide_corr_raw <-calculate_peptide_corr_distr(data_matrix = proteome_log_AgingMice, 
                                                      peptide_annotation = peptide_df_AgingMice,
                                                      protein_col = 'Gene',
                                                      feature_id_col = 'peptide_group_label')
write_csv(peptide_corr_raw, 'plots_AgingMice/interim_data_for_plots/peptide_cor_raw.csv')
Sys.time()

# peptide_corr_plot = plot_peptide_corr_distribution.corrDF(peptide_corr_raw, 
#                                plot_title = 'Peptide correlation (raw data)',
#                                filename = 'Aging_mice/graphs_Aging_mice/6b_peptide_correlation_raw.pdf')
