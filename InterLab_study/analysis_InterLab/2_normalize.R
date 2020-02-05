#!/usr/bin/env Rscript --vanilla

#load the libraries
library(proBatch3.4)
library(readr)
library(dplyr)
library(tidyr)
source('transition_fragment_conversion.R')

#load the data
essential_columns = c('peptide_group_label','filename','aggr_Fragment_Annotation', 'aggr_Peak_Area')
cols_to_get <- rep(list(col_guess()), length(essential_columns))
names(cols_to_get) <- essential_columns
cols_to_get2 = do.call(cols_only, cols_to_get)

proteome = read_delim("data_InterLab/1_original_data/all_sites_global_q_0.01_applied_to_local_global.txt", 
                      delim =  "\t", escape_double = FALSE, trim_ws = TRUE, 
                      col_types = cols_to_get2)

#transform transition-level matrix to fragment-level
fragmentome <- proteome %>% 
  transition_to_fragment()
fragmentome = fragmentome %>%
  log_transform_df(measure_col = 'Ion_intensity')

#normalize the samples
fragmentome_centered = normalize_sample_medians_df(fragmentome, 
                                                sample_id_col = 'filename',
                                                measure_col = 'Ion_intensity')
#exponentiate (reverse log-transformation)
fragmentome_centered_un_log = fragmentome_centered %>% 
  unlog_df(measure_col = 'Ion_intensity')

print(names(fragmentome_centered_un_log))


proteome_median_centered = fragment_df_to_openSWATH(fragmentome_centered_un_log, 
                                                    fragment_intensity_column = 'Ion_intensity', 
                                                    fragment_annotation_column = 'Ion_ID',
                                                    id_column = 'ions',
                                                    fragment_united_column = 'aggr_Fragment_Annotation_new',
                                                    fragment_united_int_column = 'aggr_Peak_Area_new',
                                                    un_log = NULL, 
                                                    intensities_to_exclude = c('Ion_intensity_log2', 'Ion_intensity', 'Intensity_normalized',
                                                                               'diff','median_global','median_run'))
old_names <- c("aggr_Peak_Area", "aggr_Fragment_Annotation", 
               "aggr_Peak_Area_new", "aggr_Fragment_Annotation_new")
new_names <-  c("aggr_Peak_Area_old", "aggr_Fragment_Annotation_old", 
                "aggr_Peak_Area", "aggr_Fragment_Annotation")

supporting_info_cols = c("transition_group_id", "Sequence", "FullPeptideName", "RT", 
                         "assay_rt", "Intensity", "ProteinName", "m_score", "run_id", 
                         "peak_group_rank", "Charge", "decoy", 'peptide_group_label','filename')
cols_to_get <- rep(list(col_guess()), length(supporting_info_cols))
names(cols_to_get) <- supporting_info_cols
cols_to_get2 = do.call(cols_only, cols_to_get)
proteome_supporting_info = read_delim('data_InterLab/1_original_data/all_sites_global_q_0.01_applied_to_local_global.txt',
                                      delim =  "\t", escape_double = FALSE, trim_ws = TRUE, 
                                      col_types = cols_to_get2)

proteome_median_centered = proteome_median_centered %>%
  rename_at(vars(old_names), ~ new_names) %>%
  merge(proteome_supporting_info, by = c('peptide_group_label','filename'))

#save new data frame
write_delim(proteome_median_centered, 
            path = "data_InterLab/2_interim_data/all_sites_global_q_001_applied_to_local_global_medianCentered.tsv",
            delim = '\t')
