#!/usr/bin/env Rscript --vanilla

#load the libraries
library(proBatch)
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

#save new data frame
write_delim(fragmentome_centered_un_log, 
            path = "data_InterLab/2_interim_data/all_sites_global_q_001_applied_to_local_global_medianCentered.tsv",
            delim = '\t')
