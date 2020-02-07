#!/usr/bin/env Rscript --vanilla

#load the libraries
library(tidyverse)
library(proBatch)

#load the data
#TODO: update "m_score" according to proper JPP outputpro
good_columns =c("peptide_group_label", "filename", "Intensity", "m_score",'ProteinName')
cols_to_get <- rep(list(col_guess()), length(good_columns))
names(cols_to_get) <- good_columns
cols_to_get2 = do.call(cols_only, cols_to_get)

proteome_df_AgingMice = read_delim("data_AgingMice/1_original_data/E1801171630_feature_alignment_requant.tsv",
                      "\t", escape_double = FALSE, trim_ws = TRUE, col_types = cols_to_get2)
#TODO: change for JPP version of the data

##load the sample annotation
sample_annotation_AgingMice_publication = read_csv("data_AgingMice/3_data_for_plots/sample_annotation_AgingMice.csv")


proteome_df_AgingMice = proteome_df_AgingMice %>%
  filter(!grepl('DECOY', ProteinName)) %>%
  select(-ProteinName) %>%
  mutate(FullRunName = gsub('.+tmpdir/wevan_(.+)\\.mzXML\\.gz', '\\1', filename)) %>% 
  filter(FullRunName %in% sample_annotation_AgingMice_publication$FullRunName) #12,168,375 measurements

proteome_df_AgingMice = log_transform_df(proteome_df_AgingMice, 
                                                   measure_col = 'Intensity',
                                                   log_base = 2, offset = 1)

write_csv(proteome_df_AgingMice, 'data_AgingMice/2_interim_data/raw_proteome_AgingMice_with_requants.csv')


#remove imputed values:
proteome_df_AgingMice = proteome_df_AgingMice %>% 
  mutate(Intensity = ifelse(m_score == 2, NA, Intensity))

write_csv(proteome_df_AgingMice, 'data_AgingMice/2_interim_data/raw_proteome_AgingMice.csv')



