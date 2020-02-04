#!/usr/bin/env Rscript --vanilla

#load the libraries
library(readr)
library(dplyr)
library(stringr)

#load the given annotation
sample_annotation_original <- read_delim('data_PanCancer/1_original_data/annotation_col_anova_PCP_final.txt', delim = '\t')

#adapt sample annotation
sample_annotation_PanCancer = sample_annotation_original %>% 
  mutate(X = NULL, X.1 = NULL, X.2 = NULL, X.3 = NULL) %>%
  mutate(Condition = gsub('[0-9]+_','', sample),
         BioReplicate = str_extract(sample, '^[0-9]+'),
         tissue = gsub('([a-z]+\\.)|([A-Z]+_)','', Condition),
         case_control = ifelse(grepl('ctrl', Condition), 'control', 'case')) %>%
  mutate(tissue = gsub('ctrl_', '', tissue)) %>%
  mutate(tissue = gsub( "-.*$", "", tissue )) %>%
  mutate(replicate = ifelse(grepl('-R[0-9]$', Condition), gsub('.*-(R[0-9])$', '\\1', Condition), NA)) %>%
  mutate(sample_id_final = paste(tissue, case_control, sep = '_')) %>%
  group_by(replicate, case_control, tissue) %>% 
  mutate(id = row_number()) %>% ungroup() %>%
  mutate(sample_id_final = ifelse(is.na(replicate), paste(sample_id_final, id, sep = '_'), 
                                  paste(sample_id_final, replicate, sep = '_'))) %>%
  mutate(sample_name = paste(Condition, BioReplicate, sep = '_'))
sample_annotation_PanCancer$patient_ID = sample_annotation_PanCancer %>% group_indices(sample_id_final)
sample_annotation_PanCancer$patient_ID = paste('patient_', 
                                                sprintf("%03d", sample_annotation_PanCancer$patient_ID ), sep = '')


write_csv(sample_annotation_PanCancer, 
          path = 'data_PanCancer/3_data_for_plots/sample_annotation_PanCancer.csv')
