#!/usr/bin/env Rscript --vanilla

#load the libraries
library(readr)
library(tidyverse)

#load the data
good_columns = c('filename')
cols_to_get <- rep(list(col_guess()), length(good_columns))
names(cols_to_get) <- good_columns
cols_to_get2 = do.call(cols_only, cols_to_get)

proteome = read_delim("data_InterLab/1_original_data/all_sites_global_q_0.01_applied_to_local_global.txt",
                      delim = '\t', escape_double = FALSE, trim_ws = TRUE, 
                      col_types = cols_to_get2)

#create the sample annotation
sample_annotation = proteome %>%
  select(filename) %>% distinct() %>%
  mutate(filename_new = gsub('/scratch/[0-9]+\\.tmpdir/','', filename)) %>%
  mutate(filename_new = gsub('\\.mzXML\\.gz','', filename_new)) %>% 
  mutate(site = tolower(str_extract(filename, '[sS]ite[0-9]+')),
         dilution_series = gsub('(.*_|.*DynRange_?|.*_DynRge)(S[1-5])(_.*|[A-Z].*)', 
                                '\\2', filename),
         day = gsub('.*([dD]ay_?[0-9]+).*','\\1', filename))

sample_annotation = as.data.frame(sample_annotation)

#save the data frame
write_csv(sample_annotation, 'data_InterLab/3_data_for_plots/sample_annotation_InterLab.csv')