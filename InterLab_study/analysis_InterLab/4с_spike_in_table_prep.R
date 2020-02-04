#!/usr/bin/env Rscript --vanilla

#load the libraries
library(readr)
library(proBatch)
library(tidyverse)
library(reshape2)

#load the annotation
sample_annotation_InterLab = read_csv('data_InterLab/3_data_for_plots/sample_annotation_InterLab.csv')


#load the spike-in concentrations
#load the list of spike-in peptides and corresponding groups (MOESM4)
spike_in_group_df  = read_csv('data_InterLab/1_original_data/spike_in_group.csv')

spike_in_group_df$Sequence = gsub('\\[\\+(10|08)\\]', '', spike_in_group_df$Peptide)

#load the list of spike-in peptide dilution series (is it also one of the supplementary tables?)
dilution_groups_df = read_csv('data_Interlab/1_original_data/group_to_dilutionSeries.csv')

#load the list of dilution series to groups mapped (MOESM3)
##TOREMOVE:
##dilution_groups_df = group_to_dilutionSeries
dilution_groups_df = melt(dilution_groups_df, id.vars = 'Group',
                          measure.vars = paste('S', 1:5, sep = ''),
                          value.name = 'concentration',
                          variable.name = 'dilution_series')

#merge it in one super-table:
dilution_groups_df = dilution_groups_df %>% 
  merge(spike_in_group_df, by = 'Group')
spike_ins_in_samples = sample_annotation_InterLab %>%
  merge(dilution_groups_df, by = 'dilution_series', all = T)

write_csv(spike_ins_in_samples,'data_Interlab/2_interim_data_InterLab/spike_ins_in_samples.csv')

