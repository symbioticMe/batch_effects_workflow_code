
#load the libraries
library(readr)
library(tidyverse)
library(reshape2)

spike_ins_in_samples = read_csv('data_InterLab/2_interim_data/spike_ins_in_samples.csv')

#load the data raw (non-normalized) peptide data
peptide_df_raw = read_csv("data_InterLab/2_interim_data/peptide_df_raw.csv") 
spike_in_quantities_raw = peptide_df_raw %>%
  filter(grepl('AQUA', protein_id)) %>%
  mutate(Step = 'Raw',
         filename_new =  gsub('\\.mzXML\\.gz','', run_id)) 

#load the data normalized (median centered) peptide data
peptide_df_medianCentered = read_csv("data_InterLab/2_interim_data/peptide_df_medianCentered.csv") 
spike_in_quantities_normalized = peptide_df_medianCentered %>%
  filter(grepl('AQUA', protein_id))%>%
  mutate(Step = 'Normalized',
         filename_new =  gsub('\\.mzXML\\.gz','', run_id))


spike_in_measurement = rbind(spike_in_quantities_raw, 
                             spike_in_quantities_normalized) %>%
  rename(Intensity = peptide_intensity) %>%
  select(-concentration)
spike_in_measurement$Step = factor(spike_in_measurement$Step, 
                                   levels = c("Raw", "Normalized"))

#TODO: check if the harmonization of the filenames is required
spike_in_quantities = spike_in_measurement %>%
  merge(spike_ins_in_samples, 
        by.x = c('peptide_sequence', 'filename_new'),
        by.y = c('Sequence', 'filename_new'), all.x = T)

spike_in_quantities = spike_in_quantities %>%
  group_by(peptide_id, concentration, site, Step, peptide_sequence) %>% 
  summarize(Intensity = mean(Intensity))

write_csv(spike_in_quantities,'data_InterLab/3_data_for_plots/spike_in_quantities.csv')
