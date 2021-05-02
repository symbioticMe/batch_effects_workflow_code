#!/usr/bin/env Rscript --vanilla

#load the libraries
library(tidyverse)
library(proBatch)

selected_peptide = "46213_NVGVSFYADKPEVTQEQK_2"

example_peptide_annotation %>% filter(Gene == 'Haao') %>% pull(peptide_group_label)

# load the data
sample_annotation_AgingMice = read_csv("Aging_mice/data_Aging_mice/2_interim_data_AgingMice//sample_annotation_Aging_mice.csv")
normalized_df_AgingMice = read_csv(file = "Aging_mice/data_Aging_mice/2_interim_data_AgingMice/quantile_normalized_matrix_agingMice.csv")
color_annotation_AgingMice = readRDS("Aging_mice/data_Aging_mice/2_interim_data_AgingMice/color_annotation.rda")
peptide_df = read_csv("Aging_mice/data_Aging_mice/1_original_data_AgingMice/peptide_annotation_AgingMice.csv")


# prepare color annotation
color_list_AgingMice <- color_annotation_AgingMice$list_of_colors

# plot principal componenet analysis
spikein_bovine_A1ag_plot <- plot_spike_in(normalized_df_AgingMice, color_annotation_AgingMice, 
                                          peptide_annotation = peptide_annotation,
                                          protein_col = 'Gene', spike_ins = "BOVINE_A1ag", 
                                          plot_title = 'Spike-in BOVINE protein peptides',
                                          color_by_batch = TRUE, color_scheme = color_list_AgingMice[["MS_batch"]])


# save figure
saveRDS(spikein_bovine_A1ag_plot, 
        file = "Aging_mice/graphs_Aging_mice/interim_ggplot_objects_AgingMice/4b_spikein_bovine_A1ag_plot.rds")
source('lib/helpers.R')
ggsave(spikein_bovine_A1ag_plot + theme_publication(),
       file = "Aging_mice/graphs_Aging_mice/4b_spikein_bovine_A1ag_plot.png")






