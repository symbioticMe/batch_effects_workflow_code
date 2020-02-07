#!/usr/bin/env Rscript --vanilla

#load the libraries
library(readr)
library(proBatch)
library(RColorBrewer)

#load the data
sample_annotation_AgingMice = read_csv('data_AgingMice/3_data_for_plots/sample_annotation_AgingMice.csv')

#define the colors
color_list_AgingMice = sample_annotation_to_colors(sample_annotation_AgingMice,
                                                   sample_id_col = 'FullRunName',
                                                   factor_columns = c('MS_batch','EarTag', "Strain", "Diet", "Sex", 'digestion_batch'),
                                                   numeric_columns = c("Age_Days", "order", 'DateTime'),
                                                   rare_categories_to_other = T,
                                                   numeric_palette_type = 'brewer')

color_list_AgingMice$Diet = c(HF = 'gold', CD = 'forestgreen', CDHFD = 'grey')
names_dig_batch = names(color_list_AgingMice$digestion_batch)
color_list_AgingMice$digestion_batch = brewer.pal(n = length(color_list_AgingMice$digestion_batch), name = 'Dark2')
names(color_list_AgingMice$digestion_batch) = names_dig_batch
color_list_AgingMice$EarTag['ET1506'] = 'forestgreen'
color_list_AgingMice$EarTag['ET1524'] = 'darkred'
color_list_AgingMice$MS_batch['Batch_7'] = 'purple'
#save the color list & data.frame to cache
saveRDS(color_list_AgingMice, file="data_AgingMice/3_data_for_plots/color_annotation.rda")