#!/usr/bin/env Rscript --vanilla

#load the libraries
library(tidyverse)
library(proBatch)

#load the data
sample_annotation_original = read_csv("data_AgingMice/1_original_data/sample_annotation_E1891171630.csv")

#prepare annotation: infer sample order, convert date Time columns to PoSIX etc
sample_annotation_AgingMice = date_to_sample_order(sample_annotation_AgingMice_publication,
                                          time_column = c('RunDate','RunTime'),
                                          new_time_column = 'DateTime',
                                          dateTimeFormat = c("%b_%d", "%H:%M:%S"),
                                          new_order_col = 'order',
                                          instrument_col = NULL)
sample_annotation_AgingMice$Age_Days = as.numeric(sample_annotation_AgingMice$Age_Days)

#save the data frame
write_csv(sample_annotation_AgingMice, 'data_AgingMice/3_data_for_plots/sample_annotation_AgingMice.csv')
