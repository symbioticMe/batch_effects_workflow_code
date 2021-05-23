############### proBatch analysis on Fabian Wendt & Sandra Goetze DDA data #################################
# biological (animal, time, condition) and technical (plate_row, plate_col, ms_order (given as YYMM_SEQ))
# conducted quantile normalization + combat (by robot batch)
# Chloe H. Lee, 05/01/2021
############################################################################################################

library(tidyverse)
library(proBatch)
library(ggplot2)
library(grid)
library(gridExtra)
setwd('C:/Users/Chloe/Desktop/ETH Zurich/Manuscript/DDA_data')
source("wes_palettes.R")


############################### load data ###############################################################
# Note: I relabelled batch as MS_batch for simplicity of initial analysis but probably will need to relabel it (e.g. robot_batch)
original_mtx <- read.csv('Bariatric/20210413_Progenesis_grouped-conditions_peptide-output.csv', skip = 2)
quantData <- read.csv("Bariatric/pep_mtx_raw.csv") %>% dplyr::rename(peptide_group_label = FullPeptideName)
sample_annotation <- read.csv("Bariatric/meta.csv") %>% dplyr::rename(FullRunName = X, robot_batch = batch) %>%
  mutate(order = as.numeric(order)) %>% mutate(time = as.factor(time)) 
peptide_annotation <- original_mtx %>% dplyr::select(Sequence, Modifications, Accession) %>% 
  mutate(peptide_group_label = paste0(Sequence, Modifications))
length(intersect(peptide_annotation$FullPeptideName, quantData$FullPeptideName)) #all match corresponding to quantData 

# data matrix
data_matrix <- quantData %>% column_to_rownames(var = 'peptide_group_label') %>% as.matrix()

# color annotation table 
color_list <- sample_annotation_to_colors(sample_annotation,
                                          factor_columns = c('animal', 'time', 'plate_row', 'plate_col', 'robot_batch', 'condition'),
                                          numeric_columns = c('order'))


############################### data analysis ###############################################################
# quantile normalization 
quantile_normalized_matrix = normalize_data_dm(data_matrix, normalize_func = 'quantile')
quantile_normalized_long <- matrix_to_long(quantile_normalized_matrix) %>% 
  mutate_if(is.factor, as.character)

# filter peptides with many missing values for combat 
quantile_normalized_matrix_for_combat <- quantile_normalized_matrix[rowSums(is.na(quantile_normalized_matrix)) <= 10,] 
quantile_normalized_matrix_for_combat[is.na(quantile_normalized_matrix_for_combat)] = -5.5
quantile_normalized_long_for_combat <- matrix_to_long(quantile_normalized_matrix_for_combat) %>%  mutate_if(is.factor, as.character)

# combat correction 
combat_corrected_lg_batch <- correct_with_ComBat_df(quantile_normalized_long_for_combat, sample_annotation,  batch_col = "robot_batch")
combat_corrected_matrix_batch <- long_to_matrix(combat_corrected_lg_batch,   measure_col = 'Intensity')

# ratio (normalized to t0 per animal)
time0_animals <- sample_annotation %>%   mutate(time = as.numeric(as.character(time))) %>% filter(time == 0) %>% pull(animal) %>% unique

quantile_normalized_lg_ratio <- quantile_normalized_long_for_combat %>% 
  left_join(sample_annotation, by = 'FullRunName') %>%
  mutate(time = as.numeric(as.character(time))) %>%
  filter(animal %in% time0_animals) %>%
  group_by(animal, peptide_group_label) %>%
  mutate(dIntensity = Intensity - Intensity[time == 0]) %>% 
  filter(time != 0)
combat_corrected_lg_batch_ratio <- combat_corrected_lg_batch %>% 
  left_join(sample_annotation, by = 'FullRunName') %>%
  mutate(time = as.numeric(as.character(time))) %>%
  filter(animal %in% time0_animals) %>%
  group_by(animal, peptide_group_label) %>%
  mutate(dIntensity = Intensity - Intensity[time == 0]) %>% 
  filter(time != 0)
quantile_normalized_mtx_ratio <- long_to_matrix(quantile_normalized_lg_ratio,     measure_col = 'dIntensity')
combat_corrected_mtx_batch_ratio <- long_to_matrix(combat_corrected_lg_batch_ratio,     measure_col = 'dIntensity')



############################### figures for manuscirpt #############################################
# 1. inter / intra on raw (non-normalized data) vs. quantile normalized. Not on ratio data.
# 2. PCA colored by robot batch (or robot row) after quantile vs after combat. Not on ratio data.
# 3. Hierachical clustering of ratio data before any correction and after all corrections. Only show animal.


#### 1. peptide correlation (using all peptides)
p1 <- plot_peptide_corr_distribution(data_matrix,  peptide_annotation,
                                     protein_col = 'Accession',  plot_title = 'Log transformed') + 
  theme(axis.text.x = element_text(size = 12)) 

p2 <- plot_peptide_corr_distribution(quantile_normalized_matrix,  peptide_annotation,
                                     protein_col = 'Accession',  plot_title = 'Quantile normalized') + 
  theme(axis.text.x = element_text(size = 12)) 
grid.arrange(p1, p2, ncol = 2)


##### 2. PCA plot for robot batch (where peptides missing in >10 samples removed, and NA replaced with -5.5)
pca1 = plot_PCA(quantile_normalized_matrix_for_combat, sample_annotation, color_by = 'robot_batch',
                plot_title = 'Quantile normalized', color_scheme = color_list[['robot_batch']]) +
  theme(plot.title = element_text(hjust = 0.5)) 
pca2 = plot_PCA(combat_corrected_matrix_batch, sample_annotation, color_by = 'robot_batch',
                plot_title = 'ComBat corrected by robot batch', color_scheme = color_list[['robot_batch']]) +
  theme(plot.title = element_text(hjust = 0.5)) 
grid.arrange(pca1, pca2, ncol = 2)



####### hierarchical clustering on ratio data 
selected_annotations <- c('animal')
plot_heatmap_diagnostic(quantile_normalized_mtx_ratio, sample_annotation,
                        factors_to_plot = selected_annotations,
                        cluster_cols = TRUE,
                        color_list = color_list,
                        show_rownames = FALSE, show_colnames = FALSE)

plot_heatmap_diagnostic(combat_corrected_mtx_batch_ratio, sample_annotation,
                        factors_to_plot = selected_annotations,
                        cluster_cols = TRUE,
                        color_list = color_list,
                        show_rownames = FALSE, show_colnames = FALSE)



