############### proBatch analysis on Fabian Wendt & Sandra Goetze DDA data #################################
# biological (animal, time, condition) and technical (plate_row, plate_col, ms_order (given as YYMM_SEQ))
# conducted quantile normalization + combat (by batch / plate row)
# Chloe H. Lee, 04/24/2021
############################################################################################################

library(tidyverse)
library(proBatch)
library(ggplot2)
library(grid)
library(gridExtra)
setwd('C:/Users/Chloe/Desktop/ETH Zurich/Manuscript/DDA_data')
source("wes_palettes.R")
plot_PVCA.df <- function(pvca_res,  colors_for_bars = NULL,  filename = NULL, width = NA, height = NA, 
                         units = c('cm','in','mm'), plot_title = NULL,  theme = 'classic',
                         base_size = 20){
  pvca_res = pvca_res %>%
    mutate(label = factor(label, levels=label))
  
  y_title = 'Weighted average proportion variance'
  gg  = ggplot(pvca_res, aes(x = label, y = weights, fill = category))+
    geom_bar(stat = 'identity', color = 'black')+
    ylab(y_title)
  
  
  if(is.null(colors_for_bars)){
    colors_for_bars = c('grey', wes_palettes$Rushmore[3:5])
    names(colors_for_bars) = c('residual', 'biological', 
                               'biol:techn', 'technical')
    
  } else {
    if (length(colors_for_bars) != 4){
      color_names = paste(c('residual', 'biological', 'biol:techn', 
                            'technical'), collapse = ' ')
      warning(sprintf('four colors for: %s were expected', color_names))
    }
  }
  gg = gg + scale_fill_manual(values = colors_for_bars)
  
  if (!is.null(plot_title)){
    gg = gg + ggtitle(plot_title)
  }
  
  #Change the theme
  if(!is.null(theme) && theme == 'classic'){
    gg = gg + theme_classic(base_size = base_size)
  }else{
    message("plotting with default ggplot theme, only theme = 'classic' 
            implemented")
  }
  
  gg = gg +
    theme(axis.title.x = NULL, 
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))+
    xlab(NULL)+
    theme(text = element_text(size=15))+
    guides(fill=guide_legend(override.aes=list(color=NA), title=NULL))
  
  return(gg)
}


#################### laod data ##################################
# Note: I relabelled batch as MS_batch for simplicity of initial analysis but probably will need to relabel it (e.g. robot_batch)
original_mtx <- read.csv('Bariatric/20210413_Progenesis_grouped-conditions_peptide-output.csv', skip = 2)
quantData <- read.csv("Bariatric/pep_mtx_raw.csv") %>% dplyr::rename(peptide_group_label = FullPeptideName)
sample_annotation <- read.csv("Bariatric/meta.csv") %>% dplyr::rename(FullRunName = X, MS_batch = batch) %>%
  mutate(order = as.numeric(order)) %>% mutate(time = as.factor(time))
peptide_annotation <- original_mtx %>% dplyr::select(Sequence, Modifications, Accession) %>% 
  mutate(peptide_group_label = paste0(Sequence, Modifications))
length(intersect(peptide_annotation$FullPeptideName, quantData$FullPeptideName)) #all match corresponding to quantData 

# data matrix
data_matrix <- quantData %>% column_to_rownames(var = 'peptide_group_label') %>% as.matrix()

# color annotation table 
color_list <- sample_annotation_to_colors(sample_annotation,
                                          factor_columns = c('animal', 'time', 'plate_row', 'plate_col', 'MS_batch', 'condition'),
                                          numeric_columns = c('order'))



#################### illustrating experimental design ######################################
a <- sample_annotation %>% 
  ggplot(aes_string(x = 'order', y = 'animal')) + 
  geom_point(aes(colour = factor(plate_row)), size = 4) + theme_bw()

b <- sample_annotation %>% 
  ggplot(aes_string(x = 'order', y = 'time')) + 
  geom_point(aes(colour = factor(MS_batch)), size = 4) + theme_bw()

c <- sample_annotation %>% 
  ggplot(aes_string(x = 'order', y = 'condition')) + 
  geom_point(aes(colour = factor(plate_col)), size = 4) + theme_bw()

grid.arrange(a, c, b)


###################### Initial assessment ##################################
plot_sample_mean(data_matrix, sample_annotation, order_col = 'order',
                 batch_col = 'MS_batch', color_by_batch = TRUE, ylimits = c(3, 7),
                 color_scheme = color_list[['MS_batch']], vline_color = NULL)

log_transformed_long <- matrix_to_long(data_matrix)
batch_col = 'MS_batch'
plot_boxplot(log_transformed_long, sample_annotation, batch_col = batch_col, color_scheme = color_list[[batch_col]]) +
  theme(axis.text.x = element_text(size = 8))



################## Normalization and assessment of normalized matrix ############################
quantile_normalized_matrix = normalize_data_dm(data_matrix, normalize_func = 'quantile')
quantile_normalized_long <- matrix_to_long(quantile_normalized_matrix) %>% 
  mutate_if(is.factor, as.character)

plot_sample_mean(quantile_normalized_matrix, sample_annotation,  color_by_batch = TRUE, ylimits = c(5, 6),
                 color_scheme = color_list[['MS_batch']])
plot_boxplot(quantile_normalized_long, sample_annotation, batch_col = batch_col, color_scheme = color_list[[batch_col]]) +
  theme(axis.text.x = element_text(size = 8))



################ Diagnostics on normalized matrix ###################################################
# hierarchical clustering 
selected_annotations <- c('animal', 'time', 'condition', 'plate_row', "plate_col" , 'order', 'MS_batch')
plot_hierarchical_clustering(quantile_normalized_matrix,
                             sample_annotation = sample_annotation,
                             color_list = color_list,
                             factors_to_plot = selected_annotations,
                             distance = 'euclidean', agglomeration = 'complete',
                             label_samples = FALSE)

# heatmap 
plot_heatmap_diagnostic(quantile_normalized_matrix, sample_annotation,
                        factors_to_plot = selected_annotations,
                        cluster_cols = TRUE,
                        color_list = color_list,
                        show_rownames = FALSE, show_colnames = FALSE)


# pvca 
# Note: due to low # of samples for animal:time, plate_col:plate_row, I computed PVCA for shorlisted covariates and combined variances 
# from multiple computations, averaging overlapping variances. Please refer to pvca_analysis_for_dda_data for details. 
pvca_analysis_for_dda_data <- function(matrix){
  
  technical_factors_1 = c('MS_batch', 'plate_col')
  technical_factors_2 = c('MS_batch', 'plate_row')
  biological_factors_1 = c('animal', 'condition')
  biological_factors_2 = c('time', 'condition')
  
  pvca_1 <- plot_PVCA(matrix, sample_annotation,
                      technical_factors = technical_factors_1,   biological_factors = biological_factors_1)
  pvca_2 <- plot_PVCA(matrix, sample_annotation,
                      technical_factors = technical_factors_1,   biological_factors = biological_factors_2)
  pvca_3 <- plot_PVCA(matrix, sample_annotation,
                      technical_factors = technical_factors_2,   biological_factors = biological_factors_1)
  pvca_4 <- plot_PVCA(matrix, sample_annotation,
                      technical_factors = technical_factors_2,   biological_factors = biological_factors_2)
  
  pvca_averaged_weights <- pvca_1$data %>% rbind(pvca_2$data) %>% rbind(pvca_3$data) %>% rbind(pvca_4$data) %>%
    group_by(label, category) %>% dplyr::summarise(weights = mean(weights)) %>% 
    mutate_if(is.factor, as.character) %>% dplyr::arrange(-weights) 
  labels <- as.character(pvca_averaged_weights$label)
  labels <- c(labels[!labels %in% c('resid', 'Below 1%')], c( 'Below 1%', 'resid'))
  pvca_averaged_weights <- pvca_averaged_weights %>% mutate(label = factor(label, levels = labels))
  
  plot_PVCA.df(pvca_averaged_weights) 
}
pvca_analysis_for_dda_data(quantile_normalized_matrix)



######################### Batch correction ######################################################
###### analyze peptides with missing values 
length(which(rowSums(is.na(quantile_normalized_matrix)) > 0)) #469/1987 peptides contain missing values 

no_na_per_peptides <- rowSums(is.na(quantile_normalized_matrix)) 
hist(no_na_per_peptides/ncol(quantile_normalized_matrix), main = '% of missing values per peptide')
hist(no_na_per_peptides[which(no_na_per_peptides != 0)]/ncol(quantile_normalized_matrix), main = '% of missing values (w/o 0) per peptide')

hist(no_na_per_peptides, main = '# of missing values per peptide')
hist(no_na_per_peptides[which(no_na_per_peptides != 0)], main = '# of missing values (w/o 0) per peptide')



###### Correction by ComBat
# data with peptides > 10 missing values removed and replace NA with -5.5 (suggested by Patrick)
quantile_normalized_matrix_filterna <- quantile_normalized_matrix[rowSums(is.na(quantile_normalized_matrix)) <= 10,] 
quantile_normalized_matrix_filterna[is.na(quantile_normalized_matrix_filterna)] = -5.5
quantile_normalized_long_for_combat <- matrix_to_long(quantile_normalized_matrix_filterna) %>%  mutate_if(is.factor, as.character)

# by robot batch 
combat_corrected_lg_batch <- correct_with_ComBat_df(quantile_normalized_long_for_combat, sample_annotation,  batch_col = "MS_batch")
combat_corrected_matrix_batch <- long_to_matrix(combat_corrected_lg_batch,   measure_col = 'Intensity')

# by plate row
combat_corrected_lg_platerow <- correct_with_ComBat_df(quantile_normalized_long_for_combat, sample_annotation,  batch_col = "plate_row")
combat_corrected_matrix_platerow <- long_to_matrix(combat_corrected_lg_platerow,   measure_col = 'Intensity')




##### median centering 
# data with peptides > 10 missing values removed 
quantile_normalized_matrix_filterna <- quantile_normalized_matrix[rowSums(is.na(quantile_normalized_matrix)) <= 10,] 
quantile_normalized_long_for_medcenter <- matrix_to_long(quantile_normalized_matrix_filterna) %>%  mutate_if(is.factor, as.character)


# by robot batch 
medcenter_corrected_lg_batch <- center_feature_batch_medians_df(quantile_normalized_long_for_medcenter, sample_annotation,
                                                  batch_col = "MS_batch", keep_all = 'default')
medcenter_corrected_matrix_batch <- long_to_matrix(medcenter_corrected_lg_batch,     measure_col = 'Intensity')

# by plate row 
medcenter_corrected_lg_platerow <- center_feature_batch_medians_df(quantile_normalized_long_for_medcenter, sample_annotation,
                                                          batch_col = "plate_row", keep_all = 'default')
medcenter_corrected_matrix_platerow <- long_to_matrix(medcenter_corrected_lg_platerow,     measure_col = 'Intensity')



#################### prepare data normalized to time 0 #####################################
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
combat_corrected_lg_platerow_ratio <- combat_corrected_lg_platerow %>% 
  left_join(sample_annotation, by = 'FullRunName') %>%
  mutate(time = as.numeric(as.character(time))) %>%
  filter(animal %in% time0_animals) %>%
  group_by(animal, peptide_group_label) %>%
  mutate(dIntensity = Intensity - Intensity[time == 0]) %>% 
  filter(time != 0)

quantile_normalized_mtx_ratio <- long_to_matrix(quantile_normalized_lg_ratio,     measure_col = 'dIntensity')
combat_corrected_mtx_batch_ratio <- long_to_matrix(combat_corrected_lg_batch_ratio,     measure_col = 'dIntensity')
combat_corrected_mtx_platerow_ratio <- long_to_matrix(combat_corrected_lg_platerow_ratio,     measure_col = 'dIntensity')




###################### evaluate batch corrected matrices ######################################################
# PVCA 
pvca_analysis_for_dda_data(combat_corrected_matrix_batch)
pvca_analysis_for_dda_data(medcenter_corrected_matrix_batch)
pvca_analysis_for_dda_data(combat_corrected_matrix_platerow)
pvca_analysis_for_dda_data(medcenter_corrected_matrix_platerow)


# pheatmap on ratio data 
selected_annotations <- c('animal', 'time', 'condition', 'plate_row', "plate_col" , 'order', 'MS_batch')
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
plot_heatmap_diagnostic(combat_corrected_mtx_platerow_ratio, sample_annotation,
                        factors_to_plot = selected_annotations,
                        cluster_cols = TRUE,
                        color_list = color_list,
                        show_rownames = FALSE, show_colnames = FALSE)


# hierarchical clustering on ratio data 
plot_hierarchical_clustering(combat_corrected_mtx_batch_ratio,
                             sample_annotation = sample_annotation,
                             color_list = color_list,
                             factors_to_plot = selected_annotations,
                             distance = 'euclidean', agglomeration = 'complete',
                             label_samples = FALSE)
plot_hierarchical_clustering(combat_corrected_mtx_platerow_ratio,
                             sample_annotation = sample_annotation,
                             color_list = color_list,
                             factors_to_plot = selected_annotations,
                             distance = 'euclidean', agglomeration = 'complete',
                             label_samples = FALSE)





######################### sample and peptide correlation ###################
# sample correlation - having ms batch as the batch col 
p1 <- plot_sample_corr_distribution(quantile_normalized_mtx_ratio,   sample_annotation,
                                                  batch_col = 'MS_batch',   biospecimen_id_col = 'animal',
                                                  plot_title = 'Quantile normalized',    plot_param = 'batch_replicate') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12)) + ylim(-0.3, 1) + xlab(NULL)
p2 <- plot_sample_corr_distribution(combat_corrected_mtx_batch_ratio,   sample_annotation,
                                                  batch_col = 'MS_batch',   biospecimen_id_col = 'animal',
                                                  plot_title = 'ComBat corrected (by batch)',    plot_param = 'batch_replicate') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12)) + ylim(-0.3, 1)  + xlab(NULL)
p3 <- plot_sample_corr_distribution(combat_corrected_mtx_platerow_ratio,   sample_annotation,
                                    batch_col = 'MS_batch',   biospecimen_id_col = 'animal',
                                    plot_title = 'ComBat corrected (by plate row)',    plot_param = 'batch_replicate') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12)) + ylim(-0.3, 1)  + xlab(NULL)
grid.arrange(p1, p2, p3, ncol = 3)


# sample correlation - having plate row as batch col
p1 <- plot_sample_corr_distribution(quantile_normalized_mtx_ratio,   sample_annotation,
                                    batch_col = 'plate_row',   biospecimen_id_col = 'animal',
                                    plot_title = 'Quantile normalized',    plot_param = 'batch_replicate') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12)) + ylim(-0.3, 1) + xlab(NULL)
p2 <- plot_sample_corr_distribution(combat_corrected_mtx_batch_ratio,   sample_annotation,
                                    batch_col = 'plate_row',   biospecimen_id_col = 'animal',
                                    plot_title = 'ComBat corrected (by batch)',    plot_param = 'batch_replicate') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12)) + ylim(-0.3, 1)  + xlab(NULL)
p3 <- plot_sample_corr_distribution(combat_corrected_mtx_platerow_ratio,   sample_annotation,
                                    batch_col = 'plate_row',   biospecimen_id_col = 'animal',
                                    plot_title = 'ComBat corrected (by plate row)',    plot_param = 'batch_replicate') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12)) + ylim(-0.3, 1)  + xlab(NULL)
grid.arrange(p1, p2, p3, ncol = 3)


# peptide correlation on non-ratio data
p4a <- plot_peptide_corr_distribution(quantile_normalized_matrix_filterna,  peptide_annotation,
                                      protein_col = 'Accession',  plot_title = 'Quantile normalized') + 
  theme(axis.text.x = element_text(size = 12)) 

p5a <- plot_peptide_corr_distribution(combat_corrected_matrix_batch,  peptide_annotation,
                                      protein_col = 'Accession',  plot_title = 'Combat corrected (by batch)') + 
  theme(axis.text.x = element_text(size = 12)) 
p6a <- plot_peptide_corr_distribution(combat_corrected_matrix_platerow,  peptide_annotation,
                                      protein_col = 'Accession',  plot_title = 'Combat corrected (by plate row)') + 
  theme(axis.text.x = element_text(size = 12)) 
grid.arrange(p4a, p5a, p6a, ncol = 3)



# peptide correlation on ratio data
p4 <- plot_peptide_corr_distribution(quantile_normalized_mtx_ratio,  peptide_annotation,
                                                  protein_col = 'Accession',  plot_title = 'Quantile normalized') + 
  theme(axis.text.x = element_text(size = 12)) 

p5 <- plot_peptide_corr_distribution(combat_corrected_mtx_batch_ratio,  peptide_annotation,
                               protein_col = 'Accession',  plot_title = 'Combat corrected (by batch)') + 
  theme(axis.text.x = element_text(size = 12)) 
p6 <- plot_peptide_corr_distribution(combat_corrected_mtx_platerow_ratio,  peptide_annotation,
                                     protein_col = 'Accession',  plot_title = 'Combat corrected (by plate row)') + 
  theme(axis.text.x = element_text(size = 12)) 
grid.arrange(p4, p5, p6, ncol = 3)



