############### proBatch analysis on Fabian Wendt & Sandra Goetze DDA data #################################
# biological (animal, time, condition) and technical (plate_row, plate_col, ms_order (given as YYMM_SEQ))
# conducted quantile normalization + combat (by batch / plate row)
# Chloe H. Lee, 04/24/2021
############################################################################################################

library(tidyverse)
library(proBatch)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(here)
here() #if you set up the project correctly in "bariatric_surg_mice" this will type the correct path
source("funs/wes_palettes.R")


#################### load data ##################################
# Note: I relabelled batch as robot_batch for simplicity of initial analysis
original_mtx <- read.csv('data_bariatric/20210413_Progenesis_grouped-conditions_peptide-output.csv', skip = 2)
quantData <- read.csv("data_bariatric/pep_mtx_raw.csv") %>% dplyr::rename(peptide_group_label = FullPeptideName)
sample_annotation <- read.csv("data_bariatric/meta.csv") %>% dplyr::rename(FullRunName = X, robot_batch = batch) %>%
  mutate(order = as.numeric(order)) %>% mutate(time = as.factor(time))
peptide_annotation <- original_mtx %>% dplyr::select(Sequence, Modifications, Accession) %>% 
  mutate(peptide_group_label = paste0(Sequence, Modifications))
length(intersect(peptide_annotation$peptide_group_label, quantData$peptide_group_label)) #all match corresponding to quantData 

# data matrix
data_matrix <- quantData %>% column_to_rownames(var = 'peptide_group_label') %>% as.matrix()

# color annotation table 
color_list <- sample_annotation_to_colors(sample_annotation,
                                          factor_columns = c('animal', 'time', 'plate_row', 'plate_col', 'robot_batch', 'condition'),
                                          numeric_columns = c('order'))



#################### illustrating experimental design ######################################
plate_row_design <- sample_annotation %>% 
  ggplot(aes_string(x = 'order', y = 'animal')) + 
  geom_point(aes(colour = plate_row), size = 4) + theme_bw()+
  scale_color_manual(values = color_list$plate_row)

robot_batch_design <- sample_annotation %>% 
  ggplot(aes_string(x = 'order', y = 'time')) + 
  geom_jitter(aes(colour = robot_batch), size = 3,height = 0.1) + theme_bw()+
  scale_color_manual(values = color_list$robot_batch)

plate_col_dsgn <- sample_annotation %>% 
  ggplot(aes_string(x = 'order', y = 'condition')) + 
  geom_jitter(aes(colour = plate_col), size = 3,height = 0.1)+
  scale_color_manual(values = color_list$plate_col) + theme_bw()

gg_design <- ggarrange(plate_row_design, 
                       ggarrange(robot_batch_design, plate_col_dsgn, nrow = 2), 
                       ncol = 2)


###################### Initial assessment ##################################
gg_mean <- plot_sample_mean(data_matrix, sample_annotation, order_col = 'order',
                 batch_col = 'robot_batch', color_by_batch = TRUE, ylimits = c(3, 7),
                 color_scheme = color_list[['robot_batch']], vline_color = NULL)

log_transformed_long <- matrix_to_long(data_matrix)
batch_col = 'robot_batch'
plot_boxplot(log_transformed_long, sample_annotation, batch_col = batch_col, color_scheme = color_list[[batch_col]]) +
  theme(axis.text.x = element_text(size = 8))



################## Normalization and assessment of normalized matrix ############################
quantile_normalized_matrix = normalize_data_dm(data_matrix, normalize_func = 'quantile')
quantile_normalized_long <- matrix_to_long(quantile_normalized_matrix) %>% 
  mutate_if(is.factor, as.character)

gg_mean_norm <- plot_sample_mean(quantile_normalized_matrix, sample_annotation,  
                                 color_by_batch = TRUE, ylimits = c(5, 6), batch_col = 'robot_batch',
                 color_scheme = color_list[['robot_batch']], vline_color = NULL)
gg_boxplot_norm <- plot_boxplot(quantile_normalized_long, sample_annotation, batch_col = batch_col, color_scheme = color_list[[batch_col]]) +
  theme(axis.text.x = element_text(size = 8))



################ Diagnostics on normalized matrix ###################################################
# hierarchical clustering 
selected_annotations <- c('animal', 'time', 'condition', 'plate_row', "plate_col" , 'order', 'robot_batch')

fill_missing <- round(min(data_matrix, na.rm = T))-.5

# heatmap 
plot_heatmap_diagnostic(quantile_normalized_matrix, sample_annotation,
                        factors_to_plot = selected_annotations,
                        cluster_cols = TRUE,
                        color_list = color_list,
                        show_rownames = FALSE, show_colnames = FALSE)

# pvca 
pvca_joint <- plot_PVCA(quantile_normalized_matrix, sample_annotation,
                        technical_factors = c('robot_batch', 'plate_row'),   
                        biological_factors = c('condition', 'time'),
                        fill_the_missing = fill_missing, plot_title = 'quantile normalized')

remove <- FALSE
if (!remove) {
  plot_hierarchical_clustering(quantile_normalized_matrix,
                               sample_annotation = sample_annotation,
                               color_list = color_list,
                               factors_to_plot = selected_annotations,
                               distance = 'euclidean', agglomeration = 'complete',
                               label_samples = FALSE)
  
  #PCA
  gg_PCA_row <- plot_PCA(quantile_normalized_matrix, sample_annotation, 
                         color_by = 'plate_row')
  #I can't reproduce the PCA plot!!!
  
  gg_PCA_batch <- plot_PCA(quantile_normalized_matrix, sample_annotation, 
                           color_by = 'robot_batch') #that's some clustering, but not that "perfect"
  
}


# Missing_value_analysis --------------------------------------------------
###### analyze peptides with missing values 
length(which(rowSums(is.na(quantile_normalized_matrix)) > 0)) #469/1987 peptides contain missing values 

no_na_per_peptides <- rowSums(is.na(quantile_normalized_matrix)) 
hist(no_na_per_peptides/ncol(quantile_normalized_matrix), main = '% of missing values per peptide')
hist(no_na_per_peptides[which(no_na_per_peptides != 0)]/ncol(quantile_normalized_matrix), main = '% of missing values (w/o 0) per peptide')

hist(no_na_per_peptides, main = '# of missing values per peptide')
hist(no_na_per_peptides[which(no_na_per_peptides != 0)], main = '# of missing values (w/o 0) per peptide')


# Filter out too missing --------------------------------------------------
# data with peptides > 10 missing values removed 
quantile_normalized_matrix_filterna <- quantile_normalized_matrix[rowSums(is.na(quantile_normalized_matrix)) <= 10,] 
quantile_normalized_long_fltrd <- matrix_to_long(quantile_normalized_matrix_filterna) %>%  mutate_if(is.factor, as.character)


# batch correction: median centering  -------------------------------------
###### Correction by ComBat
# by robot batch 
combat_corrected_lg_batch <- correct_with_ComBat_df(quantile_normalized_long_fltrd, sample_annotation,  batch_col = "robot_batch")
combat_corrected_matrix_batch <- long_to_matrix(combat_corrected_lg_batch,   measure_col = 'Intensity')

# by robot batch 
medcenter_corrected_lg_batch <- center_feature_batch_medians_df(quantile_normalized_long_fltrd, sample_annotation,
                                                  batch_col = "robot_batch", keep_all = 'default')
medcenter_corrected_matrix_batch <- long_to_matrix(medcenter_corrected_lg_batch,     measure_col = 'Intensity')

# by plate row 
if (!remove) {
  
  # data with peptides > 10 missing values removed and replace NA with -5.5 (suggested by Patrick)
  quantile_normalized_matrix_filterna <- quantile_normalized_matrix[rowSums(is.na(quantile_normalized_matrix)) <= 10,] 
  quantile_normalized_matrix_filterna[is.na(quantile_normalized_matrix_filterna)] = -5.5
  quantile_normalized_long_for_combat <- matrix_to_long(quantile_normalized_matrix_filterna) %>%  mutate_if(is.factor, as.character)
  
  ###### Correction by ComBat
  # by robot batch 
  combat_corrected_lg_batch <- correct_with_ComBat_df(quantile_normalized_long_for_combat, sample_annotation,  batch_col = "robot_batch")
  combat_corrected_matrix_batch <- long_to_matrix(combat_corrected_lg_batch,   measure_col = 'Intensity')
  
  # by plate row
  combat_corrected_lg_platerow <- correct_with_ComBat_df(quantile_normalized_long_for_combat, sample_annotation,  batch_col = "plate_row")
  combat_corrected_matrix_platerow <- long_to_matrix(combat_corrected_lg_platerow,   measure_col = 'Intensity')
  
  medcenter_corrected_lg_platerow <- center_feature_batch_medians_df(quantile_normalized_long_for_medcenter, sample_annotation,
                                                                     batch_col = "plate_row", keep_all = 'default')
  medcenter_corrected_matrix_platerow <- long_to_matrix(medcenter_corrected_lg_platerow,     measure_col = 'Intensity')
}




# negative control:PVCA --------------------------------------------------------
#PVCA
pvca_joint_med <- plot_PVCA(medcenter_corrected_matrix_batch, sample_annotation,
                            technical_factors = c('robot_batch', 'plate_row'),   
                            biological_factors = c('condition', 'time'),
                            fill_the_missing = fill_missing, plot_title = 'after median centering')
pvca_joint_ComBat <- plot_PVCA(combat_corrected_matrix_batch, sample_annotation,
                            technical_factors = c('robot_batch', 'plate_row'),   
                            biological_factors = c('condition', 'time'),
                            fill_the_missing = fill_missing, plot_title = 'after ComBat')


pvca_panels <- ggarrange(pvca_joint, pvca_joint_ComBat, nrow = 2, labels = c('A', 'B'))


# QC: sample correlation --------------------------------------------------
# sample correlation - having robot_batch as batch col
p1.corr <- calculate_sample_corr_distr(quantile_normalized_matrix,   sample_annotation,
                                       batch_col = 'robot_batch',   biospecimen_id_col = 'animal')
p2.corr <- calculate_sample_corr_distr(combat_corrected_matrix_batch,   sample_annotation,
                                       batch_col = 'robot_batch',   biospecimen_id_col = 'animal')


  #plot batch_replicate
  p1.corrplot <- plot_sample_corr_distribution.corrDF(p1.corr, plot_title = 'normalized data')
  
  p2.corrplot <- plot_sample_corr_distribution.corrDF(p2.corr, plot_title = 'after ComBat')
  value_df1 <- p1.corr %>% group_by(batch_replicate) %>% summarise(median = median(correlation), n = n())
  
  value_df2 <- p2.corr %>% group_by(batch_replicate) %>% summarise(median = median(correlation), n = n())
  
  p1.corrplot_f <- p1.corrplot+geom_label(data = value_df1, 
                                          aes(y = median, x = batch_replicate, label = round(median, 3)), 
                                          color = "brown", fill = "white", alpha = .3, fontface = "bold", 
                                          size = 7, vjust = .05)+
    coord_flip()
  p2.corrplot_f <- p2.corrplot+geom_label(data = value_df2, 
                                          aes(y = median, x = batch_replicate, label = round(median, 3)), 
                                          color = "brown", fill = "white", alpha = .3, fontface = "bold", 
                                          size = 7, vjust = .05)+
    coord_flip()
  corrplot_batch <- ggarrange(p1.corrplot_f, p2.corrplot_f, nrow = 2)

#plot by row plate
  p1.corr <- calculate_sample_corr_distr(quantile_normalized_matrix,   sample_annotation,
                                         batch_col = 'plate_row',   biospecimen_id_col = 'animal')
  p2.corr <- calculate_sample_corr_distr(combat_corrected_matrix_batch,   sample_annotation,
                                         batch_col = 'plate_row',   biospecimen_id_col = 'animal')
  
  p1.corrplot <- plot_sample_corr_distribution.corrDF(p1.corr, plot_title = 'normalized data')
  
  p2.corrplot <- plot_sample_corr_distribution.corrDF(p2.corr, plot_title = 'after ComBat')
  value_df1 <- p1.corr %>% group_by(batch_replicate) %>% summarise(median = median(correlation), n = n())
  
  value_df2 <- p2.corr %>% group_by(batch_replicate) %>% summarise(median = median(correlation), n = n())
  
  p1.corrplot_f <- p1.corrplot+geom_label(data = value_df1, 
                                          aes(y = median, x = batch_replicate, label = round(median, 3)), 
                                          color = "brown", fill = "white", alpha = .3, fontface = "bold", 
                                          size = 7, vjust = .05)+
    coord_flip()
  p2.corrplot_f <- p2.corrplot+geom_label(data = value_df2, 
                                          aes(y = median, x = batch_replicate, label = round(median, 3)), 
                                          color = "brown", fill = "white", alpha = .3, fontface = "bold", 
                                          size = 7, vjust = .05)+
    coord_flip()
  corrplot_plateRow <- ggarrange(p1.corrplot_f, p2.corrplot_f, nrow = 2) %>%
    annotate_figure(., top = text_grob("correlation of sample (batch - plate row)"))

  
  correlations_original <- ggarrange(corrplot_batch %>%
                                     annotate_figure(., top = text_grob("correlation of sample (batch - robot batch)")), 
                                     corrplot_plateRow %>%
                                     annotate_figure(., top = text_grob("correlation of sample (batch - plate row)")), ncol = 2)
  
  
#################### prepare data normalized to time 0 #####################################
time0_animals <- sample_annotation %>%   mutate(time = as.numeric(as.character(time))) %>% filter(time == 0) %>% pull(animal) %>% unique

quantile_normalized_lg_ratio <- quantile_normalized_long_fltrd %>% 
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

# sample correlation - having robot_batch as batch col
p1.corr <- calculate_sample_corr_distr(quantile_normalized_mtx_ratio,   sample_annotation,
                                       batch_col = 'robot_batch',   biospecimen_id_col = 'animal')
p2.corr <- calculate_sample_corr_distr(combat_corrected_mtx_batch_ratio,   sample_annotation,
                                       batch_col = 'robot_batch',   biospecimen_id_col = 'animal')
#plot batch_replicate
p1.corrplot <- plot_sample_corr_distribution.corrDF(p1.corr, plot_title = 'quantile normalized')

p2.corrplot <- plot_sample_corr_distribution.corrDF(p2.corr, plot_title = 'after Combat')
value_df1 <- p1.corr %>% group_by(batch_replicate) %>% summarise(median = median(correlation), n = n())

value_df2 <- p2.corr %>% group_by(batch_replicate) %>% summarise(median = median(correlation), n = n())

p1.corrplot_f <- p1.corrplot+geom_label(data = value_df1, 
                                        aes(y = median, x = batch_replicate, label = round(median, 3)), 
                                        color = "brown", fill = "white", alpha = .3, fontface = "bold", 
                                        size = 7, vjust = .05)+
  coord_flip()
p2.corrplot_f <- p2.corrplot+geom_label(data = value_df2, 
                                        aes(y = median, x = batch_replicate, label = round(median, 3)), 
                                        color = "brown", fill = "white", alpha = .3, fontface = "bold", 
                                        size = 7, vjust = .05)+
  coord_flip()
corrplot_ratios_batch <- ggarrange(p1.corrplot_f, p2.corrplot_f, 
                                   nrow = 2, labels = c('C', 'D'))


# Final figure assembly! --------------------------------------------------

variant_1 <- ggarrange(pvca_panels, corrplot_ratios_batch, ncol =2, widths = c(1.3, 1))
ggsave(variant_1, file = 'output_bariatric/bariatric_figure.pdf', 
       device = cairo_pdf(), width = 14, height = 8)                                                                        

# sample correlation - having robot_batch as batch col
p1.corr <- calculate_sample_corr_distr(quantile_normalized_mtx_ratio,   sample_annotation,
                                       batch_col = 'plate_row',   biospecimen_id_col = 'animal')
p2.corr <- calculate_sample_corr_distr(combat_corrected_mtx_batch_ratio,   sample_annotation,
                                       batch_col = 'plate_row',   biospecimen_id_col = 'animal')
#plot batch_replicate
p1.corrplot <- plot_sample_corr_distribution.corrDF(p1.corr, plot_title = 'quantile normalized')

p2.corrplot <- plot_sample_corr_distribution.corrDF(p2.corr, plot_title = 'after ComBat')
value_df1 <- p1.corr %>% group_by(batch_replicate) %>% summarise(median = median(correlation), n = n())

value_df2 <- p2.corr %>% group_by(batch_replicate) %>% summarise(median = median(correlation), n = n())

p1.corrplot_f <- p1.corrplot+geom_label(data = value_df1, 
                                        aes(y = median, x = batch_replicate, label = round(median, 3)), 
                                        color = "brown", fill = "white", alpha = .3, fontface = "bold", 
                                        size = 7, vjust = .05)+
  coord_flip()
p2.corrplot_f <- p2.corrplot+geom_label(data = value_df2, 
                                        aes(y = median, x = batch_replicate, label = round(median, 3)), 
                                        color = "brown", fill = "white", alpha = .3, fontface = "bold", 
                                        size = 7, vjust = .05)+
  coord_flip()
corrplot_ratios_plateRow <- ggarrange(p1.corrplot_f, p2.corrplot_f, nrow = 2)

correlations_ratios <- ggarrange(corrplot_ratios_batch %>%
                                   annotate_figure(., top = text_grob("correlation of sample (batch - robot batch)")), corrplot_ratios_plateRow %>%
                                   annotate_figure(., top = text_grob("correlation of sample (batch - plate row)")), ncol = 2)


if (!remove){
  combat_corrected_lg_platerow_ratio <- combat_corrected_lg_platerow %>% 
    left_join(sample_annotation, by = 'FullRunName') %>%
    mutate(time = as.numeric(as.character(time))) %>%
    filter(animal %in% time0_animals) %>%
    group_by(animal, peptide_group_label) %>%
    mutate(dIntensity = Intensity - Intensity[time == 0]) %>% 
    filter(time != 0)
  combat_corrected_mtx_platerow_ratio <- long_to_matrix(combat_corrected_lg_platerow_ratio,     measure_col = 'dIntensity')
}






###################### evaluate batch corrected matrices ######################################################
# PVCA 
if (!remove) {
  #this is not the way PVCA should be approached!!!
  pvca_analysis_for_dda_data(combat_corrected_matrix_batch)
  pvca_analysis_for_dda_data(medcenter_corrected_matrix_batch)
  pvca_analysis_for_dda_data(combat_corrected_matrix_platerow)
  pvca_analysis_for_dda_data(medcenter_corrected_matrix_platerow)
}


# pheatmap on ratio data 
if (!remove) {
  selected_annotations <- c('animal', 'time', 'condition', 'plate_row', "plate_col" , 'order', 'robot_batch')
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

}




######################### sample and peptide correlation ###################
p3 <- plot_sample_corr_distribution(combat_corrected_mtx_platerow_ratio,   sample_annotation,
                                    batch_col = 'robot_batch',   biospecimen_id_col = 'animal',
                                    plot_title = 'ComBat corrected (by plate row)',    plot_param = 'batch_replicate') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12)) + ylim(-0.3, 1)  + xlab(NULL)
grid.arrange(p1, p2, p3, ncol = 3)



if (F) {
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
}




