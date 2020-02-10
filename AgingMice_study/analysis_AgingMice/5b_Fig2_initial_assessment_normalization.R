library(tidyverse)
library(proBatch)
library(ggpubr)
library(ggplot2)

sample_annotation_AgingMice = read_csv("data_AgingMice/3_data_for_plots/sample_annotation_AgingMice.csv")
proteome_df_AgingMice = read_csv("data_AgingMice/2_interim_data/raw_proteome_AgingMice.csv")
color_list_AgingMice = readRDS("data_AgingMice/2_interim_data/color_annotation.rda")

batch_col = 'MS_batch'

#source('lib/helpers.R')
#panel A: meanplot
#  transform proteome df into matrix 
proteome_log_AgingMice = long_to_matrix(proteome_df_AgingMice, feature_id_col = 'peptide_group_label',
                                        measure_col = 'Intensity', sample_id_col = 'FullRunName')
meanplot <- plot_sample_mean(proteome_log_AgingMice, sample_annotation_AgingMice, 
                             sample_id_col = 'FullRunName', color_by_batch = T,
                             order_col = 'order',batch_col = batch_col,
                             plot_title = "Meanplot of Raw Data Matrix", 
                             color_scheme = color_list_AgingMice[[batch_col]])
meanplot = meanplot + 
  #theme_publication()+
  geom_point(size = -1, aes_string(fill = batch_col)) + 
  scale_fill_manual(values = color_list_AgingMice[[batch_col]])+
  guides(fill = guide_legend(nrow = 2, 
                             override.aes = list(shape = 15, size = 7)))

meanplot_reps <- plot_sample_mean(proteome_log_AgingMice, 
                                  sample_annotation_AgingMice %>% 
                                    mutate(EarTag = merge_rare_levels(EarTag, 6)), 
                                  sample_id_col = 'FullRunName', color_by_batch = T,
                                  order_col = 'order', batch_col = 'EarTag',
                                  plot_title = "Mean Intensity of each sample, Raw Data", 
                                  color_scheme = color_list_AgingMice[['EarTag']], 
                                  vline_color = NULL)
meanplot_reps_with_line = proBatch:::add_vertical_batch_borders('order', 'FullRunName', 'MS_batch', 
                                                                 'red', 
                                                                 NULL, sample_annotation_AgingMice, 
                                                                 meanplot_reps)

#panel B: one of QTL protein peptides;
Acads_peptide = '9111_LVIAGHLLR_2'
peptide_df_AgingMice = read_csv("data_AgingMice/1_original_data/peptide_annotation.csv")
allelle_annotation_df = read_csv('data_AgingMice/3_data_for_plots/allelle_annotation_df.csv')
allelle_annotation_df = allelle_annotation_df %>%
  filter(Gene == 'Acads')
peptide_df_AgingMice = peptide_df_AgingMice %>%
  filter(Gene == 'Acads')

annotation_df = sample_annotation_AgingMice %>%
  merge(allelle_annotation_df, by = 'Strain')%>%
  merge(peptide_df_AgingMice, by = 'Gene') 

QTL_data1 = proteome_df_AgingMice %>%
  filter(peptide_group_label == Acads_peptide) %>%
  merge(annotation_df %>%
          filter(peptide_group_label == Acads_peptide), 
        by = c("FullRunName", "peptide_group_label"), 
        all.x = T, all.y = T)

colors_for_alleles = c(RColorBrewer::brewer.pal(3, 'Set1')[1:2], 'grey')
names(colors_for_alleles) = c("DBA2J", "C57BL6J", "heterozygous")
best_QTL = plot_single_feature(Acads_peptide, 
                               df_long = QTL_data1,
                               sample_annotation = sample_annotation_AgingMice,
                               geom = 'point',
                               vline_color = 'grey', 
                               plot_title = 'ACADS protein peptide LVIAGHLLR') 
best_QTL = best_QTL + 
  geom_point(data = QTL_data1 %>% filter(peptide_group_label == Acads_peptide), 
             aes(color = allele, x = order, y = Intensity))+ 
  scale_color_manual(values = colors_for_alleles, breaks = names(colors_for_alleles))
best_QTL = best_QTL + theme(legend.position="top")
y_lims_QTL = layer_scales(best_QTL)$y$range$range
saveRDS(y_lims_QTL, 'plots_AgingMice/interim_data_for_plots/QTL_range.rds')
# best_QTL = best_QTL + theme_publication()
#best_QTL = readRDS("plots_AgingMice/interim_ggplot_objects_AgingMice/2b_best_QTL_raw_data.rds")

# panel C: Boxplots before correction;
boxplot_raw_Aging_mice <- plot_boxplot(proteome_df_AgingMice, sample_annotation_AgingMice,
                                       sample_id_col = 'FullRunName', measure_col = 'Intensity',
                                       order_col = 'order', batch_col = batch_col, 
                                       color_by_batch = T, color_scheme = color_list_AgingMice[[batch_col]],
                                       plot_title = "Intensity Distribution of raw data",
                                       outliers = FALSE)
boxplot_raw_Aging_mice = boxplot_raw_Aging_mice +
  geom_boxplot(fatten = 5000, lwd = .0002, outlier.shape = NA)
  # boxplot_raw_Aging_mice = boxplot_raw_Aging_mice  +  theme_publication()
y_lims_boxplot_raw = layer_scales(boxplot_raw_Aging_mice)$y$range$range


#panel D: Boxplots after correction
normalized_df_AgingMice = read_csv(file = "data_AgingMice/3/normalized_proteome_AgingMice.csv")
boxplot_normalized_Aging_mice <- plot_boxplot(normalized_df_AgingMice, sample_annotation_AgingMice,
                                              sample_id_col = 'FullRunName', measure_col = 'Intensity',
                                              order_col = 'order', batch_col = batch_col, 
                                              color_by_batch = T, color_scheme = color_list_AgingMice[[batch_col]],
                                              plot_title = "Intensity Distribution of normalized data",
                                              outliers = FALSE)
boxplot_normalized_Aging_mice = boxplot_normalized_Aging_mice + 
  geom_boxplot(fatten = 5000, lwd = .0002, outlier.shape = NA)

sample_cor_raw  <-calculate_sample_corr_distr(data_matrix = proteome_log_AgingMice, 
                                              sample_annotation = sample_annotation_AgingMice,
                                              sample_id_col = 'FullRunName',
                                              batch_col = 'MS_batch')
sample_cor_raw  = sample_cor_raw  %>% 
  mutate(batch_replicate_upd = ifelse(replicate, 'replicates', 
                                      ifelse(batch_the_same, 'same batch', 'different batches')))
sample_cor_raw$batch_replicate_upd = factor(sample_cor_raw$batch_replicate_upd,
                                            levels = c('different batches', 'same batch', 'replicates'))
plot_sample_corr_raw  <- plot_sample_corr_distribution.corrDF(sample_cor_raw , 
                                                              plot_title = 'Sample correlation',
                                                              plot_param = 'batch_replicate_upd')+ 
  theme(axis.title.x = element_blank())
y_lims_corrplot_raw = layer_scales(plot_sample_corr_raw)$y$range$range
saveRDS(y_lims_corrplot_raw, file = 'plots_AgingMice/interim_data_for_plots/y_lims_corrplot_raw.rds')




fig_1_ext2 = ggarrange(ggarrange(meanplot + 
                                   theme(plot.title = element_text(face = "bold",
                                                                   size = 16, hjust = 0.5)), 
                                 best_QTL,
                                 plot_sample_corr_raw+
                                   theme(axis.title.x = element_text(color = "white")), 
                                 ncol = 3,  widths = c(2, 2, 1),
                                 labels =c("A","B","C")),
                       ggarrange(boxplot_raw_Aging_mice, boxplot_normalized_Aging_mice, 
                                 ncol = 2, 
                                 labels =c("D","E"),
                                 common.legend = TRUE, legend = "right"), 
                       nrow = 2)
ggsave(fig_1_ext2, filename = 'plots_AgingMice/Fig2_initial_assessment.png',
       width = 25, height = 11, units = 'cm')

QTL_title <- expression(bold("QTL detection problem in ")~bolditalic("Acads")~bold(" protein"))
panel_ABC <- ggarrange(meanplot_reps + 
                         ylab('Mean Intensity')+
                         ggtitle('Intensity drift in raw data') +
                                   theme(plot.title = element_text(face = "bold",
                                                                   size = 16, hjust = 0.5))+
                         theme(plot.margin = margin(t=.5, l=.5, 
                                                    r=.5, b=.5, unit = "cm")),
                       plot_sample_corr_raw+
                         theme(axis.title.x = element_text(color = "white"))+
                         theme(plot.title = element_text(face = "bold",size = 16, hjust = 0.5))+
                         theme(plot.margin = margin(t=.5, l=.5, 
                                                    r=.5, b=.5, unit = "cm")), 
                       best_QTL+ggtitle(QTL_title)+
                         theme(plot.title = element_text(face = "bold",size = 16, hjust = 0.5))+
                         theme(plot.margin = margin(t=.5, l=.5, 
                                                    r=.5, b=.5, unit = "cm")), 
                                 ncol = 3,  widths = c(1.5, 1, 1.5),
                                 labels =c("A","B","C"),font.label = list(size=22))
ggsave(panel_ABC, filename = 'plots_AgingMice/Fig2_ABC_initial_assessment.pdf',
       width = 18, height = 5.5, units = 'in')
panel_DE_no_legend <- ggarrange(boxplot_raw_Aging_mice + 
                        rremove('legend') +
                        theme(plot.title = element_text(face = "bold",size = 16, hjust = 0.5))+
                        theme(plot.margin = margin(t=.5, l=.5, 
                                                   r=.5, b=.5, unit = "cm")), 
                      boxplot_normalized_Aging_mice +
                        rremove('legend') +
                        theme(plot.title = element_text(face = "bold",size = 16, hjust = 0.5))+
                        theme(plot.margin = margin(t=.5, l=.5, 
                                                   r=.5, b=.5, unit = "cm")), 
                      ncol = 2, 
                      labels =c("D","E"),font.label = list(size=22),
                      common.legend = TRUE, legend = "none")
ggsave(panel_DE_no_legend, filename = 'plots_AgingMice/Fig2_DE_initial_assessment.pdf',
       width = 18, height = 6.5, units = 'in')
ggsave(panel_DE_no_legend, filename = 'plots_AgingMice/Fig2_DE_initial_assessment.png',
       width = 18, height = 6.5, units = 'in')
