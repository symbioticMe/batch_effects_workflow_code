library(proBatch)
library(readr)
library(dplyr)
library(MSstats)
library(tibble)
library(ggplot2)
library(ggdendro)
library(ggpubr)
library(cowplot)
# library(dendextend)


source('plot_dendrogram_replicates.R')

sample_annotation_PanCancer = read_csv('data_PanCancer/3_data_for_plots/sample_annotation_PanCancer.csv') %>% as.data.frame()
proteome_raw = read_csv('data_PanCancer/3_data_for_plots/proteome_df_raw.csv')
proteome_normalized = read_csv('data_PanCancer/3_data_for_plots/proteome_df_normalized.csv')
proteome_corrected = read_csv('data_PanCancer/3_data_for_plots/proteome_df_corrected.csv')

#transform "long" data frames into matrices:
matrix_normalized = proteome_normalized %>% 
  long_to_matrix(feature_id_col = 'Protein', 
                 measure_col = 'Intensity', 
                 sample_id_col = 'sample_name')

matrix_batchCorrected = proteome_corrected %>% 
  long_to_matrix(feature_id_col = 'Protein', 
                 measure_col = 'Intensity', 
                 sample_id_col = 'sample_name')

#define the colors
color_list_PanCancer = sample_annotation_to_colors(sample_annotation = sample_annotation_PanCancer,
                                                   sample_id_col = 'sample_name',
                                                   numeric_columns = NULL,
                                                   factor_columns = c('Batch', 'tissue', 'case_control', 'patient_ID'))


###==========1. boxplots (panels A and B)===========
boxlot_PanCancer_raw = plot_boxplot(proteome_raw, 
                                     sample_annotation_PanCancer, 
                                     sample_id_col = 'sample_name', measure_col = "Intensity",
                                     order_col = NULL,color_by_batch = T, batch_col = 'Batch',
                                     color_scheme = color_list_PanCancer[['Batch']])+
  theme(axis.text.x = element_text(size = 3), 
        axis.title.x = element_blank(),legend.position = 'none')

boxlot_PanCancer_norm = plot_boxplot(proteome_normalized, 
                                     sample_annotation_PanCancer, 
                                     sample_id_col = 'sample_name', measure_col = "Intensity",
                                     order_col = NULL,color_by_batch = T, batch_col = 'Batch',
                                     color_scheme = color_list_PanCancer[['Batch']])+
  theme(axis.text.x = element_text(size = 3), 
        axis.title.x = element_blank(),legend.position = 'none')

panel_boxplots = ggarrange(boxlot_PanCancer_raw, boxlot_PanCancer_norm, labels = c('A', 'B'))



###=========2. Protein-level Spike-in plotting (panels C and D)================
fetubin_norm = plot_single_feature(feature_name = "1/sp|Q58D62|FETUB_BOVIN", 
                                      df_long = proteome_normalized,
                                      sample_annotation = sample_annotation_PanCancer,
                                      order_col = NULL, sample_id_col = "sample_name",
                                      batch_col = "Batch", measure_col = "Intensity", 
                                      feature_id_col = 'Protein', geom = "point",
                                      color_by_batch = T, 
                                      color_scheme = color_list_PanCancer[['Batch']], 
#                                      facet_col = 'Step', 
                                      plot_title = NULL,
                                      vline_color = NULL, theme = 'classic')+
  theme(strip.background = element_rect(fill = "white", colour = NA), 
        strip.text.x = element_text(face = 'bold', size = 16),
        axis.text.x = element_text(size = 3), 
        axis.title.x = element_blank(),
        legend.position = "none")
y_lims_fetubin = layer_scales(fetubin_norm)$y$range$range

fetubin_corr = plot_single_feature(feature_name = "1/sp|Q58D62|FETUB_BOVIN", 
                                           df_long = proteome_corrected,
                                           sample_annotation = sample_annotation_PanCancer,
                                           order_col = NULL, sample_id_col = "sample_name",
                                           batch_col = "Batch", measure_col = "Intensity", 
                                           feature_id_col = 'Protein', geom = "point",
                                           color_by_batch = T, 
                                           color_scheme = color_list_PanCancer[['Batch']], 
                                           plot_title = NULL,
                                           vline_color = NULL, theme = 'classic')+
  theme(strip.background = element_rect(fill = "white", colour = NA), 
        strip.text.x = element_text(face = 'bold', size = 16),
        axis.text.x = element_text(size = 3), 
        axis.title.x = element_blank(),
        legend.position = "none")+
  ylim(y_lims_fetubin)

panel_proteins = ggarrange(fetubin_norm, fetubin_corr, 
                           labels = c('C','D'))

###=========clustering (better after protein quantification, though)===============
factors_to_plot = c("Batch", "tissue", "case_control")
color_for_heatmap = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)

heatmap_norm <- plot_heatmap_diagnostic(matrix_normalized, 
                                         sample_annotation = sample_annotation_PanCancer,
                                         sample_id_col = "sample_name", 
                                         factors_to_plot = factors_to_plot,
                                         fill_the_missing = T, 
                                         cluster_rows = T, cluster_cols = T, 
                                         color_list = color_list_PanCancer,
                                         heatmap_color = color_for_heatmap, 
                                         color_for_missing = "black", 
                                         #             filename = 'PanCancer/graphs_PanCancer/6_PanCancer_clustering_heatmap_after_MSstats.png',
                                         plot_title = 'Before batch correction',
                                         show_rownames = F, show_colnames = F, width = 7.49, height = 6.14)



heatmap_corrected <- plot_heatmap_diagnostic(matrix_batchCorrected,
                                  sample_annotation = sample_annotation_PanCancer,
                                  sample_id_col = "sample_name", 
                                  factors_to_plot = factors_to_plot,
                                  fill_the_missing = T, 
                                  cluster_rows = T, cluster_cols = T, 
                                  color_list = color_list_PanCancer,
                                  heatmap_color = color_for_heatmap, 
                                  color_for_missing = "black", 
                                  #             filename = 'PanCancer/graphs_PanCancer/6_PanCancer_clustering_heatmap_after_MSstats.png',
                                  plot_title = 'After batch correction',
                                  show_rownames = F, show_colnames = F, width = 7.49, height = 6.14)



##======Replicate clustering======================
sample_annotation_PanCancer = sample_annotation_PanCancer %>% 
  mutate(patient_IDs = merge_rare_levels(patient_ID))
names(color_list_PanCancer)[4] = 'patient_IDs'

source('plot_dendrogram_replicates.R')

gg_norm = plot_dendrogram_replicates(matrix_normalized, 
                          sample_annotation_PanCancer, 
                          color_list_PanCancer, plot_title = NULL)

gg_corrected = plot_dendrogram_replicates(matrix_batchCorrected, 
                               sample_annotation_PanCancer, 
                               color_list_PanCancer, plot_title = NULL)


#===== Join all the plots together in one plot:

#to bring it to the state of ggplot graphs:

panel_heatmaps = ggarrange(heatmap_norm$gtable %>%
                      gtable_remove_grobs(., c('annotation_legend', 'legend'), trim = FALSE), 
                      heatmap_corrected$gtable, ncol = 2, labels = c('E','F'))
panel_replicates = ggarrange(gg_norm, gg_corrected, labels = c('G','H'))

figure_supp_Tatjana_with_boxplots = ggarrange(panel_boxplots,
                                         panel_proteins, 
                                         panel_heatmaps,
                                         panel_replicates, nrow = 4, ncol = 1, 
                                         heights = c(0.95, 0.75, 2, 1))

ggsave(figure_supp_Tatjana_with_boxplots, filename = 'plots_PanCancer/Supp_Fig2_PanCancer.png',
       width = 15, height = 18, units = 'in', device = 'png')
ggsave(figure_supp_Tatjana_with_boxplots, filename = 'plots_PanCancer/Supp_Fig2_PanCancer.pdf',
       width = 15, height = 18, units = 'in', device = cairo_pdf)
