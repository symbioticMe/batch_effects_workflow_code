library(proBatch)
library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(png)
source('../lib/plot_themes.R')
source('../lib/helpers.R')

#load the data:
sample_annotation_AgingMice = read_csv("data_AgingMice/3_data_for_plots/sample_annotation_AgingMice.csv")
proteome_df_AgingMice = read_csv("data_AgingMice/2_interim_data/raw_proteome_AgingMice.csv")

normalized_df_AgingMice = read_csv(file = "data_AgingMice/3_data_for_plots/normalized_proteome_AgingMice.csv")
color_list_AgingMice = readRDS("data_AgingMice/3_data_for_plots/color_annotation.rda")

#Panel A: correlation of samples in raw data
proteome_matrix_AgingMice = long_to_matrix(proteome_df_AgingMice, 
                                             feature_id_col = 'peptide_group_label',
                                             measure_col = 'Intensity', 
                                             sample_id_col = 'FullRunName')
corr_heatmap = plot_sample_corr_heatmap(data_matrix = proteome_matrix_AgingMice, 
                                        sample_annotation = sample_annotation_AgingMice, 
                                        sample_id_col = 'FullRunName', 
                                        factors_to_plot = c('MS_batch','digestion_batch'), 
                                        color_list = color_list_AgingMice,
                                        plot_title = 'Sample correlation (raw data)',
                                        show_colnames = F, show_rownames = F,
                                        annotation_legend = F)
gg_corr_heatmap = ggarrange(corr_heatmap$gtable)


#Panel B: PCA colored by digestion batch and Diet
normalized_matrix_AgingMice = long_to_matrix(normalized_df_AgingMice, 
                                             feature_id_col = 'peptide_group_label',
                                             measure_col = 'Intensity', 
                                             sample_id_col = 'FullRunName')


##Panel B left: PCA on digestion batch
pca_digesitonBatch_batch_AgingMice <- plot_PCA(normalized_matrix_AgingMice, 
                                   sample_annotation_AgingMice, 
                                   fill_the_missing = -1,
                                   color_by = 'digestion_batch', 
                                   plot_title = "PCA colored by digestion batch",
                                   color_scheme = color_list_AgingMice[["digestion_batch"]])
pca_digesitonBatch_batch_AgingMice1 <- addSmallLegend(pca_digesitonBatch_batch_AgingMice, 
                                                     pointSize = 1.5, textSize = 10, spaceLegend = .2)

##Panel B right: PCA on Diet
pca_Diet_AgingMice <- plot_PCA(normalized_matrix_AgingMice, 
                                   sample_annotation_AgingMice, 
                                   fill_the_missing = -1,
                                   color_by = 'Diet', 
                                   plot_title = "PCA colored by Diet",
                                   color_scheme = color_list_AgingMice[["Diet"]])
pca_Diet_AgingMice1 <- addSmallLegend(pca_Diet_AgingMice, 
                                                      pointSize = 1.5, textSize = 10, spaceLegend = .2)

pca_plot = ggarrange(pca_digesitonBatch_batch_AgingMice1 + 
                       ggtitle(NULL)+
                       theme_publication_mild(),
                     pca_Diet_AgingMice1+ 
                       ggtitle(NULL)+ theme_publication_mild(),
                     ncol = 2, nrow = 1, align = 'h', widths = c(1.1, 1))

panel_A_B = ggarrange(gg_corr_heatmap+
                        theme(plot.margin = margin(t=1, l=0.5, 
                                                   r=2.25, b=1, unit = "cm")), 
  pca_plot + 
                        ggtitle('Principal Components #1 and #2')+
                        theme(plot.title = element_text(face = "bold",size = 14, hjust = 0.5))+
                        theme(plot.margin = margin(t=1, l=0.25, 
                                                   r=0.5, b=2.5, unit = "cm")),
                      ncol = 2,
                      nrow = 1,font.label = list(size=24), 
                      widths = c(0.45, 1), labels = c('A','B'))

#Panel C: normalization reduces bias
iRTs_illustrative = c("7994_LFLQFGAQGSPFLK_2", "1146_ADVTPADFSEWSK_3")
peptides_illustrative = c('1384_AVDSLVPIGR_2', '13610_VLSIGDGIAR_2')
p_iRT_1_raw = plot_single_feature(peptides_illustrative[1], df_long =  proteome_df_AgingMice, 
                              sample_annotation = sample_annotation_AgingMice,
                              batch_col = "MS_batch",
                              facet_col = 'peptide_group_label',
                              plot_title = 'Raw data')
iRT1_limits = ggplot_build(p_iRT_1_raw)$layout$panel_params[[1]]$y.range
p_iRT_1_raw = p_iRT_1_raw  +
  ylim(iRT1_limits)  +
  ylab('Intensity\n(log2 scale)')
p_iRT_2_raw = plot_single_feature(peptides_illustrative[2], df_long =  proteome_df_AgingMice %>%
                                mutate(Step = 'normalized'), 
                              sample_annotation = sample_annotation_AgingMice,
                              batch_col = "MS_batch", 
                              facet_col = 'peptide_group_label')
iRT2_limits = ggplot_build(p_iRT_2_raw)$layout$panel_params[[1]]$y.range
p_iRT_2_raw = p_iRT_2_raw +
  ylim(iRT2_limits) +
  ylab('Intensity\n(log2 scale)')
p_iRT_raw = ggarrange(p_iRT_1_raw+ 
                    rremove('xlab')+
                    theme(plot.title = element_text(face = "bold",
                                                    size = 16, hjust = 0.5)),
                  p_iRT_2_raw, 
                  ncol = 1, nrow = 2, align = 'v')

p_iRT_1 = plot_single_feature(peptides_illustrative[1], df_long =  normalized_df_AgingMice, 
                              sample_annotation = sample_annotation_AgingMice,
                              batch_col = "MS_batch",
                              facet_col = 'peptide_group_label',
                              plot_title = 'Normalized')+
  ylim(iRT1_limits) +
  theme(axis.title.y=element_blank(),
        axis.text.y = element_blank())
p_iRT_2 = plot_single_feature(peptides_illustrative[2], df_long =  normalized_df_AgingMice %>%
                                mutate(Step = 'normalized'), 
                              sample_annotation = sample_annotation_AgingMice,
                              batch_col = "MS_batch", 
                              facet_col = 'peptide_group_label')+
  ylim(iRT2_limits)+
  theme(axis.title.y=element_blank(),
        axis.text.y = element_blank())
p_iRT = ggarrange(p_iRT_1+ 
                    rremove('xlab')+
                    theme(plot.title = element_text(face = "bold",
                                                    size = 16, hjust = 0.5)),
                  p_iRT_2, 
                  ncol = 1, nrow = 2, align = 'v')


raw_norm_comparison = ggarrange(p_iRT_raw,
                  p_iRT, 
                  ncol = 2, nrow = 1, align = 'v', widths = c(2.35, 2))

#CV (not included in the plot)
CV_raw = proteome_df_AgingMice %>% 
  group_by(peptide_group_label) %>% 
  summarize(n = sum(!is.na(Intensity)), sd = sd(Intensity, na.rm = T), 
            mean = mean(Intensity, na.rm = T)) %>% 
  mutate(CV = sd/mean, Step = 'raw')
CV_norm = normalized_df_AgingMice %>% 
  group_by(peptide_group_label) %>% 
  summarize(n = sum(!is.na(Intensity)), sd = sd(Intensity, na.rm = T), 
            mean = mean(Intensity, na.rm = T)) %>% 
  mutate(CV = sd/mean, Step = 'normalized')
CV_df = rbind(CV_norm, CV_raw) %>% filter(n > 300)
summary(CV_norm$CV)
summary(CV_raw$CV)
CV_df$Step = factor(CV_df$Step, levels = c('raw', 'normalized'))
gg_CV = ggplot(CV_df, aes(x = Step, y = CV))+
  geom_boxplot()+
  scale_y_log10()+
  theme_bw()

#Panel D: span demonstration
loess_fit_75_df = read_csv(file = "data_AgingMice/3_data_for_plots/adjusted_fit_df_agingMice.csv")
loess_fit_25_df = read_csv(file = "data_AgingMice/3_data_for_plots/adjusted_fit_25_agingMice_no_requants.csv")

gg_fit_raw1 = plot_with_fitting_curve(peptides_illustrative[1], 
                                 df_long = loess_fit_25_df,
                                 sample_annotation = sample_annotation_AgingMice,
                                 fit_df = loess_fit_25_df,
                                 measure_col = 'preTrendFit_Intensity',
                                 color_by_batch = T, 
                                 facet_col = 'peptide_group_label', 
                                 color_scheme = color_list_AgingMice[['MS_batch']],
                                 plot_title = NULL) +
  ylab('Intensity\n(log2 scale)')
gg_fit_raw2 = plot_with_fitting_curve(peptides_illustrative[2], 
                                      df_long = loess_fit_25_df,
                                      sample_annotation = sample_annotation_AgingMice,
                                      fit_df = loess_fit_25_df,
                                      measure_col = 'preTrendFit_Intensity',
                                      color_by_batch = T, 
                                      facet_col = 'peptide_group_label', 
                                      color_scheme = color_list_AgingMice[['MS_batch']],
                                      plot_title = NULL) +
  ylab('Intensity\n(log2 scale)')
p_fit_25 = ggarrange(gg_fit_raw1+ 
                       rremove('legend')+ 
                       rremove('xlab')+
                        ggtitle('Span = 0.25 (normalized data)')+
                        theme(plot.title = element_text(face = "bold",
                                                        size = 16, hjust = 0.5)),
                      gg_fit_raw2 +
                       rremove('legend'), 
                      ncol = 1, nrow = 2, align = 'v')
gg_fit_norm1 = plot_with_fitting_curve(peptides_illustrative[1], 
                                 df_long = loess_fit_75_df,
                                 sample_annotation = sample_annotation_AgingMice,
                                 fit_df = loess_fit_75_df,
                                 measure_col = 'preTrendFit_Intensity',
                                 color_by_batch = T, 
                                 facet_col = 'peptide_group_label', 
                                 color_scheme = color_list_AgingMice[['MS_batch']],
                                 plot_title = NULL)+
  theme(axis.title.y=element_blank(),
        axis.text.y = element_blank())
gg_fit_norm2 = plot_with_fitting_curve(peptides_illustrative[2], 
                                       df_long = loess_fit_75_df,
                                       sample_annotation = sample_annotation_AgingMice,
                                       fit_df = loess_fit_75_df,
                                       measure_col = 'preTrendFit_Intensity',
                                       color_by_batch = T, 
                                       facet_col = 'peptide_group_label', 
                                       color_scheme = color_list_AgingMice[['MS_batch']],
                                       plot_title = NULL)+
  theme(axis.title.y=element_blank(),
        axis.text.y = element_blank())
p_fit_75 = ggarrange(gg_fit_norm1+ 
                       rremove('legend')+ 
                       rremove('xlab')+
                       ylab('Intensity')+
                       ggtitle('Span = 0.75 (normalized data)')+
                       theme(plot.title = element_text(face = "bold",
                                                       size = 16, hjust = 0.5)),
                     gg_fit_norm2+ 
                       ylab('Intensity')+
                       rremove('legend'), 
                     ncol = 1, nrow = 2, align = 'v')
p_span = ggarrange(p_fit_25+ 
                    theme(plot.title = element_text(face = "bold",
                                                    size = 16, hjust = 0.5)),
                   p_fit_75, 
                  ncol = 2, nrow = 1, align = 'v', widths = c(2.5, 2))
panel_bottom2 = ggarrange(raw_norm_comparison+
                            theme(plot.margin = margin(t=1.25, l=0.85, 
                                                       r=1.5, b=.5, unit = "cm")), 
                         p_span+
                           theme(plot.margin = margin(t=1.25, l=0.75, 
                                                      r=.5, b=.5, unit = "cm")), ncol = 2, align = 'v', font.label = list(size=24), 
                         widths = c(2, 2), labels = c('C','D'))

supp_fig3 = ggarrange(panel_A_B, panel_bottom2, nrow =2, heights = c(2.3, 2.5))
ggsave(supp_fig3, filename = 'plots_AgingMice/SuppFig3_pipeline_extras_v1.pdf',
       width = 18, height = 13)

