library(proBatch)
library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)
library(png)

sample_annotation_AgingMice = read_csv("data_AgingMice/3_data_for_plots/sample_annotation_AgingMice.csv")
normalized_df_withRequants = read_csv(file = "data_AgingMice/3_data_for_plots/normalized_proteome_AgingMice_with_requants.csv")
color_list_AgingMice = readRDS("data_AgingMice/3_data_for_plots/color_annotation.rda")


normalized_matrix_AgingMice = long_to_matrix(normalized_df_withRequants, 
                                             feature_id_col = 'peptide_group_label',
                                             measure_col = 'Intensity', 
                                             sample_id_col = 'FullRunName', 
                                             qual_col = 'm_score', qual_value = 2)

#selected annotation factors:
factors_to_plot_AgingMice <- c("MS_batch",  "digestion_batch", "Diet", 'DateTime')
breaks_for_heatmap = seq(8, 24, length.out = 101)

n_samples = ncol(normalized_matrix_AgingMice)
pep_completeness = apply(normalized_matrix_AgingMice, MARGIN = 1, FUN = function(x) sum(!is.na(x))/n_samples)
peptides_75 = names(pep_completeness)[pep_completeness > .75 & pep_completeness]

heatmap_normalized_AgingMice_75 <- plot_heatmap_diagnostic(normalized_matrix_AgingMice[peptides_75,], 
                                                           sample_annotation = sample_annotation_AgingMice,
                                                           factors_to_plot = factors_to_plot_AgingMice,
                                                           fill_the_missing = -1,
                                                           sample_id_col = 'FullRunName', 
                                                           color_list = color_list_AgingMice,
                                                           distance = "euclidean", method = "complete",
                                                           breaks = c(-1, breaks_for_heatmap),
                                                           show_rownames = F, show_colnames = F, 
                                                           cluster_cols = T,
                                                           legend = F,
                                                           annotation_legend = F,
                                                           plot_title = 'Peptides with >75% complete measurements')
heatmap_normalized_AgingMice_75_p = heatmap_normalized_AgingMice_75$gtable
ggsave(filename = 'plots_AgingMice/png_for_plots/heatmap_75_normalized_clustered.png', 
       plot = heatmap_normalized_AgingMice_75_p, width = 7.18, height = 7.08, units = 'in')
heatmap_img = readPNG("plots_AgingMice/png_for_plots/heatmap_75_normalized_clustered.png")
heatmap_75 <- ggplot() +
  background_image(heatmap_img) +
  # This ensures that the image leaves some space at the edges
  theme(plot.margin = margin(t=.5, b = .5, unit = "cm"))
heatmap_75_gg = ggarrange(heatmap_75, labels = 'A',
                          font.label = list(size=22))

#Panel B: heatmap of sample correlation for selected samples
back_to_back_samples <- read_csv("data_AgingMice/1_original_data/back_to_back_samples.csv")
replicate_filenames = back_to_back_samples %>% pull(FullRunName)
factors_to_show = c("MS_batch", "EarTag")
# set color annotation
breaksList <- seq(0.8, 1, by = 0.01) # color scale of pheatmap 
heatmap_colors = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)-1)

normalized_matrix_withRequants = long_to_matrix(normalized_df_withRequants, 
                                             feature_id_col = 'peptide_group_label',
                                             measure_col = 'Intensity', 
                                             sample_id_col = 'FullRunName')

corr_with_requants = plot_sample_corr_heatmap(normalized_matrix_withRequants,
                                              samples_to_plot = replicate_filenames, 
                                              sample_annotation = sample_annotation_AgingMice %>% 
                                                select(c('FullRunName', factors_to_show)),
                                              sample_id_col = 'FullRunName',
                                              plot_title = 'With requants', 
                                              factors_to_plot = factors_to_show,
                                              color_list = color_list_AgingMice, 
                                              heatmap_color = heatmap_colors, 
                                              breaks = breaksList, show_rownames = F,
                                              cluster_rows=F, cluster_cols=F,
                                              annotation_names_col = TRUE, 
                                              legend = FALSE,
                                              annotation_legend = FALSE, 
                                              show_colnames = F)

corr_no_requants = plot_sample_corr_heatmap(normalized_matrix_AgingMice,
                                            samples_to_plot = replicate_filenames, 
                                            sample_annotation = sample_annotation_AgingMice %>% 
                                              select(c('FullRunName', factors_to_show)),
                                            sample_id_col = 'FullRunName',
                                            plot_title = 'No requants', 
                                            factors_to_plot = factors_to_show,
                                            color_list = color_list_AgingMice, 
                                            heatmap_color = heatmap_colors, 
                                            breaks = breaksList, show_rownames = F, 
                                            cluster_rows=F, cluster_cols=F,
                                            annotation_names_col = TRUE, 
                                            annotation_legend = FALSE, 
                                            show_colnames = F)
corr_heatmap = ggarrange(ggarrange(corr_with_requants$gtable) +
                           theme(plot.margin = margin(l=.5,r=.3, unit = "cm")), 
                         ggarrange(corr_no_requants$gtable) +
                           theme(plot.margin = margin(l=3,r=1, unit = "cm")), 
                         ncol = 2, nrow = 1, widths = c(4.65,7.25))

#Panel C: correlation distribution plots
sample_cor_withRequants <-calculate_sample_corr_distr(data_matrix = normalized_matrix_withRequants, 
                                             sample_annotation = sample_annotation_AgingMice,
                                             repeated_samples = replicate_filenames,
                                             sample_id_col = 'FullRunName',
                                             batch_col = 'MS_batch')
sample_cor_withRequants = sample_cor_withRequants %>% 
  mutate(batch_replicate_upd = ifelse(replicate, 'replicates', 
                                      ifelse(batch_the_same, 'same batch', 'different batches')))
sample_cor_withRequants$batch_replicate_upd = factor(sample_cor_withRequants$batch_replicate_upd,
                                                   levels = c('different batches', 'same batch', 'replicates'))
sample_cor_withRequants$Step = 'With requants'

sample_cor_noRequants <-calculate_sample_corr_distr(data_matrix = normalized_matrix_AgingMice, 
                                                      sample_annotation = sample_annotation_AgingMice,
                                                      repeated_samples = replicate_filenames,
                                                      sample_id_col = 'FullRunName',
                                                      batch_col = 'MS_batch')
sample_cor_noRequants = sample_cor_noRequants %>% 
  mutate(batch_replicate_upd = ifelse(replicate, 'replicates', 
                                      ifelse(batch_the_same, 'same batch', 'different batches')))
sample_cor_noRequants$batch_replicate_upd = factor(sample_cor_noRequants$batch_replicate_upd,
                                            levels = c('different batches', 'same batch', 'replicates'))
sample_cor_noRequants$Step = 'No requants'
sample_cor_selected = rbind(sample_cor_noRequants, sample_cor_withRequants)
sample_cor_selected$Step = factor(sample_cor_selected$Step,
                                  levels = c('With requants', 'No requants'))
corr_distr <- plot_sample_corr_distribution.corrDF(sample_cor_selected, 
                                                   plot_title = NULL,
                                                   # plot_title = 'Selected sample correlation distribution',
                                                   plot_param = 'batch_replicate_upd')+ 
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size = 1.2)+
  theme(axis.title.x = element_blank(),
        strip.text = element_text(face = 'bold', size = 14),
        strip.background=element_rect(color = NA,  fill=NA),
        panel.spacing = unit(9, "lines"),
        axis.text.x = element_text(angle = 45,hjust=1, size = 10),
        plot.margin = margin(t=.5, l=0, 
                             r=0, b=1, unit = "cm"))


#combine together
panel_B_C = ggarrange(corr_heatmap+ 
                        theme(plot.margin = margin(t = .5, 
                                                   r=1, b = .75, unit = "cm")), 
                      corr_distr + 
                        theme(plot.margin = margin(t = .5, 
                                                   r=3.4, b = 1, unit = "cm")), 
                      labels = c('B','C'), nrow = 2, 
                      heights = c(2.75,2.5),
                      font.label = list(size=22))+
  theme(plot.margin = margin(t=0, l=1.5, 
                             r=0, b=0, unit = "cm"))


gg_Fig4 = ggarrange(heatmap_75_gg, 
                    panel_B_C,
          ncol = 2, nrow = 1, widths = c(2.35, 3))
ggsave(gg_Fig4, 
       filename = 'plots_AgingMice/Fig4_missing_values.pdf', 
       dev = cairo_pdf, width = 17.79, height = 8.66)
