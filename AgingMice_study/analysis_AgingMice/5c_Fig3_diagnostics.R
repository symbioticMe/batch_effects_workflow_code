library(proBatch)
library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(png)
library(here)
source('../lib/plot_themes.R')

#if "AgingMice_study" is a project, the following command will set wd in it:
library(here)
#check if the path is "/my-path-to-code/batch_effects_workflow_code/AgingMice_study"
here()

#load the data:
sample_annotation_AgingMice = read_csv("data_AgingMice/3_data_for_plots/sample_annotation_AgingMice.csv")
normalized_df_AgingMice = read_csv(file = "data_AgingMice/3_data_for_plots/normalized_proteome_AgingMice.csv")
color_list_AgingMice = readRDS("data_AgingMice/3_data_for_plots/color_annotation.rda")

normalized_matrix_AgingMice = long_to_matrix(normalized_df_AgingMice, 
                                             feature_id_col = 'peptide_group_label',
                                             measure_col = 'Intensity', 
                                             sample_id_col = 'FullRunName')


#Panel A: PCA on MS batch
pca_MS_batch_AgingMice <- plot_PCA(normalized_matrix_AgingMice, 
                                   sample_annotation_AgingMice, 
                                   fill_the_missing = -1,
                                   color_by = 'MS_batch', 
                                   plot_title = "PCA colored by MS batch",
                                   color_scheme = color_list_AgingMice[["MS_batch"]], base_size = 16)

pca_EarTag_AgingMice <- plot_PCA(normalized_matrix_AgingMice, 
                                 sample_annotation_AgingMice %>% 
                                   mutate(EarTag = merge_rare_levels(EarTag, rare_thr = 6)),
                                 fill_the_missing = -1,
                                 color_by = 'EarTag',
                                 color_scheme = color_list_AgingMice[["EarTag"]],
                                 plot_title = "PCA colored by EarTag", base_size = 16)
pca_plot = ggarrange(pca_MS_batch_AgingMice + 
                       guides(colour = guide_legend(override.aes = list(size=5))) + 
                       ggtitle(NULL)+
                       theme(plot.margin = margin(t=.5, l=.5, 
                                                  r=.75, b=.5, unit = "cm"),
                             axis.text = element_text(face = 'bold')),
                     pca_EarTag_AgingMice + 
                       guides(colour = guide_legend(override.aes = list(size=5))) +
                       theme(plot.margin = margin(t=.5, l=.75, 
                                                  r=.75, b=.5, unit = "cm"),
                             axis.text = element_text(face = 'bold'))+ 
                       ggtitle(NULL),
                     ncol = 2, nrow = 1, align = 'h')


#Panel B: hierarchical
factors_to_plot_AgingMice <- c("MS_batch",  "digestion_batch", 'DateTime', 'Diet')
plot_hierarchical_clustering(normalized_matrix_AgingMice, 
                             sample_annotation = sample_annotation_AgingMice,
                             color_list = color_list_AgingMice,
                             factors_to_plot = factors_to_plot_AgingMice,
                             fill_the_missing = NULL,
                             distance = "euclidean", agglomeration = "complete",
                             label_samples = FALSE,
#                             filename = "plots_AgingMice/png_for_plots/4b_hierarchical_normalized_selected_factors_noMissing.png",
                             width=7.135/1.5, height=4.74/1.5, units = 'in')
hier_image <- readPNG("plots_AgingMice/png_for_plots/4b_hierarchical_normalized_selected_factors_noMissing.png")

p_hierarchical <- ggplot() +
  background_image(hier_image) +
  # This ensures that the image leaves some space at the edges
  theme(plot.margin = margin(t=1, l=1, r=1, b=1, unit = "cm"))

#Panel C: PVCA
#for calculation of PVCA, which requires a high-performance computer, see script 
# 4a_PVCA_analysis.R
pvca_normalized_50 <- read_csv("plots_AgingMice/interim_data_for_plots/pvca_normalized_50.csv")
pvca_norm_50 <- plot_PVCA.df(pvca_normalized_50,
                             colors_for_bars = NULL,
                             plot_title = 'Principal Variance Compoment Analysis (PVCA), peptides with >50% complete measurements',
                             theme = 'classic', base_size = 18)


#Panel D: iRT peptide examples
iRTs_illustrative = c("7994_LFLQFGAQGSPFLK_2", "1146_ADVTPADFSEWSK_3")
p_iRT_1 = plot_single_feature(iRTs_illustrative[1], df_long =  normalized_df_AgingMice, 
                                            sample_annotation = sample_annotation_AgingMice,
                                            batch_col = "MS_batch",
                              facet_col = 'peptide_group_label',
                            plot_title = 'Representative peptides with batch-specific bias') +
  ylab('Intensity\n(log2 scale)')
p_iRT_2 = plot_single_feature(iRTs_illustrative[2], df_long =  normalized_df_AgingMice %>%
                                mutate(Step = 'normalized'), 
                              sample_annotation = sample_annotation_AgingMice,
                              batch_col = "MS_batch", 
                              facet_col = 'peptide_group_label') +
  ylab('Intensity\n(log2 scale)')
p_iRT = ggarrange(p_iRT_1+ 
                    rremove('xlab')+
                    theme(plot.title = element_text(face = "bold",
                                                           size = 16, hjust = 0.5)),
                  p_iRT_2, 
                  ncol = 1, nrow = 2, align = 'v')


#combine together:
panel_A_C = ggarrange(pca_plot + 
                        ggtitle('Principal Components #1 and #2')+
                        theme(plot.title = element_text(face = "bold",size = 18, hjust = 0.5))+
                        theme(plot.margin = margin(t=1, l=0, 
                                                   r=0, b=1, unit = "cm")), 
                      pvca_norm_50 + theme_publication_mild()+
                        theme(axis.text.y = element_text(size = 16))+
                        theme(axis.text.x = element_text(size = 15))+
                        theme(legend.position=c(0.6,0.7))+
                        theme(plot.margin = margin(t=1, l=1.5, 
                                                   r=1, b=1, unit = "cm")), 
                      ncol = 1,
                      nrow = 2,font.label = list(size=24), 
                      heights = c(0.8, 1), labels = c('A','C'))
panel_B_D = ggarrange(p_hierarchical + 
                        ggtitle('Hierarchical clustering') +
                        theme(plot.title = element_text(face = "bold",size = 18, 
                                                        hjust = 0.6)), 
                      p_iRT+
                        theme(plot.margin = margin(t=1, l=0, 
                                                   r=1.2, b=1.5, unit = "cm")), 
                      ncol = 1,
                      nrow = 2,font.label = list(size=24), 
                      heights = c(0.8, 1), align = "v", labels = c('B','D'))

fig_2 = ggarrange(panel_A_C +
                    theme(plot.margin = margin(r=1, unit = "cm")), 
                  panel_B_D, ncol = 2, nrow = 1, widths = c(2, 1))

ggsave(fig_2, filename = 'plots_AgingMice/Fig3_diagnostics_v2.pdf',
       width = 20, height = 11)

