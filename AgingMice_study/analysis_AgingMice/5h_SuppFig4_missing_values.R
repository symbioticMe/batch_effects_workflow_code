library(readr)
library(proBatch)
library(dplyr)

library(ggplot2)
library(ggpubr)


#Haao peptides
Haao_peptides = c('13503_VLEQGQHR_2', 
                  '3353_ELQAGTSLSLFGDSYETQVIAHGQGSSK_3',
                  '4754_AQGSVALSVTQDPAR_3')

normalized_df_withRequants = read_csv(file = "data_AgingMice/3_data_for_plots/normalized_proteome_AgingMice_with_requants.csv")
sample_annotation_AgingMice = read_csv("data_AgingMice/3_data_for_plots/sample_annotation_AgingMice.csv")
requant_peptides = plot_single_feature(Haao_peptides, 
                                       df_long =  normalized_df_withRequants, 
                              sample_annotation = sample_annotation_AgingMice,
                              batch_col = "MS_batch",qual_col = 'm_score', qual_value = 2, 
                              facet_col = 'peptide_group_label',
                              plot_title = 'Representative peptides with batch-specific requant bias')


#PVCA with varying degree of missingness:
#TODO: remove legends
#TODO: align bottom axes
pvca_normalized_100 <- read_csv("plots_AgingMice/interim_data_for_plots/pvca_naOmit_normalized.csv")
pvca_norm_100 <- plot_PVCA.df(pvca_normalized_100,
                             colors_for_bars = NULL,
                             plot_title = 'PVCA, > 100% complete',
                             theme = 'classic')

pvca_normalized_70 <- read_csv("plots_AgingMice/interim_data_for_plots/pvca_normalized_70.csv")
pvca_norm_70 <- plot_PVCA.df(pvca_normalized_70,
                             colors_for_bars = NULL,
                             plot_title = 'PVCA, > 70% complete',
                             theme = 'classic')

pvca_normalized_50 <- read_csv("plots_AgingMice/interim_data_for_plots/pvca_normalized_50.csv")
pvca_norm_50 <- plot_PVCA.df(pvca_normalized_50,
                             colors_for_bars = NULL,
                             plot_title = 'PVCA, > 50% complete',
                             theme = 'classic')

pvca_normalized_90 <- read_csv("plots_AgingMice/interim_data_for_plots/pvca_normalized_90.csv")
pvca_norm_90 <- plot_PVCA.df(pvca_normalized_90,
                             colors_for_bars = NULL,
                             plot_title = 'PVCA, > 90% complete',
                             theme = 'classic')

panel_bottom = ggarrange(pvca_norm_100 + ylim(c(0, .35))+
                         theme(plot.margin = margin(t = 1, l =.5, r=.5,
                                                    b=2, unit = "cm")), 
                         pvca_norm_90 + ylim(c(0, .35)), 
                         pvca_norm_70 + ylim(c(0, .35)),
                       ncol = 3, nrow = 1, 
                       labels = c('A','B','C'),
                       font.label = list(size=22))

gg_FigSupplX = ggarrange(requant_peptides, 
                      panel_bottom +
                        theme(plot.margin = margin(b = 1.2, unit = 'cm')), 
                      ncol = 1, nrow = 2,
                      heights = c(1,1))
ggsave(gg_Fig4_2, 
       filename = 'plots_AgingMice/Supp_Fig4_missing_values.pdf', 
       dev = cairo_pdf, width = 10.95, height = 7.92)
