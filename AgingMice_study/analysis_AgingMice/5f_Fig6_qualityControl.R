library(readr)
library(proBatch)
library(dplyr)

library(ggplot2)
library(ggpubr)
#library(grid)
library(RColorBrewer)
library(grImport)
#library(png)

sample_annotation_AgingMice = read_csv("data_AgingMice/3_data_for_plots/sample_annotation_AgingMice.csv")
color_list_AgingMice = readRDS("data_AgingMice/3_data_for_plots/color_annotation.rda")
peptide_df_AgingMice = read_csv("data_AgingMice/1_original_data/peptide_annotation.csv")

raw_proteome = read_csv("data_AgingMice/2_interim_data/raw_proteome_AgingMice.csv")
batchCorrected_df_AgingMice  = read_csv("data_AgingMice/3_data_for_plots/batchCorrected_proteome_AgingMice.csv")

#Panel A: Corrected QTL:
allelle_annotation_df = read_csv('data_AgingMice/3_data_for_plots/allelle_annotation_df.csv')

QTL_proteins = c('Acads')
Acads_peptide_3 = '9111_LVIAGHLLR_2'

allelle_annotation_df = allelle_annotation_df %>% 
  filter(Gene == 'Acads')

annotation_df = sample_annotation_AgingMice %>%
  merge(allelle_annotation_df, by = 'Strain')%>%
  merge(peptide_df_AgingMice, by = 'Gene') 

QTL_data_batchCorr = batchCorrected_df_AgingMice %>%
  filter(peptide_group_label == Acads_peptide_3) %>%
  merge(annotation_df, 
        by = intersect(names(batchCorrected_df_AgingMice), names(annotation_df)),
        all.x = T)
#df has columns: Gene, Stain, locus (Blue, Red, heterozygous)
colors_for_alleles = c(RColorBrewer::brewer.pal(3, 'Set1')[1:2], 'grey')
names(colors_for_alleles) = c("DBA2J", "C57BL6J", "heterozygous")
plot_title_batchCorr_3 = expression(bold("QTL detection problem in ")~bolditalic("Acads")~bold(" protein"))
best_QTL_batchCorr_3 = plot_single_feature(Acads_peptide_3, 
                                           df_long = QTL_data_batchCorr,
                                           sample_annotation = sample_annotation_AgingMice,
                                           geom = 'point',
                                           vline_color = 'grey', 
                                           plot_title = plot_title_batchCorr_3, base_size = 16) 
best_QTL_batchCorr_3_1 = best_QTL_batchCorr_3 + 
  geom_point(data = QTL_data_batchCorr %>% filter(peptide_group_label == Acads_peptide_3), 
             aes(color = allele, x = order, y = Intensity)) + 
  scale_color_manual(values = colors_for_alleles, breaks = names(colors_for_alleles))

y_lims_QTL = readRDS('plots_AgingMice/interim_data_for_plots/QTL_range.rds')
best_QTL_batchCorr_3_scaledRaw = best_QTL_batchCorr_3_1+
  ylim(c(y_lims_QTL))+ 
  theme(legend.position="top")




#  transform proteome df into matrix 
proteome_log_AgingMice = long_to_matrix(raw_proteome, feature_id_col = 'peptide_group_label',
                                        measure_col = 'Intensity', sample_id_col = 'FullRunName')
batchCorrected_matrix_AgingMice = long_to_matrix(batchCorrected_df_AgingMice, 
                                                 feature_id_col = 'peptide_group_label',
                                                 measure_col = 'Intensity', 
                                                 sample_id_col = 'FullRunName',
                                                 qual_col = NULL)
rm(raw_proteome)
rm(batchCorrected_df_AgingMice)
#Panel B: Sample correlation


sample_cor_batchCor <-calculate_sample_corr_distr(data_matrix = batchCorrected_matrix_AgingMice, 
                                                  sample_annotation = sample_annotation_AgingMice,
                                                  sample_id_col = 'FullRunName',
                                                  batch_col = 'MS_batch')
sample_cor_batchCor  = sample_cor_batchCor  %>%
  mutate(batch_replicate_upd = ifelse(replicate, 'replicates', 
                                    ifelse(batch_the_same, 'same\nbatch', 'different\nbatches')))
sample_cor_batchCor$batch_replicate_upd = factor(sample_cor_batchCor$batch_replicate_upd,
                                                 levels = c('different\nbatches', 'same\nbatch', 'replicates'))
sample_cor_batchCor$Step = 'After correction'

y_lims_corrplot_raw = readRDS(file = 'plots_AgingMice/interim_data_for_plots/y_lims_corrplot_raw.rds')
corr_distr <- plot_sample_corr_distribution.corrDF(sample_cor_batchCor, 
                                                   plot_title = NULL,
                                                   # plot_title = 'Selected sample correlation distribution',
                                                   plot_param = 'batch_replicate_upd',
                                                   base_size = 16)+ 
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size = 1.2)+
  theme(axis.title.x = element_blank(),
        strip.text = element_text(face = 'bold', size = 16),
        strip.background=element_rect(color = NA,  fill=NA),
        panel.spacing = unit(2, "lines"),
        axis.text.x = element_text(angle = 45,hjust=1, size = 10),
        plot.margin = margin(t=.5, l=0, 
                             r=0, b=1, unit = "cm"))+
  ylim(c(y_lims_corrplot_raw[1],1))

panel_top2 = egg::ggarrange(best_QTL_batchCorr_3_scaledRaw+
                         ggtitle(NULL)+
                         theme(plot.margin = margin(t = .21, l =.5, r=1.5,
                                                    b=1.2, unit = "cm")), 
                       corr_distr+
                         theme(plot.margin = margin(t=.21,l=1, r = 1.2, unit = "cm")), 
                      ncol = 2, nrow = 1, widths = c(2, 1), 
                      labels = c('A','B'),
                      #font.label = list(size=22)) #for ggpubr ggarrange
                      label.args = list(gp = grid::gpar(fontface = "bold", fontsize = 22)))

#panel C: peptide correlation heatmap
peptide_count = rowSums(!is.na(proteome_log_AgingMice))
peptides_complete_80 = rownames(proteome_log_AgingMice)[peptide_count > 0.8 * ncol(proteome_log_AgingMice)]
rm(peptide_count)

peptide_df_AgingMice_80 = peptide_df_AgingMice %>%
  filter(peptide_group_label %in% peptides_complete_80)
interesting_genes = c('Haao', 'Dhtkd1','Nnt','Glo1')

raw_80 = proteome_log_AgingMice[peptides_complete_80, ]
protein_corrplot_plot_raw <- plot_protein_corrplot(raw_80,
                                                   protein_name = interesting_genes,
                                                   peptide_annotation = peptide_df_AgingMice_80,
                                                   plot_title = 'Raw data',
                                                   protein_col = 'Gene',
                                                   breaks = seq(-.6, 1, by = 0.016),
                                                   factors_to_plot = 'Gene', 
                                                   show_colnames = F, 
                                                   show_rownames = F, 
                                                   legend = F, annotation_legend = F)

batch_80 = batchCorrected_matrix_AgingMice[peptides_complete_80, ]
protein_corrplot_plot_corrected <- plot_protein_corrplot(batch_80,
                                                         protein_name = interesting_genes,
                                                         peptide_annotation = peptide_df_AgingMice_80,
                                                         plot_title = 'After correction',
protein_col = 'Gene',
breaks = seq(-.6, 1, by = 0.016),
factors_to_plot = 'Gene', show_colnames = F, show_rownames = F)


prot_corr_heatmap <- ggarrange(ggarrange(protein_corrplot_plot_raw$gtable)+
                                 theme(plot.margin = margin(l=30, t = 10)), 
                               ggarrange(protein_corrplot_plot_corrected$gtable)+
                                 theme(plot.margin = margin(l=45, t = 10)), 
                               ncol = 2, nrow = 1, widths = c(3,4)) %>%
  annotate_figure(top = text_grob("Correlation of peptides within/between protein(s)", 
                                  face = "bold", size = 14))

#TODO: re-plot after the calculation on Euler is complete:
peptide_correlation <- readPNG("plots_AgingMice/png_for_plots/6d_peptide_correlation.png")

p_pep_corr <- ggplot() +
  background_image(peptide_correlation) +
  # This ensures that the image leaves some space at the edges
  theme(plot.margin = margin(t=1, l=1, r=1, b=1, unit = "cm"))


panel_bottom = ggarrange(prot_corr_heatmap, p_pep_corr,  
                         ncol = 2, nrow = 1, widths = c(2, 1), 
                         labels = c('C','D'),
                         font.label = list(size=22))

gg_Fig4_2 = ggarrange(panel_top2, 
                    panel_bottom +
                      theme(plot.margin = margin(b = 1.2, unit = 'cm')), 
                    ncol = 1, nrow = 2,
                    heights = c(1,1))
ggsave(gg_Fig4_2, 
       filename = 'plots_AgingMice/Fig6_quality_control1.pdf', 
       dev = cairo_pdf, width = 10.95, height = 7.92)

