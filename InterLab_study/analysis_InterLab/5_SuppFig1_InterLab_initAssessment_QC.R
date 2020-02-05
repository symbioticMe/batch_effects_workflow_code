library(readr)
library(proBatch3.4)
library(dplyr)
library(ggpubr)

protein_df_raw <- read_csv("data_InterLab/3_final_data_InterLab/protein_df_raw.csv")
sample_annotation_InterLab <- read_csv("data_InterLab/3_final_data_InterLab/sample_annotation_InterLab.csv")

protein_df_raw = log_transform_df(protein_df_raw, log_base = 10, measure_col = 'response')
protein_df_raw = protein_df_raw%>% 
  mutate(filename_new = gsub('\\.mzXML\\.gz', '', run_id))
raw_proteome_matrix = protein_df_raw  %>% 
  long_to_matrix(., feature_id_col = 'protein_id', measure_col = 'response', 
                 sample_id_col = 'filename_new')

#extract proteins that are measured in at least half of the samples
n_missing = rowSums(is.na(raw_proteome_matrix))
good_proteins = names(n_missing)[n_missing < 50]
raw_proteome_matrix = raw_proteome_matrix[good_proteins,]

sample_corr_heatmap = plot_sample_corr_heatmap(raw_proteome_matrix, 
                                               sample_annotation = sample_annotation_InterLab, 
                                               sample_id_col = 'filename_new', 
                                               factors_to_plot = 'site', 
                                               plot_title = 'Sample correlation (raw data)',
                                               show_colnames = F, show_rownames = F,
                                               annotation_legend = F)

protein_df_medianCentered <- read_csv("data_InterLab/3_data_for_plots/protein_df_medianCentered.csv")
protein_df_medianCentered = log_transform_df(protein_df_medianCentered, log_base = 10, measure_col = 'response')
protein_df_medianCentered = protein_df_medianCentered %>% 
  mutate(filename_new = gsub('\\.mzXML\\.gz', '', run_id))
norm_proteome_matrix = protein_df_medianCentered %>% 
  long_to_matrix(., feature_id_col = 'protein_id', measure_col = 'response', 
                 sample_id_col = 'filename_new')

norm_proteome_matrix = norm_proteome_matrix[good_proteins,]

sample_corr_heatmap_normalized = plot_sample_corr_heatmap(norm_proteome_matrix, 
                                               sample_annotation = sample_annotation_InterLab, 
                                               sample_id_col = 'filename_new', 
                                               factors_to_plot = 'site', 
                                               plot_title = 'Sample correlation (normalized data)',
                                               show_colnames = F, show_rownames = F)
#ToDo: extract, which part here comes from helpers?
source('lib/helpers.R')
boxplot_proteins_normalized = plot_boxplot(df_long = protein_df_medianCentered, 
                                           sample_annotation = sample_annotation_InterLab,
                                           sample_id_col = 'filename_new',
                                           batch_col = 'site', 
                                           order_col = NULL,
                                           measure_col =  'response', 
                                           color_by_batch = T,
                                           color_scheme = colors_for_sites$site,
                                           facet_col = NULL, outliers = F) +
  ylab('Log10 protein abundance') +
  theme(axis.text.x = element_text(size = rel(0.25)), 
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12), 
        axis.title.x=element_blank(),
        legend.position = 'none')

boxplot_proteins_raw = plot_boxplot(df_long = protein_df_raw, 
                                           sample_annotation = sample_annotation_InterLab,
                                           sample_id_col = 'filename_new',
                                           batch_col = 'site', 
                                           order_col = NULL,
                                           measure_col =  'response', 
                                           color_by_batch = T, color_scheme = colors_for_sites$site,
                                           facet_col = NULL, outliers = F) +
  ylab('Log10 protein abundance') +
  theme_cute_legend()+ 
  guides(fill = guide_legend(nrow = 1)) + 
  theme(axis.text.x = element_text(size = rel(0.25)), 
        axis.title.x=element_blank())


#===== Plot the selected spike-in peptides before and after normalization ====
spike_in_quantities = read_csv('data_InterLab/3_data_for_plots/spike_in_quantities.csv')

selected_spikeIns = c('APAELEVECATQLR','SGGLLQLWK')
gg_spike_ins = ggplot(spike_in_quantities %>% filter(peptide_sequence %in% selected_spikeIns), 
                aes(x = concentration, y = Intensity, 
                    group = site, color = site))+
  geom_line()+
  facet_grid(peptide_id ~ Step, scales = 'free')+
  theme_calc()+
  scale_y_log10(labels = scales::comma, minor_breaks = c(10^(log10(5) + 4:5), 10^(log10(2.5) + 4:6)))+
  scale_x_log10()+
  scale_color_manual(values = colors_for_sites$site) + 
  theme(panel.border=element_rect(colour="black",size=1, fill = NA),
        strip.text.x = element_text(size = rel(1.5), face = 'bold'), 
        strip.background.x = element_blank(), 
        panel.grid.minor.x = element_line(colour = 'lightgrey', size = .1),
        panel.spacing = unit(1, "lines"),
        axis.title = element_text(size = 12, color = 'black'))+
  xlab('Concentration (fmol on column)') +
  ylab('Peptide peak area')

#===== Calculate  and plot the CV of proteins in raw and normalized data====
protein_df_raw$Step = 'Raw'
protein_df_medianCentered$Step = 'Normalized'
summary_proteome = rbind(protein_df_raw, protein_df_medianCentered)
summary_proteome = summary_proteome  %>% filter(protein_id %in% solid_proteins) %>% group_by(Step, protein_id) %>% 
  summarize(CV = 100*sd(beforeLog_response, na.rm = T)/mean(beforeLog_response, na.rm = T), 
            mean = mean(beforeLog_response, na.rm = T))
max_intensity = max(protein_df_medianCentered$response, protein_df_raw$response)
summary_proteome = summary_proteome  %>% 
  mutate(mean_range = cut(log10(mean + 1), 
                          breaks=c(0, seq(4.25, 6.5, by = .25), 
                                   max_intensity)))

gg_CV = ggplot(summary_proteome, aes(x = mean_range, y = CV, fill = Step))+
  geom_boxplot(width = 0.52, outlier.size = 0.1)+
  theme_classic()+
  scale_fill_manual(values = c('red', 'grey'))+
  scale_y_log10(breaks = c(12.5, 25, 50, 100, 200, 400))+
  theme(legend.title = element_blank()) +
  ylab('CV (coefficient of variation)')+
  scale_x_discrete(labels= c('<4.25','4.25-4.5','4.5-4.75','4.75-5','5-5.25','5.25-5.5','5.5-5.75','5.75 - 6','6 - 6.25','6.25 - 6.5', '>6.5')) +
  xlab('Protein intensity')+
  theme(legend.position=c(.9581,.9581), legend.justification = c(1,1), 
        legend.key.size = unit(.52, "in"), 
        legend.text = element_text(size=10),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 10))

##==== Assemble the panels with plots ============================
###===== Assemble the panel A: boxplot of the raw data matrix ====
panel_boxplot_raw = ggarrange(boxplot_proteins_raw  +
                                theme(plot.margin = margin(l=1, t = 1, r = 1,
                                                           b = 1,
                                                           unit = "cm"),
                                      axis.title.y = element_text(size = 14),
                                      axis.text.y = element_text(size = 12)), 
                              labels= 'A', 
                              font.label = list(size=22))

###===== Assemble the B and C panels (sample correlation heatmap and spike-ins)
panel_middle = ggarrange(ggarrange(sample_corr_heatmap$gtable)+
                           theme(plot.margin = margin(l=30, r = 20)),
                         gg_spike_ins+ theme(legend.position = 'none', 
                                             plot.background = element_rect(color = NA)), 
                         ncol = 2, nrow = 1, labels = c('B', 'C'),
                         font.label = list(size=22))

###===== Assemble the panel D: boxplot of the normalized data matrix ====
panel_boxplot_norm = ggarrange(boxplot_proteins_normalized  +
                                theme(plot.margin = margin(l=1, t = 1, r = 1,
                                                           b = 1,
                                                           unit = "cm"),
                                      axis.title.y = element_text(size = 14),
                                      axis.text.y = element_text(size = 12)), 
                              labels= 'D', 
                              font.label = list(size=22))

###===== Assemble the panel E: CV distribution comparison, binned by protein abundance ====
panel_CV = ggarrange(gg_CV + xlab('Log10 protein abundance') +
                       theme(plot.margin = margin(l=1, t = 1, r = 1,
                                                  b = 1,
                                                  unit = "cm")),
                     labels= 'E', 
                     font.label = list(size=22))

###==== Assemble the figure together: panels A, B&C, D and E =====
fig_assembled = ggarrange(panel_boxplot_raw, 
                          panel_middle,
                          panel_boxplot_norm,
                          panel_CV, ncol = 1, nrow = 4)
ggsave(fig_assembled, filename = 'plots_InterLab/Supp_InterLab.pdf',
       width = 12.52, height = 22.44, units = 'in',device = cairo_pdf)
ggsave(fig_assembled, filename = 'plots_InterLab/Supp_InterLab.png',
       width = 12.52, height = 22.44, units = 'in', device = png)
