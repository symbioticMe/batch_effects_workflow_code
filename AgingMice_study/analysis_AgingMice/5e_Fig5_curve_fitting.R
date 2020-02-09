library(readr)
library(proBatch)
library(reshape2)
library(ggplot2)

sample_annotation_AgingMice = read_csv("data_AgingMice/3_data_for_plots/sample_annotation_AgingMice.csv")
color_list_AgingMice = readRDS("data_AgingMice/3_data_for_plots/color_annotation.rda")
loess_fit_75_df = read_csv(file = "data_AgingMice/3_data_for_plots/adjusted_fit_75_AgingMice.csv")
batchCorrected_df_AgingMice = read_csv(file = "data_AgingMice/3_data_for_plots/batchCorrected_proteome_AgingMice.csv")


iRTs_illustrative = c("7994_LFLQFGAQGSPFLK_2", "1146_ADVTPADFSEWSK_3")
selected_peptides = c("11772_SVTAFFNWLR_2", "34963_LFLQFGAQGSPFLK_3",
                      "5141_GLVLIAFSQYLQK_2", "5248_GQETSTNPIASIFAWSR_2" )
interesting_peptides = c(iRTs_illustrative, selected_peptides)

loess_fit_75_forPlot = loess_fit_75_df %>% 
  filter(peptide_group_label %in% interesting_peptides) %>%
  select(peptide_group_label, FullRunName, Intensity, preTrendFit_Intensity) %>%
  melt(id.vars = c('peptide_group_label', 'FullRunName'), 
       value.name = "Intensity") %>%
  mutate(Step = ifelse(variable == 'Intensity', 
                       'After LOESS, with medians', 'Normalized, with LOESS fit')) %>%
  select(-variable)

batchCorrected_df_AgingMice = batchCorrected_df_AgingMice %>%
  filter(peptide_group_label %in% interesting_peptides) %>%
  mutate(Step = 'Corrected data') 
batchCorrected_df_AgingMice = batchCorrected_df_AgingMice %>%
  select(peptide_group_label, FullRunName, Intensity, Step)

illustrative_df = rbind(loess_fit_75_forPlot, batchCorrected_df_AgingMice)
illustrative_df = illustrative_df %>%
  mutate(Step = factor(Step, levels = c('Normalized, with LOESS fit',
                                        'After LOESS, with medians', 'Corrected data')))

#add fit
fit_df = loess_fit_75_df %>% 
  filter(peptide_group_label %in% interesting_peptides) %>% 
  mutate(Step = 'Normalized, with LOESS fit',
         Intensity = preTrendFit_Intensity) %>%
  select(peptide_group_label, FullRunName, Intensity, fit, Step)%>%
  mutate(Step = factor(Step, levels = c('Normalized, with LOESS fit',
                                        'After LOESS, with medians', 'Corrected data')))

#add the means for the corrected
peptide_mean_df = loess_fit_75_df %>% 
  filter(peptide_group_label %in% interesting_peptides) %>%
  select(Intensity, peptide_group_label, FullRunName) %>%
  merge(sample_annotation_AgingMice %>%
          select(FullRunName, order, MS_batch), by = 'FullRunName') %>%
  group_by(peptide_group_label, MS_batch) %>%
  mutate(mean = median(Intensity, na.rm = T),
         Step = 'After LOESS, with medians')%>%
  mutate(Step = factor(Step, levels = c('Normalized, with LOESS fit',
                                        'After LOESS, with medians', 'Corrected data')))

gg_fit = plot_with_fitting_curve(iRTs_illustrative, 
                                   df_long = illustrative_df,
                                   sample_annotation = sample_annotation_AgingMice,
                                   fit_df = fit_df,
                                   facet_col = 'Step', 
                                 color_by_batch = T, 
                                 color_scheme = color_list_AgingMice[['MS_batch']],
                                 plot_title = NULL)+
  geom_line(data = peptide_mean_df %>% 
              filter(peptide_group_label %in% iRTs_illustrative),
            aes(y = mean, x=order, group = MS_batch, color = MS_batch), size=1.5) + 
  scale_color_manual(values = color_list_AgingMice[['MS_batch']])+
  annotate("segment", x=-Inf, xend=+Inf, y=-Inf, yend=-Inf)
ggsave(gg_fit, filename = 'plots_AgingMice/Fig3_batch_correction.pdf',
       width = 12.09, height = 4.4, unit = 'in', device = cairo_pdf)

ggsave(gg_fit + theme(legend.position = "none"), 
       filename = 'plots_AgingMice/Fig3_correction_noLegend.pdf',
       width = 12.09, height = 4.4, unit = 'in', device = cairo_pdf)
