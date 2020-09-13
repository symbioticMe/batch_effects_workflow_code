library(readr)
library(dplyr)
library(proBatch)
library(reshape2)
library(ggplot2)
library(ggpubr)

sample_annotation_AgingMice = read_csv("data_AgingMice/3_data_for_plots/sample_annotation_AgingMice.csv")
color_list_AgingMice = readRDS("data_AgingMice/3_data_for_plots/color_annotation.rda")
loess_fit_75_df = read_csv(file = "data_AgingMice/3_data_for_plots/adjusted_fit_df_agingMice.csv")
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
y_axis_fit = layer_scales(gg_fit)$y$range$range
ggsave(gg_fit, filename = 'plots_AgingMice/Fig5_batch_correction.pdf',
       width = 12.09, height = 4.4, unit = 'in', device = cairo_pdf)
#extract legend
leg = get_legend(gg_fit)


y_axis_fit_up = c(14.5, 17.5)
y_axis_fit_down = c(17, 22.5)
up_peptide = "1146_ADVTPADFSEWSK_3"
down_peptide = "7994_LFLQFGAQGSPFLK_2"
gg_fit_A_up = plot_with_fitting_curve(up_peptide, 
                                   df_long = illustrative_df %>% 
                                     filter(Step == 'Normalized, with LOESS fit'),
                                   sample_annotation = sample_annotation_AgingMice,
                                   fit_df = fit_df  %>% 
                                     filter(Step == 'Normalized, with LOESS fit'),
                                   facet_col = 'Step', 
                                   color_by_batch = T, 
                                   color_scheme = color_list_AgingMice[['MS_batch']],
                                   plot_title = NULL,
                                   base_size = 12)+
  theme(strip.background.x= element_blank(),
        strip.text.x = element_text(size = rel(1.42), face = 'bold'),
        strip.text.y = element_blank() ,
        strip.background.y = element_blank(),
        legend.position = 'none',
        axis.title.x = element_blank())+
  ylim(y_axis_fit_up)
gg_fit_A_down = plot_with_fitting_curve(down_peptide, 
                                      df_long = illustrative_df %>% 
                                        filter(Step == 'Normalized, with LOESS fit'),
                                      sample_annotation = sample_annotation_AgingMice,
                                      fit_df = fit_df  %>% 
                                        filter(Step == 'Normalized, with LOESS fit'),
                                      facet_col = 'Step', 
                                      color_by_batch = T, 
                                      color_scheme = color_list_AgingMice[['MS_batch']],
                                      plot_title = NULL,  
                                      base_size = 12)+
  theme(strip.background.x= element_blank(),
        strip.text.x = element_blank() ,
        legend.position = 'none')+
  ylim(y_axis_fit_down)
gg_fit_A = ggarrange(gg_fit_A_up, gg_fit_A_down, nrow = 2)




gg_fit_B_up = plot_with_fitting_curve(up_peptide, 
                                      df_long = illustrative_df %>% 
                                        filter(Step == 'After LOESS, with medians'),
                                      sample_annotation = sample_annotation_AgingMice,
                                      fit_df = fit_df  %>% 
                                        filter(Step == 'After LOESS, with medians'),
                                      facet_col = 'Step', 
                                      color_by_batch = T, 
                                      color_scheme = color_list_AgingMice[['MS_batch']],
                                      plot_title = NULL,                                     base_size = 12)+
  geom_line(data = peptide_mean_df %>% 
              filter(peptide_group_label %in% up_peptide),
            aes(y = mean, x=order, group = MS_batch, color = MS_batch), size=1.5)+
  theme(strip.background.x= element_blank(),
        strip.text.x = element_text(size = rel(1.42), face = 'bold'),
        strip.text.y = element_blank() ,
        strip.background.y = element_blank(),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  ylim(y_axis_fit_up)
gg_fit_B_down = plot_with_fitting_curve(down_peptide, 
                                        df_long = illustrative_df %>% 
                                          filter(Step == 'After LOESS, with medians'),
                                        sample_annotation = sample_annotation_AgingMice,
                                        fit_df = fit_df  %>% 
                                          filter(Step == 'After LOESS, with medians'),
                                        facet_col = 'Step', 
                                        color_by_batch = T, 
                                        color_scheme = color_list_AgingMice[['MS_batch']],
                                        plot_title = NULL,                                     base_size = 12)+
  geom_line(data = peptide_mean_df %>% 
              filter(peptide_group_label %in% down_peptide),
            aes(y = mean, x=order, group = MS_batch, color = MS_batch), size=1.5)+
  theme(strip.background.x= element_blank(),
        strip.text.x = element_blank() ,
        legend.position = 'none',
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  ylim(y_axis_fit_down)
gg_fit_B = ggarrange(gg_fit_B_up, gg_fit_B_down, nrow = 2)

gg_fit_C_up = plot_with_fitting_curve(up_peptide, 
                                      df_long = illustrative_df %>% 
                                        filter(Step == 'Corrected data'),
                                      sample_annotation = sample_annotation_AgingMice,
                                      fit_df = fit_df  %>% 
                                        filter(Step == 'Corrected data'),
                                      facet_col = 'Step', 
                                      color_by_batch = T, 
                                      color_scheme = color_list_AgingMice[['MS_batch']],
                                      plot_title = NULL,                                     base_size = 12)+
  facet_grid(peptide_group_label~Step)+
  theme(strip.background.x= element_blank(),
        strip.text.x = element_text(size = rel(1.52), face = 'bold'),
        strip.text.y = element_text(size = rel(.93)),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  ylim(y_axis_fit_up)
gg_fit_C_down = plot_with_fitting_curve(down_peptide, 
                                        df_long = illustrative_df %>% 
                                          filter(Step == 'Corrected data'),
                                        sample_annotation = sample_annotation_AgingMice,
                                        fit_df = fit_df  %>% 
                                          filter(Step == 'Corrected data'),
                                        facet_col = 'Step', 
                                        color_by_batch = T, 
                                        color_scheme = color_list_AgingMice[['MS_batch']],
                                        plot_title = NULL, 
                                        base_size = 12)+
  facet_grid(peptide_group_label~Step)+
  theme(strip.background.x= element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = rel(.93)),
        legend.position = 'none',
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  ylim(y_axis_fit_down)


gg_fit_C = ggarrange(gg_fit_C_up, gg_fit_C_down, nrow = 2)

gg_total = ggarrange(
  gg_fit_A+
    theme(plot.margin = margin(t=.31, unit = "cm")), 
  gg_fit_B+
    theme(plot.margin = margin(t=.31, l = .51, unit = "cm")), 
  gg_fit_C +
    theme(plot.margin = margin(t=.31, l = .51, unit = "cm")), 
  labels = LETTERS[1:3], ncol = 3,
  widths = c(11, 10.5, 10.8),font.label = list(size=18)
)
gg_total1 = ggarrange(gg_total, 
                      as_ggplot(leg)+
                        theme(plot.margin = margin(t=1, l=1.5, 
                                                   r=2, b=1, unit = "cm")), 
                      widths = c(8, 1))
ggsave(gg_total1, filename = 'plots_AgingMice/Fig5_batch_correction1.pdf',
       width = 15.09, height = 5.4, unit = 'in', device = cairo_pdf)

