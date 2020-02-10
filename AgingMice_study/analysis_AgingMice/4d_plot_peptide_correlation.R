library(readr)
library(ggplot2)

source('plot_split_violinplot.R')
source('../lib/plot_themes.R')
peptide_corr_raw <- read_csv('plots_AgingMice/interim_data_for_plots/peptide_cor_raw.csv')
print('file reading completed!')
Sys.time()
peptide_corr_batchCorr <- read_csv('plots_AgingMice/interim_data_for_plots/peptide_cor_corrected.csv')
print('file reading completed!')
Sys.time()

peptide_corr_raw$Step = 'Raw'
peptide_corr_batchCorr$Step = 'Corrected'

peptide_corr = rbind(peptide_corr_batchCorr, peptide_corr_raw)
peptide_corr$Step = factor(peptide_corr$Step, levels = c('Raw', 'Corrected'))

peptide_corr_plot = plot_split_violin_with_boxplot(peptide_corr, 
                                                   y_col = 'correlation', 
                                                   col_for_color = 'Step', 
                                                   col_for_box = 'same_protein', 
                                                   hlineintercept = 0) +
  xlab(NULL)

ggsave(peptide_corr_plot + theme_publication(),
       filename = 'plots_AgingMice/Fig6D_peptide_correlation_beforeAfter.pdf',
       device = cairo_pdf, width = 5.61, height = 4.47, units = 'in')

ggsave(peptide_corr_plot + theme_publication(),
       filename = 'plots_AgingMice/Fig6D_peptide_correlation_beforeAfter.eps',
       device = "eps", width = 5.61, height = 4.47, units = 'in')


corr_plot = ggplot_gtable(ggplot_build(peptide_corr_plot))
saveRDS(corr_plot, 'plots_AgingMice/interim_data_for_plots/peptide_corrplot.rds')


