library(MSstats)
library(readr)
library(dplyr)
library(proBatch)

corrected_transitome_proBatch = read_csv('data_PanCancer/2_interim_data/batchCorrected_transitome_PanCancer.csv')

start.time = Sys.time()
corrected_transitome_unlogged = corrected_transitome_proBatch %>% unlog_df()
corrected_data_proBatch_MSstats = dataProcess(corrected_transitome_unlogged, 
                                              logTrans=2, normalization=FALSE, 
                                              #betweenRunInterferenceScore=FALSE, 
                                              fillIncompleteRows=TRUE)


proteome_corrected = corrected_data_proBatch_MSstats$RunlevelData %>% 
  mutate(sample_name = paste(GROUP_ORIGINAL, SUBJECT_ORIGINAL, sep = '_'))%>%
  rename(Intensity = LogIntensities)

write_csv(proteome_corrected, 'data_PanCancer/3_data_for_plots/proteome_df_corrected.csv')

end.time = Sys.time()
print(end.time - start.time)
