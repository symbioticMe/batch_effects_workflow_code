library(MSstats)
library(readr)
library(dplyr)
library(proBatch)

normalized_transitome_PanCancer = read_csv('data_PanCancer/2_interim_data/normalized_transitome_PanCancer.csv')


start.time = Sys.time()
normalized_transitome_unlogged = normalized_transitome_PanCancer %>% unlog_df()
normalized_data_proBatch_MSstats = dataProcess(normalized_transitome_unlogged, 
                                               logTrans=2, normalization=FALSE, 
                                               #betweenRunInterferenceScore=FALSE, 
                                               fillIncompleteRows=TRUE)
proteome_normalized = normalized_data_proBatch_MSstats$RunlevelData %>% 
  mutate(sample_name = paste(GROUP_ORIGINAL, SUBJECT_ORIGINAL, sep = '_'))%>%
  rename(Intensity = LogIntensities)
write_csv(proteome_normalized,  'data_PanCancer/3_data_for_plots/proteome_df_normalized.csv')
end.time = Sys.time()
print(end.time - start.time)
