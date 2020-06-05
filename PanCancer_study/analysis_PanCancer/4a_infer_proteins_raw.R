library(MSstats)
library(readr)
library(dplyr)
library(proBatch)

#load the data
raw_transitome_unlogged = read_csv('data_PanCancer/2_interim_data/raw_transitome_PanCancer.csv')

##========Protein quantification=============
start.time = Sys.time()
raw_data_proBatch_MSstats = dataProcess(raw_transitome_unlogged, 
                                        logTrans=2, normalization=FALSE, 
                                        #betweenRunInterferenceScore=FALSE, 
                                        fillIncompleteRows=TRUE)
proteome_raw = raw_data_proBatch_MSstats$RunlevelData %>% 
  mutate(sample_name = paste(GROUP_ORIGINAL, SUBJECT_ORIGINAL, sep = '_'))%>%
  rename(Intensity = LogIntensities)
write_csv(proteome_raw, 'data_PanCancer/3_data_for_plots/proteome_df_raw.csv')
end.time = Sys.time()
print(end.time - start.time)

