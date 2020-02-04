library(readr)
library(dplyr)
library(proBatch)

raw_transitome_PanCancer = read_csv('data_PanCancer/1_original_data/feature_alignment_modB_trs.csv')


raw_transitome_PanCancer[raw_transitome_PanCancer$Intensity %in% 0,"Intensity"]<-NA

raw_transitome_PanCancer = raw_transitome_PanCancer %>% 
  mutate(transition = paste(FragmentIon, ProductCharge, sep = '_'),
         sample_name = paste(Condition, Run, sep = '_'))

#raw_transitome_PanCancer = log_transform_df(raw_transitome_PanCancer,
#                                          log_base = 2, offset = 1)

write_csv(raw_transitome_PanCancer,
          'data_PanCancer/2_interim_data/raw_transitome_PanCancer.csv')