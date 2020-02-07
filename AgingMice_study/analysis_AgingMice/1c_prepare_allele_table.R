library(dplyr)
library(readr)
library(reshape2)

Protein_Gene_allelle_mapping <- read.csv("data_AgingMice/1_original_data/Protein_Gene_LOCATION_AND_ALLELE_SYMBOL_v2.csv")
peptide_df_AgingMice = read_csv("data_AgingMice/1_original_data/peptide_annotation.csv")
sample_annotation_AgingMice = read_csv("data_AgingMice/3_data_for_plots/sample_annotation_AgingMice.csv")


rownames(Protein_Gene_allelle_mapping) = Protein_Gene_allelle_mapping$Gene
allelle_annotation_wide = Protein_Gene_allelle_mapping[unique(peptide_df_AgingMice$Gene), 
                                                                     c('Gene', setdiff(unique(sample_annotation_AgingMice$Strain), 'Mix'))]
allelle_annotation_df = melt(data = allelle_annotation_wide, id.vars = 'Gene', measure.vars = names(allelle_annotation_wide)[-1])
names(allelle_annotation_df)[2] = 'Strain'
allelle_annotation_df$allele = ifelse(allelle_annotation_df$value == 'red', 
                                      'DBA2J', 
                                      ifelse(allelle_annotation_df$value == 'blue','C57BL6J', 'heterozygous'))


write_csv(allelle_annotation_df, 'data_AgingMice/3_data_for_plots/allelle_annotation_df.csv')
