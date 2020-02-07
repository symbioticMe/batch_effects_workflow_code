library(dplyr)
library(readr)

Protein_Gene_LOCATION_AND_ALLELE_SYMBOL_v2 <- read.csv("~/Desktop/QTL_and_allele_heterozygosity/plot_intensity_w_allele_heterozygosity/Protein_Gene_LOCATION_AND_ALLELE_SYMBOL_v2.csv")
peptide_df_AgingMice = read_csv("Aging_mice/data_Aging_mice/1_original_data_AgingMice/peptide_annotations_fixed.csv")
sample_annotation_AgingMice = read_csv("Aging_mice/data_Aging_mice/2_interim_data_AgingMice//sample_annotation_Aging_mice.csv")


rownames(Protein_Gene_LOCATION_AND_ALLELE_SYMBOL_v2) = Protein_Gene_LOCATION_AND_ALLELE_SYMBOL_v2$Gene
allelle_annotation_wide = Protein_Gene_LOCATION_AND_ALLELE_SYMBOL_v2[unique(peptide_df_AgingMice$Gene), 
                                                                     c('Gene', setdiff(unique(sample_annotation_AgingMice$Strain), 'Mix'))]
allelle_annotation_df = melt(data = allelle_annotation_wide, id.vars = 'Gene', measure.vars = names(allelle_annotation_wide)[-1])
names(allelle_annotation_df)[2] = 'Strain'
allelle_annotation_df$allele = ifelse(allelle_annotation_df$value == 'red', 
                                      'DBA2J', 
                                      ifelse(allelle_annotation_df$value == 'blue','C57BL6J', 'heterozygous'))


write_csv(allelle_annotation_df, 'Aging_mice/data_Aging_mice/2_interim_data_AgingMice/allelle_annotation_df.csv')
