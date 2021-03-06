#!/usr/bin/env Rscript --vanilla

#load the libraries
library(readr)
library(aLFQ)


print(sessionInfo())
print(Sys.time())
#load the data
data.raw <- import(ms_filenames= "data_InterLab/1_original_data/all_sites_global_q_0.01_applied_to_local_global.txt", 
                   ms_filetype = "openswath", concentration_filename=NA, averageruns=FALSE, 
                   sumruns=FALSE, mprophet_cutoff=0.01, openswath_superimpose_identifications=FALSE, 
                   openswath_replace_run_id=FALSE, openswath_filtertop=FALSE, openswath_removedecoys=TRUE)
print('Import complete!')
print(Sys.time())

#infer the proteins
prots <- ProteinInference(data.raw, peptide_method = "top", peptide_topx = 3,
                          peptide_strictness = "loose",peptide_summary = "sum", transition_topx = 5,
                          transition_strictness = "loose",transition_summary = "sum", fasta = NA, model = NA,
                          combine_precursors = FALSE, consensus_proteins = FALSE, consensus_peptides = FALSE, 
                          consensus_transitions = FALSE)
print(Sys.time())

write_csv(prots, path = 'data_InterLab/3_data_for_plots/protein_df_raw.csv')
print(Sys.time())