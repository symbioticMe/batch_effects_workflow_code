#!/usr/bin/env Rscript --vanilla

#load the libraries
library(readr)
library(aLFQ)

#load the data
data.normalized <- import(ms_filenames=  "data_Interlab/1_original_data/all_sites_global_q_0.01_applied_to_local_global.txt",
                          ms_filetype = "openswath", concentration_filename=NA,
                          averageruns=FALSE, sumruns=FALSE, mprophet_cutoff=0.01, 
                          openswath_superimpose_identifications=FALSE, 
                          openswath_replace_run_id=FALSE, openswath_filtertop=FALSE, 
                          openswath_removedecoys=TRUE)


#infer the proteins
peptides.median <- PeptideInference(data.normalized, transition_topx = 5,
                                    transition_strictness = "loose", transition_summary = "sum", 
                                    consensus_proteins = FALSE, consensus_peptides = FALSE, 
                                    consensus_transitions = FALSE)

#save the data frame (LONG or WIDE format? What are the column names?)
write_csv(peptides.median, path = 'data_Interlab/2_interim_data/peptide_df_raw.csv')