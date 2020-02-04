#!/usr/bin/env Rscript --vanilla

#load the libraries
library(readr)
library(aLFQ)

#load the data
data.median <- import(ms_filenames=  "data_InterLab/2_interim_data/all_sites_global_q_001_applied_to_local_global_medianCentered.tsv",
                      ms_filetype = "openswath", 
                      concentration_filename=NA, averageruns=FALSE, sumruns=FALSE, 
                      mprophet_cutoff=0.01,
                      openswath_superimpose_identifications=FALSE, 
                      openswath_replace_run_id=FALSE,
                      openswath_filtertop=FALSE, openswath_removedecoys=TRUE)

#infer the proteins
prots.median <- ProteinInference(data.median, peptide_method = "top", peptide_topx = 3,
                                 peptide_strictness = "loose", 
                                 peptide_summary = "sum", transition_topx = 5,
                                 transition_strictness = "loose", 
                                 transition_summary = "sum", fasta = NA, model = NA,
                                 combine_precursors = FALSE, 
                                 consensus_proteins = FALSE, consensus_peptides = FALSE, 
                                 consensus_transitions = FALSE)

#save the data frame
write_csv(prots.median, path = 'data_InterLab/3_data_for_plots/protein_df_medianCentered.csv')