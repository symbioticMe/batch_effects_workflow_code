This study was used as an example to define the workflow, so all steps are executed
exactly as they were in the original study, with minor polishing of the code to
`proBatch` format and phylosophy.

In this study, we go through the following steps:
1) Initial assessment of raw data (protein level);
  - Boxplots of raw samples (protein level);
  - Heatmap of sample correlation (protein level)
2) Normalization using median-centering of samples;
3) Quality control of normalized data (peptide and protein level)
 - Comparison of spike-in peptide quantification in raw and median-centered data
 - Boxplots of normalized samples (protein-level);
 - Coefficient of variation of proteins before and after normalization.

Steps 1) and 3) are data vizualization steps, that use the following data tables:

Therefore, we need the following data frames:
1) Sample annotation (`sample_annotation_InterLab.csv`)
2) Raw data:
  - peptide level data frame (`peptide_df_raw.csv`)
  - protein level data frame (`protein_df_raw.csv`)
3) Normalized data:
  - peptide level data frame (`peptide_df_medianCentered.csv`)
  - protein level data frame (`protein_df_medianCentered.csv`)
4) Spike-in representation in samples (`spike_ins_in_samples.csv`)

The starting tables, that are supplementary files in InterLab study (?need to check!),
are stored in `data/1_original_data` (further referred to as `1_original_data`):
1) OpenSWATH output proteome file: `all_sites_global_q_0.01_applied_to_local_global.txt`;
2) Spike-in tables:
  - `group_to_dilutionSeries.csv` 
  - `spike_in_group.csv`

These tables need to be transformed through some intermediate tables (stored in 
`data/2_interim_data`) into tables, required for plots (stored in `data/3_data_for_plots`)


Thus, the following steps are needed:
1) Sample_annotation_preparation (from `1_original_data/all_sites_global_q_0.01_applied_to_local_global.txt` to
  `3_data_for_plots/sample_annotation_InterLab.csv`);
2) Normalization:
  1. Extract fragmentome from OpenSWATH output file, used in normalization (from `1_original_data/all_sites_global_q_0.01_applied_to_local_global.txt` to `2_interim_data/all_sites_fragments_global_q_001_applied_to_local_global_full.tsv`);
  2. Normalize fragmentome (from `2_interim_data/all_sites_fragments_global_q_001_applied_to_local_global_full.tsv` to `2_interim_data/all_sites_fragments_global_q_001_applied_to_local_global_medianCentered_full.tsv`);
3) Prepare proteomic data before and after normalization for plotting: 
  1. Infer proteins from original OpenSWATH output file (from `1_original_data/all_sites_global_q_0.01_applied_to_local_global.txt` to `3_data_for_plots/proteins_df_raw.csv`)
  2. Infer proteins from normalized fragmentome (from 
  `2_interim_data/all_sites_fragments_global_q_001_applied_to_local_global_medianCentered_full.tsv` to `3_data_for_plots/protein_df_medianCentered.csv`)
4) Prepare spike-in peptide data for plotting:
  a. Infer peptides from the original OpenSWATH output file (from `1_original_data/all_sites_global_q_0.01_applied_to_local_global.txt` to `2_interim_data/peptide_df_raw.csv`);
  b. Infer peptides from normalized fragmentome (from 
  `2_interim_data/all_sites_fragments_global_q_001_applied_to_local_global_medianCentered_full.tsv` to `2_interim_data/peptide_df_medianCentered.csv`)
  c. Prepare_spike_in_table (from `1_original_data/group_to_dilutionSeries.csv`, `1_original_data/spike_in_group.csv` to `2_interim_data/spike_ins_in_samples.csv`
  d. Merge three dataframes (`2_interim_data/spike_ins_in_samples.csv`, `2_interim_data/peptide_df_raw.csv`and `2_interim_data/peptide_df_medianCentered.csv` to `3_data_for_plots/spike_ins_comparison_df.csv`).


Thus, the code is organized as follows:
1) `1_prepare_sample_annotation.R`;
2) Raw data frame preparation:
  a)`2a_infer_peptides_raw.R`
  b)`2a_infer_proteins_raw.R`
3) Normalize and prepare for plotting:
  a)`3a_normalize_fragment_level.R`
  b)`3b_infer_peptides_normalized.R`
  c)`3c_infer_proteins_normalized.R`
4) Prepare spike-in tables:
  `4_prepare_spike_in_tables.R`
5) Plotting the data for SuppFigure1:
  `5_SuppFig1_InterLab_initAssessment_QC.R`
