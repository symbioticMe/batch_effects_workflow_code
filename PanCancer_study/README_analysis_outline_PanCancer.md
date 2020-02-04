This is the code, accompanying the batch effects analysis in the PanCancer study.
This study was used as an example to define the workflow, so all steps are executed
exactly as they were in the original study, with minor polishing of the code to
`proBatch` format and phylosophy.

In this study, we go through the following steps:
1) Initial assessment of raw data (protein level): Boxplots of raw samples 
2) Normalization using quantile normalization;
3) Intermediate diagnostics:
  - Control of normalization by plotting boxplots of normalized data;
  - Diagnostics of individual features (spike-in protein in normalized, but not batch-corrected data);
  - Diagnostics of whole proteome using clustering enhanced by the heatmap (normalized, but not batch-corrected data.)
4) Correction of batch effect by per-feature median centering;
5) Quality control (protein level)
 - Comparison of spike-in protein abundance in normalized and batch-corrected data;
 - Comparison of clustering enhanced by heatmaps in normalized and batch-corrected data;
 - Comparison of replicate clustering in normalized and batch-corrected data.

Steps 1) and 4) are data vizualization steps, that use the following data tables:

Therefore, we need the following data frames:
1) Sample annotation (`sample_annotation_PanCancer.csv`)
2) Raw data:
  - transition-level data frame (`raw_transitome_PanCancer.csv`)
  - protein-level data frame (`proteome_df_raw.csv`)
3) Normalized data:
  - transition-level data frame (`normalized_transitome_PanCancer.csv`)
  - protein-level data frame (`proteome_df_normalized.csv`)
4) Batch-corrected data:
  - transition-level data frame (`batchCorrected_transitome_PanCancer.csv`);
  - protein-level data frame (`proteome_df_corrected.csv`).
  

The starting tables from PanCancer (are they supplements? need to check!),
are stored in `data_PanCancer/1_original_data` (further referred to as `1_original_data`):
1) OpenSWATH output proteome file: `feature_alignment_modB_trs.csv`;
2) Sample annotation: `annotation_col_anova_PCP_final.txt`

These tables need to be transformed through some intermediate tables (stored in 
`data_PanCancer/2_interim_data`) into tables, required for plots (stored in `data_PanCancer/3_data_for_plots`)


Thus, the following steps are needed:
1) Data preparation
  a.  Sample_annotation_preparation (from `1_original_data/annotation_col_anova_PCP_final.txt` to
  `3_data_for_plots/sample_annotation_PanCancer.csv`);
  b. Prepare raw data for normalization and protein inference (from `1_original_data/feature_alignment_modB_trs.csv` to `2_interim_data/raw_transitome_PanCancer.csv`);
2) Normalization on fragment level: (from `2_interim_data/raw_transitome_PanCancer.csv` to `2_interim_data/normalized_transitome_PanCancer.csv`);
3) Correct batch effects on fragment level: (from `2_interim_data/normalized_transitome_PanCancer.csv` to `2_interim_data/batchCorrected_transitome_PanCancer.csv`);
4) Infer proteins from transitions as these are used for plotting: 
  a. Infer proteins from raw transitome (from `2_interim_data/raw_transitome_PanCancer.csv` to `3_data_for_plots/proteome_df_raw.csv`)
  b. Infer proteins from normalized transitome (from 
  `2_interim_data/normalized_transitome_PanCancer.csv` to `3_data_for_plots/proteome_df_normalized.csv`);
  c. Infer proteins from normalized transitome (from 
  `2_interim_data/normalized_transitome_PanCancer.csv` to `3_data_for_plots/proteome_df_normalized.csv`);
5) Plot the initial assessment, diagnostic and quality control plots (using the files 
in `data_PanCancer/3_data_for_plots`, storing the results in `plots_PanCancer`)


Thus, the code is organized as follows:
1) Prepare the data:
  a. `1a_prepare_sample_annotation.R`;
  b. `1b_prepare_raw_data_PanCancer.R`;
2) Normalization:
  `2_normalize.R`
3) Correct batch effects:
  `3_correct_batch_effects.R`
4) Infer proteins:
  a. `4a_infer_proteins_raw.R`;
  b. `4a_infer_proteins_normalized.R`
  c. `4a_infer_proteins_corrected.R`
5) Plotting the data for SuppFigure1:
  `5_SuppFig2_PanCancer_initAssessment_QC.R`
