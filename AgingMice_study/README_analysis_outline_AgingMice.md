This is the code, accompanying the batch effects analysis in the Aging Mice study.
This study was used as a pilot dataset for the workflow, with certain approaches 
designed and first tested with this dataset. Thus, almost all of the code is aligned
with `proBatch` package.

In this study, we go through the following steps:
1) Initial assessment of raw data (peptide level): 
  a. Mean intensity of samples, sorted in MS run order;
  b. Distribution of sample correlations;
  c. Example of QTL that cannot be detected due to batch effects;
  d. Boxplots of raw samples;
2) Normalization using quantile normalization;
  - Boxplots of normalized dataset
3) Intermediate diagnostics:
  a. Plot of two Principal Components, colored by:
    - MS batch;
    - Replicate ID;
  b. Hierarchical clustering of samples;
  c. Principal Variance Component Analysis;
  d. Single-feature plot examples.
4) Correction of batch effects by LOESS curve fitting + median centering;
  - Illustration of correction on two representative peptides;
5) Quality control of the batch adjustment:
  a. Example of QTL, rescued by batch effect adjustment;
  b. Distribution of sample correlations after batch effect adjustment;
  c. Heatmap of correlation of peptides of the same or different proteins;
  d. Distribution of peptide correlation before and after correction.
  
Additionally, we examine the batch effects and missing values relationship. For this, we check:
  - Hierarchical clustering;
  - Correlation heatmap of 10 selected replicated samples with and without requant values;
  - Correlation distribution of 10 selected replicated samples with and without requant values;

The following data tables are used to illustrate each step:

Therefore, we need the following data frames:
1) Sample annotation (`sample_annotation_AgingMice.csv`)
2) Raw data:
  - peptide-level data without requant values (`raw_proteome_AgingMice.csv`)
  - peptide-level data with requant values (`raw_proteome_AgingMice_with_requants.csv`)
3) Normalized data:
  - peptide-level data frame (`normalized_proteome_AgingMice.csv`)
4) Batch-corrected data:
  - peptide-level data frame (`batchCorrected_proteome_AgingMice.csv.csv`)
5) Peptide and Allele annotation data:
  - `peptide_annotation_df.csv`;
  - `allelle_annotation_df.csv`.
6) Data frames with "requant" values: TODO: check, whether we can plot everything for normalized data only.
- `raw_proteome_AgingMice_with_requants.csv`;
- `quantile_normalized_matrix_agingMice_with_requants.csv` - TODO: what happens, if we cluster unnormalized data?
supplemented with `back_to_back_samples.csv`;

The starting tables from Aging Mice (are they supplements? need to check! with Evan),
are stored in `data_PanCancer/1_original_data` (further referred to as `1_original_data`):
1) OpenSWATH output proteome file: `E1801171630_feature_alignment_requant.tsv`;
2) Sample annotation: `sample_annotation_E1801171630.csv`;
3) genotype annotation table: `Protein_Gene_LOCATION_AND_ALLELE_SYMBOL_v2.csv`;
4) peptide to Gene mapping: `peptide_annotation.csv`;

These tables need to be transformed through some intermediate tables (stored in 
`data_AgingMice/2_interim_data`) into tables, required for plots (stored in `data_AgingMice/3_data_for_plots`). 
Some diagnostic plots, that have large memory requirements,
are calculated on high-capacity machine and the results stored in `plots_AgingMice/interim_data_for_plots`


Thus, the following steps are needed:
1) Data preparation
  a.  Sample_annotation_preparation (from `1_original_data/sample_annotation_E1801171630.csv` to
  `3_data_for_plots/sample_annotation_AgingMice.csv`);
  b. Prepare raw data for the batch adjustment - normalization and batch correction (from `1_original_data/E1801171630_feature_alignment_requant.csv` to `2_interim_data/raw_proteome_AgingMice.csv`);
  c. Prepare color palette for all the plots (from `3_data_for_plots/sample_annotation_AgingMice.csv` to `3_data_for_plots/color_annotation.rda`)
2) Normalization on peptide level: 
  - for the data without "requant" (from `3_data_for_plots/raw_proteome_Aging_mice.csv` to `2_interim_data/normalized_proteome_AgingMice.csv`);
  - for the data with "requant" values (from `3_data_for_plots/raw_proteome_AgingMice_with_requants.csv` to `3_data_for_plots/normalized_proteome_AgingMice_with_requants.csv`)
3) Correct batch effects on peptide level: 
  - correction of trends, used to illustrate the fit: (from `2_interim_data/` to `3_data_for_plots/adjusted_fit_df_agingMice.csv`);
  - correction of discrete effects (from `3_data_for_plots/adjusted_fit_df_agingMice.csv` to `3_data_for_plots/batchCorrected_proteome_AgingMice.csv`);

Additionally, some resource-demanding visualizations are prepared separately:
1) Principal Variance Component Analysis:
  - for intermediate diagnostics, PVCA of the peptides, represented in at least 
  50\% of the samples is calculated (from `3_data_for_plots/normalized_proteome_AgingMice_with_requants.csv` to `pvca_normalized_50.csv`);
  - for missing value effect illustration, additionally, several PVCA calculations are prepared:
  TODO: decide if we are going to do this in practice!
2) Peptide correlation between and within proteins:
 - `peptide_cor_raw.csv`
 - `peptide_cor_corrected.csv`


Thus, the code is organized as follows:
1) Prepare the data:
  a. `1a_prepare_sample_annotation.R`;
  b. `1b_prepare_raw_proteome.R`;
  c. `1c_prepare_allelle_table.R`
2) Normalization:
  - `2a_normalize_AgingMice.R`;
  - `2b_normalize_with_requants.R`
3) Correct batch effects:
  - `3a_fit_LOESS_curve.R`
  - `3b_correct_medianCentering.R`
4) Resource-heavy code used in plots (using the files 
in `data_AgingMice/3_data_for_plots`, storing the results in `plots_PanCancer/:
  a. `4a_PVCA_analysis`;
  b. `4b_peptide_correlation_raw.R`;
  c. `4c_peptide_correlation_batchCorrected.R`
5) Plotting the figures (using the files 
in `data_AgingMice/3_data_for_plots`, storing the results in `plots_AgingMice/interim_data_for_plots`):
 - `5a_Fig2_initial_assessment_normalization.R`;
 - `5b_Fig3_diagnostic_plots.R`;
 - `5c_Fig4_missing_value_effect.R`;
 - `5d_Fig5_curve_fitting.R`;
 - `5e_Fig6_quality_control.R`;
 - `5f_SuppFig3a_sample_correlation_heatmap.R`;
 - `5f_SuppFig3b_span_demonstration.R`;
 - `5f_SuppFig3c_QC_hierarchical_clustering.R`; 
 - `5g_SuppFig4_missing_values_extras.R`;
