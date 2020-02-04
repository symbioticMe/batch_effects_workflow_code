library(tidyr)

transition_to_fragment <- function(proteome_long, 
                                  fragment_annotation_column = 'aggr_Fragment_Annotation',
                                  fragment_intensity_column = 'aggr_Peak_Area'){
  print(names(proteome_long))
  print(fragment_annotation_column %in% names(proteome_long))
  print(fragment_intensity_column %in% names(proteome_long))
  ion_id_columns = paste('ionID', 1:6, sep = '_')
  ion_int_columns = paste('ionIntensity', 1:6, sep = '_')
  transitome = proteome_long %>%
    separate_(fragment_annotation_column, into = ion_id_columns, sep = ';', remove = F) %>%
    separate_(fragment_intensity_column, into = ion_int_columns, sep = ';', remove = F) %>%
    data.table::as.data.table() %>%
    data.table::melt( measure = list(ion_int_columns, ion_id_columns), 
                      variable.name = 'ions', 
                      value.name = c('Ion_intensity', 'Ion_ID')) %>%
    mutate(Ion_intensity = as.numeric(Ion_intensity))
  return(transitome)
}

fragment_df_to_openSWATH <- function(transitome, 
                                  fragment_annotation_column = 'Ion_ID',
                                  fragment_intensity_column = 'Ion_intensity',
                                  id_column = 'ions',
                                  fragment_united_column = 'aggr_Fragment_Annotation_new',
                                  fragment_united_int_column = 'aggr_Peak_Area_new',
                                  un_log = 2,
                                  intensities_to_exclude = NULL){
  if(!is.null(un_log) & is.numeric(un_log)){
    transitome = transitome %>%
      mutate(rlang::UQ(rlang::sym(fragment_intensity_column)) := un_log^rlang::UQ(rlang::sym(fragment_intensity_column)))
  }
  if (!is.null(intensities_to_exclude)){
    transitome = transitome %>%
      select(-one_of(intensities_to_exclude))
  }
  id_columns = setdiff(names(transitome), 
                       c(fragment_annotation_column, fragment_intensity_column, id_column))
  id_for_dcast = paste(id_columns, collapse = ' + ')
  formula_for_dcast = paste(c(id_for_dcast, id_column), collapse = ' ~ ')
  transitome_casted = transitome %>%
    data.table::as.data.table() %>%
    data.table::dcast(formula_for_dcast, 
          value.var = c(fragment_annotation_column, fragment_intensity_column))
  ion_id_columns = names(transitome_casted)[grepl(fragment_annotation_column, names(transitome_casted))]
  ion_int_columns = names(transitome_casted)[grepl(fragment_intensity_column, names(transitome_casted))]
  transitome_united = transitome_casted %>%
    as.data.frame() %>%
    unite_(fragment_united_column, ion_id_columns, sep = ';') %>%
    unite_(fragment_united_int_column, ion_int_columns, sep = ';')
  transitome_united[[fragment_united_column]] = gsub('(;NA)?','', 
                                                     transitome_united[[fragment_united_column]])
  
  transitome_united[[fragment_united_int_column]] = gsub('(;NA)?','', 
                                                     transitome_united[[fragment_united_int_columnhel]])
  return(transitome_united)
}
