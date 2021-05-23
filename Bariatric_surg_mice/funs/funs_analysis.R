plot_PVCA.df <- function(pvca_res,  colors_for_bars = NULL,  filename = NULL, width = NA, height = NA, 
                         units = c('cm','in','mm'), plot_title = NULL,  theme = 'classic',
                         base_size = 20){
  pvca_res = pvca_res %>%
    mutate(label = factor(label, levels=label))
  
  y_title = 'Weighted average proportion variance'
  gg  = ggplot(pvca_res, aes(x = label, y = weights, fill = category))+
    geom_bar(stat = 'identity', color = 'black')+
    ylab(y_title)
  
  
  if(is.null(colors_for_bars)){
    colors_for_bars = c('grey', wes_palettes$Rushmore[3:5])
    names(colors_for_bars) = c('residual', 'biological', 
                               'biol:techn', 'technical')
    
  } else {
    if (length(colors_for_bars) != 4){
      color_names = paste(c('residual', 'biological', 'biol:techn', 
                            'technical'), collapse = ' ')
      warning(sprintf('four colors for: %s were expected', color_names))
    }
  }
  gg = gg + scale_fill_manual(values = colors_for_bars)
  
  if (!is.null(plot_title)){
    gg = gg + ggtitle(plot_title)
  }
  
  #Change the theme
  if(!is.null(theme) && theme == 'classic'){
    gg = gg + theme_classic(base_size = base_size)
  }else{
    message("plotting with default ggplot theme, only theme = 'classic' 
            implemented")
  }
  
  gg = gg +
    theme(axis.title.x = NULL, 
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))+
    xlab(NULL)+
    theme(text = element_text(size=15))+
    guides(fill=guide_legend(override.aes=list(color=NA), title=NULL))
  
  return(gg)
}

# Note: due to low # of samples for animal:time, plate_col:plate_row, I computed PVCA for shorlisted covariates and combined variances 
# from multiple computations, averaging overlapping variances. Please refer to pvca_analysis_for_dda_data for details. 
pvca_analysis_for_dda_data <- function(matrix){
  
  technical_factors_1 = c('MS_batch', 'plate_col')
  technical_factors_2 = c('MS_batch', 'plate_row')
  biological_factors_1 = c('animal', 'condition')
  biological_factors_2 = c('time', 'condition')
  
  pvca_1 <- plot_PVCA(matrix, sample_annotation,
                      technical_factors = technical_factors_1,   biological_factors = biological_factors_1)
  pvca_2 <- plot_PVCA(matrix, sample_annotation,
                      technical_factors = technical_factors_1,   biological_factors = biological_factors_2)
  pvca_3 <- plot_PVCA(matrix, sample_annotation,
                      technical_factors = technical_factors_2,   biological_factors = biological_factors_1)
  pvca_4 <- plot_PVCA(matrix, sample_annotation,
                      technical_factors = technical_factors_2,   biological_factors = biological_factors_2)
  
  pvca_averaged_weights <- pvca_1$data %>% rbind(pvca_2$data) %>% rbind(pvca_3$data) %>% rbind(pvca_4$data) %>%
    group_by(label, category) %>% dplyr::summarise(weights = mean(weights)) %>% 
    mutate_if(is.factor, as.character) %>% dplyr::arrange(-weights) 
  labels <- as.character(pvca_averaged_weights$label)
  labels <- c(labels[!labels %in% c('resid', 'Below 1%')], c( 'Below 1%', 'resid'))
  pvca_averaged_weights <- pvca_averaged_weights %>% mutate(label = factor(label, levels = labels))
  
  plot_PVCA.df(pvca_averaged_weights) 
}
gg_PVCA <- pvca_analysis_for_dda_data(quantile_normalized_matrix)