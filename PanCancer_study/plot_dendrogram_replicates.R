plot_dendrogram_replicates <- function(matrix_proBatch_Tatjana, 
                            sample_annotation = sample_annotation_PanCancer, 
                            color_list = color_list_PanCancer,
                            biospecimen_id_col = 'patient_IDs',
                            sample_id_col ='sample_name', plot_title = NULL) {
  
  colors_for_patients = color_list[[biospecimen_id_col]]
  dendrogram <- matrix_proBatch_Tatjana %>% t() %>% 
    dist(method = 'manhattan') %>% hclust() %>% as.dendrogram()
  dendro_data <- dendro_data(dendrogram)
  labels_dendro <- label(dendro_data)
  labels_dendro = merge(labels_dendro, sample_annotation, 
                        by.x = 'label',by.y = sample_id_col) 
  gg_dendro <-  ggplot(segment(dendro_data)) +
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend), color = 'grey40')
  gg_dendro = gg_dendro + geom_text(data=labels_dendro,
                                    aes_string(label='label', x='x', y=0, colour=biospecimen_id_col), angle = 90, 
                                    vjust = "top", hjust = 'right', size = 2)+
    scale_color_manual(values = colors_for_patients, guide=FALSE)+ theme_dendro()+
    theme(legend.position = "none")+ylim(c(-50, 250))
  if (!is.null(plot_title)){
    gg_dendro = gg_dendro + ggtitle(plot_title)
  }
  return(gg_dendro)
}