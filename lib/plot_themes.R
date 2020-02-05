theme_publication <- function() {
  
  new_theme = theme(plot.title = element_text(face = "bold",
                                              size = rel(2), hjust = 0.5),
                    axis.title = element_text(face = "bold",size = rel(1.5)),
                    axis.title.y = element_text(angle=90,vjust =2),
                    axis.title.x = element_text(vjust = -0.2),
                    axis.text = element_text(size = rel(1.2)), 
                    axis.line = element_line(colour="black"),
                    axis.ticks = element_line())
  return(new_theme)  
}

theme_publication_mild <- function() {
  
  new_theme = theme(plot.title = element_text(face = "bold",
                                              size = 16, hjust = 0.5),
                    axis.title = element_text(face = "bold",size = 14),
                    axis.title.y = element_text(angle=90,vjust =2),
                    axis.title.x = element_text(vjust = -0.2),
                    axis.text = element_text(size = 8), 
                    axis.line = element_line(colour="black"),
                    axis.ticks = element_line(),
                    plot.margin = margin(t=1, l=1, r=.251, b=.251, unit = "cm"))
  return(new_theme)  
}

theme_cute_legend <- function(){
  new_theme = theme(legend.key = element_rect(colour = NA),
                    legend.position = "bottom",
                    legend.direction = "horizontal",
                    legend.key.size= unit(0.5, "cm"),
                    legend.title = element_text(face="italic"))
  return(new_theme)
}