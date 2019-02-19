theme.strip.bar2 <- function(){
  mytheme <- theme_bw() + theme(axis.text = element_text(size = 12, colour = 'black'),
                                axis.text.x = element_text(size = 12, colour = 'black'),
                                axis.title = element_text(size = 12, colour = 'black'),
                                axis.text.y = element_text(size = 12, colour = 'black'),
                                legend.text = element_text(size = 12, colour = 'black'),
                                legend.title = element_text(size = 12, colour = 'black'),
                                strip.text = element_text(size = 12, colour = 'black'),
                                strip.background = element_rect(fill = 'white'),
                                plot.title = element_text(hjust = 0.5))
  return(mytheme)
}


theme.pie <- function(){
  mytheme <- theme_bw() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid  = element_blank(),
          strip.background = element_rect(colour = "black", fill = "white"),
          strip.text = element_text(color = 'black', size = 12.5),
          legend.text = element_text(size = 12, colour = 'black'),
          legend.title = element_text(size = 14, colour = 'black'))
  return(mytheme)
}

# theme
theme_Publication <- function(base_size=15, base_family="Helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size= unit(0.5, "cm"),
            # legend.margin = unit(0.5, "cm"),
            legend.margin = margin(5,5,5,5),
            legend.title = element_text(face="bold"),
            #plot.margin=unit(c(10,5,5,5),"mm"),
            plot.margin=unit(c(10,5,5,10),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
}