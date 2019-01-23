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
