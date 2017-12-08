
# script to reproduce the association tests (phyloGLM) and
# corresponding figure (Fig. 2) from Hassani et al., 2017
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

library(ggplot2)
library(reshape)
library(scales)
library(grid)
library(gridExtra)

colAddAlpha <- function(col, alpha=1) { 
    
    x <- as.vector(col2rgb(col)) / 255
    
    return(rgb(x[1], x[2], x[3], alpha))

}

# colors

alpha <- 1
c_yellow <-          rgb(255 / 255, 255 / 255,   0 / 255, alpha)
c_blue <-            rgb(  0 / 255, 000 / 255, 255 / 255, alpha)
c_orange <-          rgb(255 / 255,  69 / 255,   0 / 255, alpha)
c_green <-           rgb(  50/ 255, 220 / 255,  50 / 255, alpha)
c_dark_green <-      rgb( 50 / 255, 200 / 255, 100 / 255, alpha)
c_very_dark_green <- rgb( 50 / 255, 150 / 255, 100 / 255, alpha)
c_sea_green <-       rgb( 46 / 255, 129 / 255,  90 / 255, alpha)
c_black <-           rgb(  0 / 255,   0 / 255,   0 / 255, alpha)
c_grey <-            rgb(180 / 255, 180 / 255,  180 / 255, alpha)
c_dark_brown <-      rgb(101 / 255,  67 / 255,  33 / 255, alpha)
c_red <-             rgb(200 / 255,   0 / 255,   0 / 255, alpha)
c_dark_red <-        rgb(255 / 255, 130 / 255,   0 / 255, alpha)

# ggplot2 theme

main_theme <- theme(panel.background=element_blank(),
              panel.grid=element_blank(),
              axis.line.x=element_line(color="black"),
              axis.line.y=element_line(color="black"),
              axis.ticks=element_line(color="black"),
              axis.text=element_text(colour="black", size=10),
              legend.position="top",
              legend.background=element_blank(),
              legend.key=element_blank(),
              text=element_text(family="sans"))

