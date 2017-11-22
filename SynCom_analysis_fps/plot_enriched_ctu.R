#R scripts used in Hassani et al., to plot heatmap of significantly abundant CTU (example of FlowPots)
#M. Amine Hassani - ahassani@bot.uni-kiel.de

#required packages for the analysis
pkg=c("ggplot2","data.table","reshape2", "magrittr")
lapply(pkg, library, character.only = TRUE)

#clean
rm(list = ls()[grep("*.*", ls())])

#define plotting form
theme_change <- theme(
		plot.background = element_blank(),
		panel.grid.minor = element_blank(),
		panel.grid.major = element_blank(),
		panel.background = element_blank())
color_palette <- colorRampPalette(c("#2163d4","#21abd4" ,"#e6ee71","#eeb811","#ee1c11"))(5)
color_palette1 <- colorRampPalette(c("#ee1c11","#eeb811","#e6ee71","#21abd4" ,"#2163d4"))(5)

#define order in x- and y-axis
order1<-c("HCDepletedMatrix","HCDepletedRoot","HSDepletedMatrix","HSDepletedRoot")
order2<-c("Root369","Root265","Soil538","Root682","Soil805","Root344","Root151","Soil761","Root101","Root60","Root322","Root61","Root1280",
"Root329","Root9","Root76","Root179","Root627","Soil773","Root170","Root189","Root16D2","Root318D1","Root70","Root219","Root267","Root487D2Y",
"Root700","Root710","Root381","Root552","Root635")

#upload data
mat=read.table( "enriched_ctu_fps.txt", sep="\t", row.names=1, header=T)

#reformat the matrix
DT <- mat %>% 
as.data.table(keep.rownames=TRUE) %>% 
melt(id.vars = "rn", value.name="sz")

#categorize rel. abund. values in classes
v <- c(-4,-2,-1,1,2,5)
DT$cat<- findInterval(DT$sz, v)

#prepare the plot
p=ggplot(DT, aes(y = rn, x =variable , fill = factor(cat)))
p$data$variable <- ordered(p$data$variable, levels= order1)
p$data$rn <- ordered(p$data$rn, levels= rev(order2))
p1=p+geom_tile(color = "white") +scale_fill_manual(values = color_palette1)+theme_change
print(p1)
