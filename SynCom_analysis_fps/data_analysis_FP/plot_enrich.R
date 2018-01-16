	#***require packages for the analysis the analysis
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "ape", "GUniFrac", "ape", "phytools", "metagenomeSeq", "ggtern","gtools")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

tab1=read.table( "HCDepletedMatrix.txt", sep="\t", row.names=1, header=T)
tab2=read.table( "HCDepletedRoot.txt", sep="\t", row.names=1, header=T)
tab3=read.table( "HSDepletedMatrix.txt", sep="\t", row.names=1, header=T)
tab4=read.table( "HSDepletedRoot.txt", sep="\t", row.names=1, header=T)


tab=smartbind(tab1, tab2, tab3, tab4HCDepletedMatrix.txt)
