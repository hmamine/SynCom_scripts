#***require packages for the analysis the analysis
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "ape", "GUniFrac", "ape", "phytools", "metagenomeSeq", "ggtern","RColorBrewer")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])
#		***Upload and prepare phyloseq objects***
myPalette <-c("#EF5656","#47B3DA","#F7A415","#2BB065")
#myPalette <-c("#EF5656","#F7A415","#2BB065")

mat=read.table( "table_input.txt", sep="\t", row.names=1, header=T)
OTU=otu_table(mat, taxa_are_rows=T) 
tax=read.table("taxa_table.txt", sep="\t", row.names=1, header=1)
tax=as.matrix(tax)
TAXA=tax_table(tax)
sd=read.table("sample_data.txt", sep="\t", row.names=1, header=1)
SD=sample_data(sd)
physeq= phyloseq(OTU, TAXA, SD) 


ph1=transform_sample_counts(physeq, function(x) x / sum(x) )

otumat2=as(otu_table(ph1), "matrix")
taxmat2=as(tax_table(ph1), "matrix")
DT_otumat2=data.table(otumat2, keep.rownames=T, key="rn")
DT_taxmat2=data.table(taxmat2, keep.rownames=T, key="rn")
comb_mat=merge(DT_taxmat2, DT_otumat2)
DT=melt(comb_mat)

p=ggplot(DT, aes(x=value, y=rn, color=Phylum))
p$data$rn <- factor(p$data$rn, levels = p$data$rn[order(-p$data$value)])


theme_change <- theme(
		plot.background = element_blank(),
		panel.grid.minor = element_blank(),
		panel.grid.major = element_blank(),
		panel.background = element_blank(),
		panel.border = element_blank(),
		axis.line = element_blank(),
		axis.ticks = element_blank(),
		axis.text.y=element_blank())

p1=p+geom_point(size=2)+theme_change+scale_colour_manual(values=myPalette)+geom_vline(xintercept = 0.001)

pdf("fig_input.pdf", useDingbats=FALSE,paper="A4")
plot(p1)
dev.off()
