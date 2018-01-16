#***require packages for the analysis the analysis
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "ape", "GUniFrac", "ape", "phytools", "metagenomeSeq", "ggtern","PMCMR")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

alphaInd = c("Shannon", "Observed")
color_palette<-c("#000000","#806600","#803300","666666")

#		***Upload and prepare phyloseq objects***
mat=read.table( "otu_table.txt", sep="\t", row.names=1, header=T)
mat=as.matrix(mat)
OTU=otu_table(mat, taxa_are_rows=T) 
tax=read.table("taxa_table.txt", sep="\t", row.names=1, header=1)
tax=as.matrix(tax)
TAXA=tax_table(tax)
sd=read.table("sample_data.txt", sep="\t", row.names=1, header=1)
SD=sample_data(sd)
TREE=read.tree("16s_trim_v5v6_KG_DepExp_Derep.tree")

physeq= phyloseq(OTU, TAXA, SD,TREE) 

#Remove samples with less than 1000 reads
physeq1 = prune_samples(sample_sums(physeq) > 1000, physeq)

#Select for full and hc depleted samples
cond=c("input","afull","hc")
physeq2=subset_samples(physeq1, config2 %in% cond)

comm=c("Back","Sens")
physeq2=subset_taxa(physeq2, tag %in% comm)

#Rarefication to even depth based on smallest sample size
rf=min(sample_sums(physeq2))
physeq.rrf=rarefy_even_depth(physeq2, rf, replace=TRUE)
	p_nat=plot_richness(physeq.rrf,"config3","config2" ,measures=alphaInd)
	p_nat$layers <- p_nat$layers[-1]

#neworder<-c("input","matrixpure","matrix3per","rootpure","root3per")

	#p_nat$data$config3 <- ordered(p_nat$data$config3, levels= neworder)
	p1=p_nat+geom_boxplot(data=p_nat$data, aes(x=config3, y=value, color=config2), width=1)+
	geom_point(size=4)+theme_bw()+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
	scale_colour_manual(values=color_palette)
	#ggtitle(paste0("Alpha diversity rarefication to ",rf, " sequencing depth"))

subp2 =	ggplot(data=p_nat$data[p_nat$data$variable=="Observed",], aes(x=config3, y=value, color=config2))+ geom_boxplot(width=1)+
	geom_point(size=2)+theme_bw()+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none")+
	scale_colour_manual(values=color_palette) + facet_wrap(~variable) + scale_y_continuous(limits=c(0, 130), breaks=c(0,25,50,75,100,125)) 
kruskal.test(data=subp2$data, value ~ config3)
posthoc.kruskal.conover.test(data=subp2$data, value ~ config3, method="BH")

subp1 =	ggplot(data=p_nat$data[p_nat$data$variable=="Shannon",], aes(x=config3, y=value, color=config2))+ geom_boxplot(width=1)+
	geom_point(size=2)+theme_bw()+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none")+
	scale_colour_manual(values=color_palette) + facet_wrap(~variable) + scale_y_continuous(limits=c(0, 6), breaks=c(0,0.5,1.5,2.5,3.5,4.5,5.5))
kruskal.test(data=subp1$data, value ~ config3)
posthoc.kruskal.conover.test(data=subp1$data, value ~ config3)

gridExtra::grid.arrange(subp2, subp1, ncol=2)
















