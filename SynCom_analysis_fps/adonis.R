#R scripts to reproduce beta diversity analyses used in Hassani et al., (example of FlowPots)
#M. Amine Hassani - ahassani@bot.uni-kiel.de

#required packages for the analysis
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "ape", "GUniFrac", "ape", "phytools", "metagenomeSeq", "ggtern","RColorBrewer","gdata")
lapply(pkg, require, character.only = TRUE)

#clean
rm(list = ls()[grep("*.*", ls())])

#defining a plotting form
theme_change <- theme(
plot.background = element_blank(),
panel.grid.minor = element_blank(),
panel.grid.major = element_blank(),
panel.background = element_blank())

#***Upload and prepare phyloseq objects***
mat=read.table( "count_table_fps.txt", sep="\t", row.names=1, header=T)
mat=as.matrix(mat)
OTU=otu_table(mat, taxa_are_rows=T) 
tax=read.table("taxonomy_table.txt", sep="\t", row.names=1, header=1)
tax=as.matrix(tax)
TAXA=tax_table(tax)
sd=read.table("sample_data_fps.txt", sep="\t", row.names=1, header=1)
SD=sample_data(sd)
TREE=read.tree("16s_trim_v5v7_KG_DepExp_Derep.tree")
physeq= phyloseq(OTU, TAXA, SD,TREE) 

#remove samples with less than 1000 reads
#physeq = prune_samples(sample_sums(physeq) > 1000, physeq)

#subset samples to root and matrix
cond=c("matrix","root")
physeq=subset_samples(physeq, config1 %in% cond)

#count data normalization 
otumat=as(otu_table(physeq), "matrix")
mp=newMRexperiment((otumat))
physeq=phyloseq(otu_table(MRcounts(cumNorm(mp, p=cumNormStat(mp)), norm=TRUE, log=TRUE),taxa_are_rows=TRUE), TAXA, SD,TREE )

#subset samples to full community and highly sensitive depleted samples
cond=c("afull","hs")
physeq.hs=subset_samples(physeq, config2 %in% cond)
comm=c("Back","Comp")
physeq.hs=subset_taxa(physeq.hs, tag %in% comm)

#subset to matrix samples
submat=c("matrix")
physeq.hsmat=subset_samples(physeq.hs, config1 %in% submat)

#prepare data table to compute distances
sdtab_hsmat=data.table(as(sample_data(physeq.hsmat), "data.frame"),keep.rownames=T, key="rn")
otumat_hsmat=t(as(otu_table(physeq.hsmat), "matrix") )

#BrayCurtis distance (BC)
DB.hsmat=vegdist(otumat_hsmat, method="bray")

#weighted UniFrac distance (wUF)
#DB.hsmat <-UniFrac(physeq.hsmat, weighted=TRUE)

#permutational multivariate analysis of variance using distance matrices
adB.hsmatp=with(sdtab_hsmat, adonis(DB.hsmat ~ config3))

#storing data
tmp<-data.frame(adB.hsmatp$aov.tab$R2)
colnames(tmp)<-"hs.matrix"
x1<-as.data.frame(adB.hsmatp$aov.tab[6])
p1=x1[1,1]
tmp= rbind(tmp, p1)
rownames(tmp)<-c("explained","residuals","total", "pvalue")
hs.matrix=tmp

#subset root samples
subrot=c("root")
physeq.hsrot=subset_samples(physeq.hs, config1 %in% subrot)

#prepare data table to compute distances
sdtab_hsrot=data.table(as(sample_data(physeq.hsrot), "data.frame"),keep.rownames=T, key="rn")
otumat_hsrot=t(as(otu_table(physeq.hsrot), "matrix") )

#BrayCurtis distance (BC)
DB.hsrot=vegdist(otumat_hsrot, method="bray")

#weighted UniFrac distance (wUF)
#DB.hsrot <-UniFrac(physeq.hsrot, weighted=TRUE)

#permutational multivariate analysis of variance using distance matrices
adB.hsrotp=with(sdtab_hsrot, adonis(DB.hsrot ~ config3))

#storing data
tmp<-data.frame(adB.hsrotp$aov.tab$R2)
colnames(tmp)<-"hs.root"
x2<-as.data.frame(adB.hsrotp$aov.tab[6])
p2=x2[1,1]
tmp= rbind(tmp, p2)
rownames(tmp)<-c("explained","residuals","total", "pvalue")
hs.root=tmp

#subset full community and highly competitive depleted samples
cond=c("afull","hc")
physeq.hc=subset_samples(physeq, config2 %in% cond)
comm=c("Back","Sens")
physeq.hc=subset_taxa(physeq.hc, tag %in% comm)

#subset matrix samples
physeq.hcmat=subset_samples(physeq.hc, config1 %in% submat)

#prepare data table to compute distances
sdtab_hcmat=data.table(as(sample_data(physeq.hcmat), "data.frame"),keep.rownames=T, key="rn")
otumat_hcmat=t(as(otu_table(physeq.hcmat), "matrix") )

#BrayCurtis distance
DB.hcmat=vegdist(otumat_hcmat, method="bray")

#weighted UniFrax distance
#DB.hcmat <-UniFrac(physeq.hcmat, weighted=TRUE)

#permutational multivariate analysis of variance using distance matrices
adB.hcmatp=with(sdtab_hcmat, adonis(DB.hcmat ~ config3))

#storing data
tmp<-data.frame(adB.hcmatp$aov.tab$R2)
colnames(tmp)<-"hc.matrix"
x3<-as.data.frame(adB.hcmatp$aov.tab[6])
p3=x3[1,1]
tmp= rbind(tmp, p3)
rownames(tmp)<-c("explained","residuals","total", "pvalue")
hc.matrix=tmp

#subset to root samples
physeq.hcrot=subset_samples(physeq.hc, config1 %in% subrot)

#prepare data table to compute distances
sdtab_hcrot=data.table(as(sample_data(physeq.hcrot), "data.frame"),keep.rownames=T, key="rn")
otumat_hcrot=t(as(otu_table(physeq.hcrot), "matrix") )

#BrayCurtis distance
DB.hcrot=vegdist(otumat_hcrot, method="bray")
#weighted UniFrac distance
#DB.hcrot <-UniFrac(physeq.hcrot, weighted=TRUE)

#permutational multivariate analysis of variance using distance matrices
adB.hcrotp=with(sdtab_hcrot, adonis(DB.hcrot ~ config3))

#storing data
tmp<-data.frame(adB.hcrotp$aov.tab$R2)
colnames(tmp)<-"hc.root"
x4<-as.data.frame(adB.hcrotp$aov.tab[6])
p4=x4[1,1]
tmp= rbind(tmp, p4)
rownames(tmp)<-c("explained","residuals","total", "pvalue")
hc.root=tmp

#ploting
lst<-list (adB.hcmatp,adB.hsmatp,adB.hcrotp,adB.hsrotp)
for (i in lst) {print(i) }

tab<-cbind(hc.matrix,hs.matrix,hc.root,hs.root)
tab1=t(tab[-4:-3,])
DT=melt(data.table(tab1, keep.rownames=T, key="rn"))

	p=ggplot(DT, aes (x= rn, y=value, fill=variable))
	neworder<-c("hc.matrix","hs.matrix","hc.root","hs.root")
	p$data$rn<-ordered(p$data$rn, levels= neworder)
	p1=p+theme_bw()+geom_bar(stat="identity")+geom_label(aes(label= value))+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none")+scale_y_continuous(limits=c(0, 1))

gridExtra::grid.arrange(p1,nrow=2)
