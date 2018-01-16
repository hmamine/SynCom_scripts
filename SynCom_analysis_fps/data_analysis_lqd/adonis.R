#require packages for the analysis the analysis
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "ape", "GUniFrac", "ape", "phytools", "metagenomeSeq", "ggtern","RColorBrewer","gdata")
lapply(pkg, require, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

subSt=c("strtsb")
subSk=c("unstrtsb")

#shapes-colors-forms
theme_change <- theme(
plot.background = element_blank(),
panel.grid.minor = element_blank(),
panel.grid.major = element_blank(),
panel.background = element_blank())

###***Upload and prepare phyloseq objects***
mat=read.table( "ctu_table_lqd.txt", sep="\t", row.names=1, header=T)
mat=as.matrix(mat)
OTU=otu_table(mat, taxa_are_rows=T) 
tax=read.table("taxa_table_lqd.txt", sep="\t", row.names=1, header=1)
tax=as.matrix(tax)
TAXA=tax_table(tax)
sd=read.table("sample_data_lqd.txt", sep="\t", row.names=1, header=1)
SD=sample_data(sd)
TREE=read.tree("16s_trim_v5v7_KG_DepExp_Derep.tree")

physeq= phyloseq(OTU, TAXA, SD,TREE) 

#physeq = prune_samples(sample_sums(physeq) > 1000, physeq)
cond=c("afull","hc", "hs")
physeq=subset_samples(physeq, config2 %in% cond)

#normalization of count reads 
otumat=as(otu_table(physeq), "matrix")
mp=newMRexperiment((otumat))
physeq=phyloseq(otu_table(MRcounts(cumNorm(mp, p=cumNormStat(mp)), norm=TRUE, log=TRUE),taxa_are_rows=TRUE), TAXA, SD,TREE )

#*********sensistive-depleted
cond=c("afull","hs")
physeq.hs=subset_samples(physeq, config2 %in% cond)
comm=c("Back","Comp")
physeq.hs=subset_taxa(physeq.hs, tag %in% comm)

#***shaking
physeq.hssk=subset_samples(physeq.hs, config5 %in% subSk)

sdtab_hssk=data.table(as(sample_data(physeq.hssk), "data.frame"),keep.rownames=T, key="rn")
otumat_hssk=t(as(otu_table(physeq.hssk), "matrix") )

#comuting Bray Curtis or Weighted UniFrac distances
#DBC.hssk=vegdist(otumat_hssk, method="bray")
DwUF.hssk <-UniFrac(physeq.hssk, weighted=TRUE)
adB.hssk=with(sdtab_hssk, adonis(DwUF.hssk ~ config2))

tmp<-data.frame(adB.hssk$aov.tab$R2)
colnames(tmp)<-"13hs.shaking"
x1<-as.data.frame(adB.hssk$aov.tab[6])
p1=x1[1,1]
tmp= rbind(tmp, p1)
rownames(tmp)<-c("explained","residuals","total", "pvalue")
hs.shaking=tmp


#***standing
physeq.hsst=subset_samples(physeq.hs, config5 %in% subSt)

sdtab_hsst=data.table(as(sample_data(physeq.hsst), "data.frame"),keep.rownames=T, key="rn")
otumat_hsst=t(as(otu_table(physeq.hsst), "matrix") )

#comuting Bray Curtis or Weighted UniFrac distances
#DBC.hsst=vegdist(otumat_hsst, method="bray")
DwUF.hsst <-UniFrac(physeq.hsst, weighted=TRUE)
adB.hsst=with(sdtab_hsst, adonis(DwUF.hsst ~ config2))

tmp<-data.frame(adB.hsst$aov.tab$R2)
colnames(tmp)<-"13hs.standing"
x1<-as.data.frame(adB.hsst$aov.tab[6])
p1=x1[1,1]
tmp= rbind(tmp, p1)
rownames(tmp)<-c("explained","residuals","total", "pvalue")
hs.standing=tmp


#*********competitive-depleted
cond=c("afull","hc")
physeq.hc=subset_samples(physeq, config2 %in% cond)
comm=c("Back","Sens")
physeq.hc=subset_taxa(physeq.hc, tag %in% comm)

#***shaking
physeq.hcsk=subset_samples(physeq.hc, config5 %in% subSk)

sdtab_hcsk=data.table(as(sample_data(physeq.hcsk), "data.frame"),keep.rownames=T, key="rn")
otumat_hcsk=t(as(otu_table(physeq.hcsk), "matrix") )

#comuting weighted UniFrac distances
DwUF.hcsk <-UniFrac(physeq.hcsk, weighted=TRUE)
adB.hcsk=with(sdtab_hcsk, adonis(DwUF.hcsk ~ config2))

tmp<-data.frame(adB.hcsk$aov.tab$R2)
colnames(tmp)<-"13hc.shaking"
x3<-as.data.frame(adB.hcsk$aov.tab[6])
p3=x3[1,1]
tmp= rbind(tmp, p3)
rownames(tmp)<-c("explained","residuals","total", "pvalue")
hc.shaking=tmp

#***standing
physeq.hcst=subset_samples(physeq.hc, config5 %in% subSt)

sdtab_hcst=data.table(as(sample_data(physeq.hcst), "data.frame"),keep.rownames=T, key="rn")
otumat_hcst=t(as(otu_table(physeq.hcst), "matrix") )

#comuting weighted UniFrac distances
DwUF.hcst <-UniFrac(physeq.hcst, weighted=TRUE)
adB.hcst=with(sdtab_hcst, adonis(DwUF.hcst ~ config2))

tmp<-data.frame(adB.hcst$aov.tab$R2)
colnames(tmp)<-"13hc.standing"
x3<-as.data.frame(adB.hcst$aov.tab[6])
p3=x3[1,1]
tmp= rbind(tmp, p3)
rownames(tmp)<-c("explained","residuals","total", "pvalue")
hc.standing=tmp


#***plots
lst<-list (adB.hcsk,adB.hcst,adB.hssk,adB.hsst)
for (i in lst) {print(i) }

tab<-cbind(hc.shaking,hs.shaking,hc.standing,hs.standing)
tab1=t(tab[-4:-3,])
DT=melt(data.table(tab1, keep.rownames=T, key="rn"))

	p=ggplot(DT, aes (x= rn, y=value, fill=variable))
	neworder<-c("13hc.shaking","13hs.shaking","13hc.standing","13hs.standing")
	p$data$rn<-ordered(p$data$rn, levels= neworder)
	p1=p+theme_bw()+geom_bar(stat="identity")+geom_label(aes(label= value))+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none")+scale_y_continuous(limits=c(0, 1))


gridExtra::grid.arrange(p1,nrow=2)


