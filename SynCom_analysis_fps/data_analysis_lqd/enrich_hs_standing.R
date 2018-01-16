	#***require packages for the analysis the analysis
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "ape", "GUniFrac", "ape", "phytools", "metagenomeSeq", "ggtern")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

orderFam<-c("Flavobacteriaceae","Paenibacillaceae","Bacillaceae","Nocardioidaceae","Nocardiaceae","Mycobacteriaceae",
"Streptomycetaceae","Intrasporangiaceae","Cellulomonadaceae","Microbacteriaceae","Promicromonosporaceae","Micrococcaceae",
"Oxalobacteraceae","Comamonadaceae","Alcaligenaceae","Xanthomonadaceae","Moraxellaceae","Pseudomonadaceae",
"Sphingomonadaceae","Caulobacteraceae","Phyllobacteriaceae","Rhizobiaceae","Hyphomicrobiaceae","Bradyrhizobiaceae","Methylobacteriaceae")

myPalette <-c("#EF5656","#47B3DA","#F7A415","#2BB065")
myPalette1 <-c("#EF5656","#F7A415","#8d877e","#2BB065")
myPalette2 <-c("#EF5656","#8d877e","#2BB065")

shape1=c(21,19)

#log2 fold change
folch=1
#adjusted p-value
alpha=0.05

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

#subset to full, 13HC and 13HS depleted samples
cond=c("afull","hc", "hs")
physeq=subset_samples(physeq, config2 %in% cond)

#*********sensistive-depleted
cond=c("afull","hs")
physeq.hs=subset_samples(physeq, config2 %in% cond)
comm=c("Back","Comp")
physeq.hs=subset_taxa(physeq.hs, tag %in% comm)

#***standing
dep=c("HSDepletedStanding")
submat=c("strtsb")
physeq.hsst=subset_samples(physeq.hs, config5 %in% submat)


#***differentially abundant CTU 
m = as(otu_table(physeq.hsst), "matrix") + 1L
t = data.frame(as(tax_table(physeq.hsst), "matrix"))
T=AnnotatedDataFrame(t)
s = as(sample_data(physeq.hsst), "data.frame")
S =AnnotatedDataFrame(s)
obj = newMRexperiment(m,phenoData=S,featureData=T) 

p=cumNormStatFast(obj)
objTrim=cumNorm(obj, p=p)

config2 = pData(obj)$config2

settings = zigControl(maxit = 20, verbose = TRUE)

dsg1=model.matrix(~0+config2, data =s)

res1 = fitZig(obj = objTrim, mod = dsg1, control = settings)

zigFit1 = res1$fit

finalMod1 = res1$fit$design

c.mat1 = makeContrasts ( config2afull	-	config2hs, levels = finalMod1)

fit1 = contrasts.fit(zigFit1, c.mat1)

fit1 = eBayes(fit1)

DT_1=fit1$coefficients

DT_1p=fit1$p.value

	DT_1=data.table(DT_1, keep.rownames=T, key="rn")
	DT_1p=data.table(DT_1p, keep.rownames=T, key="rn")

	setnames(DT_1, "config2afull - config2hs", "fc.FvsS")
	setnames(DT_1p, "config2afull - config2hs", "p.FvsS")
	

	DT_1=merge(DT_1,DT_1p)
	DT_tax=data.table(t, keep.rownames=T, key="rn")
	com1=merge(DT_tax, DT_1)

com1$ptest <- ifelse(com1$p.FvsS < alpha, TRUE, F)
com1$fctest <- ifelse(com1$fc.FvsS > folch, 1, ifelse(com1$fc.FvsS< -folch, 1, 0))
com1$p.fc.test <- com1$ptest*com1$fctest
com1$shape.fc <- ifelse(com1$fc.FvsS > 0, "full", "other")

com1$ph <- ifelse(com1$Phylum == "Actinobacteria", 1, ifelse(com1$Phylum == "Bacteroidetes", 2, 
ifelse(com1$Phylum == "Firmicutes", 3,4)))

com1$ph.p.fc.test<- com1$ph*com1$p.fc.test

com1$ph.p.fc.test <- ifelse(com1$ph.p.fc.test ==0,"Non_Significant", ifelse( com1$ph.p.fc.test==1, "Actinobacteria", ifelse(com1$ph.p.fc.test ==2, "Bacteroidetes", ifelse(com1$ph.p.fc.test== 3, "Firmicutes", "Proteobacteria"))))

	#com2$rn <- factor(com2$rn, levels=rev(unique(neworder)))
	#com2$Family <- factor(com2$Family, levels=rev(unique(orderFam)))

	pC=ggplot(com1,aes(-log2(p.FvsS),Family))
	pC$data$Family<-ordered(pC$data$Family, levels=orderFam)	

	p2=pC+geom_jitter(aes(colour=ph.p.fc.test,size= abs(com1$fc.FvsS),shape=shape.fc))+theme_bw()+
	theme(axis.text.x = element_text(angle=90, vjust=1), axis.text.y=element_text(size=4),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
	geom_vline(xintercept = -log2(0.05))+
	scale_colour_manual(values=myPalette1)+scale_shape_manual(values=shape1)+
	scale_size(range = c(1,8), breaks=waiver())+
	ggtitle(paste0("Sensitive-depleted-",submat))

	p3=p2+geom_text(aes(label=rn, colour=ph.p.fc.test))

#pdf(paste0("Sensitive-depleted-",cond,".pdf"), useDingbats=FALSE,paper="A4")	 
p2
#dev.off()
#pdf(paste0("Sensitive-depleted-",submat,".pdf"), useDingbats=FALSE,width = 13, height = 13)
print(p3)
#dev.off()
#Save CTU signficantly increased or decreased in relative abundance
tmp=com1[com1$ph.p.fc.test != 'Non_Significant']
tmp=tmp[,rn,fc.FvsS]
rownames(tmp)<-tmp$rn
colnames(tmp)<-c(paste0(dep), "rn")
#write.table(tmp,paste0(dep,".txt"), sep="\t")



















