#R scripts used in Hassani et al., to reproduce the analysis of CTU that significantly increased or decreased in rel. abund. upon pertubation (example of FlowPots)
#M. Amine Hassani - ahassani@bot.uni-kiel.de

#required packages for the analysis
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "ape", "GUniFrac", "ape", "phytools", "metagenomeSeq", "ggtern")
lapply(pkg, library, character.only = TRUE)

#clean
rm(list = ls()[grep("*.*", ls())])

#define variable
orderFam<-c("Flavobacteriaceae","Paenibacillaceae","Bacillaceae","Nocardioidaceae","Nocardiaceae","Mycobacteriaceae",
"Streptomycetaceae","Intrasporangiaceae","Cellulomonadaceae","Microbacteriaceae","Promicromonosporaceae","Micrococcaceae",
"Oxalobacteraceae","Comamonadaceae","Alcaligenaceae","Xanthomonadaceae","Moraxellaceae","Pseudomonadaceae",
"Sphingomonadaceae","Caulobacteraceae","Phyllobacteriaceae","Rhizobiaceae","Hyphomicrobiaceae","Bradyrhizobiaceae","Methylobacteriaceae")
myPalette <-c("#EF5656","#47B3DA","#F7A415","#2BB065")
myPalette1 <-c("#EF5656","#F7A415","#8d877e","#2BB065")
myPalette2 <-c("#EF5656","#8d877e","#2BB065")
shape1=c(21,19)

#setting threshold value for log2 fold change in the relative abundance of CTU
folch=1
#setting threshold value for p-value (adjusted for multiple testing)
alpha=0.05

#***upload and prepare phyloseq objects***
mat=read.table( "count_table_fps.txt", sep="\t", row.names=1, header=T)
OTU=otu_table(mat, taxa_are_rows=T) 
tax=read.table("taxonomy_table.txt", sep="\t", row.names=1, header=1)
tax=as.matrix(tax)
TAXA=tax_table(tax)
sd=read.table("sample_data_fps.txt", sep="\t", row.names=1, header=1)
SD=sample_data(sd)
TREE=read.tree("16s_trim_v5v7_KG_DepExp_Derep.tree")
physeq= phyloseq(OTU, TAXA, SD,TREE ) 

#remove samples with less than 1000 reads
#physeq = prune_samples(sample_sums(physeq) > 1000, physeq)

#subset samples to root and matrix samples
cond=c("matrix","root")
physeq=subset_samples(physeq, config1 %in% cond)

#***subset to full and highly competitive depleted samples
hcdep=c("afull","hc")
physeq.hs=subset_samples(physeq, config2 %in% hcdep)
comm=c("Back","Sens")
physeq.hs=subset_taxa(physeq.hs, tag %in% comm)

#***subet to matrix samples
dep=c("HCDepletedMatrix")
condmat=c("matrix")
physeq.hsmat=subset_samples(physeq.hs, config1 %in% condmat)

#testing for differential abundance CTU
m = as(otu_table(physeq.hsmat), "matrix") + 1L
t = data.frame(as(tax_table(physeq.hsmat), "matrix"))
T=AnnotatedDataFrame(t)
s = as(sample_data(physeq.hsmat), "data.frame")
S =AnnotatedDataFrame(s)
obj = newMRexperiment(m,phenoData=S,featureData=T) 
p=cumNormStatFast(obj)
objTrim=cumNorm(obj, p=p)
config3 = pData(obj)$config3
settings = zigControl(maxit = 20, verbose = TRUE)
dsg1=model.matrix(~0+config3, data =s)

#using fonction fitZig to test differential abundance 
res1 = fitZig(obj = objTrim, mod = dsg1, control = settings)
zigFit1 = res1$fit
finalMod1 = res1$fit$design

#making contrast for comparaison of abundance
c.mat1 = makeContrasts ( config3matrix.full	-	config3matrix.hc, levels = finalMod1)
fit1 = contrasts.fit(zigFit1, c.mat1)
fit1 = eBayes(fit1)
DT_1=fit1$coefficients
DT_1p=fit1$p.value

#storing output in data table
DT_1=data.table(DT_1, keep.rownames=T, key="rn")
DT_1p=data.table(DT_1p, keep.rownames=T, key="rn")
setnames(DT_1, "config3matrix.full - config3matrix.hc", "fc.FvsC")
setnames(DT_1p, "config3matrix.full - config3matrix.hc", "p.FvsC")
DT_1=merge(DT_1,DT_1p)
DT_tax=data.table(t, keep.rownames=T, key="rn")
com1=merge(DT_tax, DT_1)

#applying p-value and log2 fold change thresholds
com1$ptest <- ifelse(com1$p.FvsC < alpha, TRUE, F)
com1$fctest <- ifelse(com1$fc.FvsC > folch, 1, ifelse(com1$fc.FvsC< -folch, 1, 0))
com1$p.fc.test <- com1$ptest*com1$fctest

com1$shape.fc <- ifelse(com1$fc.FvsC > 0, "full", "other")
com1$ph <- ifelse(com1$Phylum == "Actinobacteria", 1, ifelse(com1$Phylum == "Bacteroidetes", 2, 
ifelse(com1$Phylum == "Firmicutes", 3,4)))

com1$ph.p.fc.test<- com1$ph*com1$p.fc.test

com1$ph.p.fc.test <- ifelse(com1$ph.p.fc.test ==0,"Non_Significant", ifelse( com1$ph.p.fc.test==1, "Actinobacteria", ifelse(com1$ph.p.fc.test ==2, "Bacteroidetes", ifelse(com1$ph.p.fc.test== 3, "Firmicutes", "Proteobacteria"))))

#preparing the plots
pC=ggplot(com1,aes(-log2(p.FvsC),Family))
pC$data$Family<-ordered(pC$data$Family, levels=orderFam)	

p2=pC+geom_jitter(aes(colour=ph.p.fc.test,size= abs(com1$fc.FvsC),shape=shape.fc))+theme_bw()+
theme(axis.text.x = element_text(angle=90, vjust=1), axis.text.y=element_text(size=4),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
geom_vline(xintercept = -log2(0.05))+
scale_colour_manual(values=myPalette2)+scale_shape_manual(values=shape1)+
scale_size(range = c(1,8), breaks=waiver())+
ggtitle(paste0("Competitive-depleted-",condmat))

p3=p2+geom_text(aes(label=rn, colour=ph.p.fc.test))

#saving plot as pdf files
#pdf(paste0("Competitive-depleted-",cond,".pdf"), useDingbats=FALSE,paper="A4")	 
print(p2)
dev.new()
#pdf(paste0("Competitive-depleted-",condmat,".pdf"), useDingbats=FALSE,width = 13, height = 13)
print(p3)
#dev.off()

#save CTU that signficantly increased or decreased in relative abundance
tmp=com1[com1$ph.p.fc.test != 'Non_Significant']
tmp=tmp[,rn,fc.FvsC]
rownames(tmp)<-tmp$rn
colnames(tmp)<-c(paste0(dep), "rn")
write.table(tmp,paste0(dep,".txt"), sep="\t")
