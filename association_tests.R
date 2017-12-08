
# script to reproduce the association tests (phyloGLM) and
# corresponding figure (Fig. 2) from Hassani et al., 2017
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de


options(warn=-1)

# cleanup

rm(list=ls())

# load plotting functions

source("plotting_functions.R")

# load plotting functions

library("phylolm")
library("ape")
library("ggplot2")
library("scales")
library("grid")

# directories

data.dir <- "./"
figures.dir <- "./"

### files
 
# mapping.txt           metadata table for At-SPHERE isolates (from Bai et al., 2015)
# wgs_taxonomy.txt      WGS taxonomic prediction for At-SPHERE isolates
# ABBA.txt              table containing BGCs per strain and antagonistic activity
# ABBA_amphora.newick   newick tree of tested At-SPHERE strains 

mapping.file <- paste(data.dir, "mapping.txt", sep="")
taxonomy.file <- paste(data.dir, "wgs_taxonomy.txt", sep="")
phenotypes.file <- paste(data.dir, "ABBA.txt", sep="")
tree.file <- paste(data.dir, "ABBA_amphora.newick", sep="")

# load data

mapping <- read.table(mapping.file, sep="\t", header=T, check.names=F, colClasses="character")
taxonomy <- read.table(taxonomy.file, sep="\t", header=T)
phenotypes <- read.table(phenotypes.file, sep="\t", header=T, check.names=F)
tree <- read.tree(tree.file)

### PhyloGLM tests

# set False Discovery Rate alpha threshold

alpha <- 0.05

# Phylogenetic logistic regression described in Ives and Garland (2010)

method <- "logistic_MPLE"
# method <- "poisson_GEE"
# method <- "logistic_IG10"

# extract genomic features (BGCs) from master table
# MPLE allows for continuous phe

traits <- as.matrix(phenotypes[, 8:33])

# extract phenotypes (antagonistic activity against target strains)
# MPLE model requires discrete phenotypes

phenotypes <- as.matrix(phenotypes[, 36:232])
phenotypes <- (phenotypes > 0 ) * 1

### fit the model for each phenotype

# initiallize p. value and coefficient matrices

pvals <- matrix(ncol=ncol(phenotypes), nrow=ncol(traits), NA)
coefs <- matrix(ncol=ncol(phenotypes), nrow=ncol(traits), NA)

colnames(pvals) <- colnames(coefs) <- colnames(phenotypes)
rownames(pvals) <- rownames(coefs) <- colnames(traits)

cs <- matrix(ncol=4, nrow=ncol(traits), NA)
colnames(cs) <- c("Estimate", "StdErr", "z.value", "p.value")
rownames(cs) <- colnames(traits)

# iterate over each phenotype and fit the model if not empty

for (j in 1:ncol(phenotypes)) {

    if (sum(phenotypes[, j])!=0) {
    
        for (i in 1:ncol(traits)) {
        
        print(paste(j, "/", ncol(phenotypes), " ",  i, "/", ncol(traits), sep=""))
        
            # fit model and extract summary statistics

            dat <- data.frame(trait=phenotypes[, j], predictor=traits[, i])
            rownames(dat) <- rownames(phenotypes)
            fit <- phyloglm(trait ~ predictor, data=dat, phy=tree, method=method)
        
            s <- summary(fit)
            cs[i, ] <- s$coefficients[2, ]
        
        }
        
        # store p. values and coefficients

        cs <- data.frame(cs)

        pvals[, j] <- cs$p.value
        coefs[, j] <- cs$Estimate

    } else {
    
        pvals[, j] <- 1
        coefs[, j] <- 0
    
    }
}

# apply multiple testing correction

pvals.adj <- apply(pvals, 2, FUN=p.adjust)

# set to 0 coefficients for non-significant traits

coefs.adj <- coefs
coefs.adj[pvals.adj > alpha] <- 0

# take only significant predictors

coefs.adj.sig <- coefs.adj[rowSums(abs(coefs.adj)) > 0, ]
coefs.adj.sig <- coefs.adj.sig[rowSums(abs(coefs.adj.sig) > 0) > 1, ]

