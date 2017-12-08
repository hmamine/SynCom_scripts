
# script to reproduce the association tests (phyloGLM) and
# corresponding figure (Fig. 2) from Hassani et al., 2017
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# subset mapping table and get only phenotypes with significant associations

idx <- colnames(coefs.adj.sig) %in% mapping$ID
colnames(coefs.adj.sig)[idx] <- mapping$Strain[match(colnames(coefs.adj.sig)[idx], mapping$ID)]

### plot heatmap

# create data frame with coefficients

df <- melt(coefs.adj.sig)
colnames(df) <- c("trait", "phenotype", "coef")

colors <- rev(c("#ff0000", "#ff5555", "#ffaaaa",
                "#ffffff",
                "#aaccff", "#5599ff", "#0066ff"))

b <- c(100, 1, 0.5, 0.1, 0, -0.5, -1, -100)

df$coef_bins <- cut(df$coef, breaks=b, right=F)

# perform agglomerative hierarchical clustering on the association
# coefficients (gives the order of the traits) and on the phenotypes
# (gives the order of the columns) using WPGMC

m <- "median"

h <- hclust(dist(coefs.adj.sig), method=m)
df$trait <- factor(df$trait, levels=rev(rownames(coefs.adj.sig)[h$order]))

h <- hclust(dist(t(coefs.adj.sig)), method=m)
df$phenotype <- factor(df$phenotype, levels=colnames(coefs.adj.sig)[h$order])

# saturate outliers to brightest colors

df$coef[df$coef > 2.5] <- 2.5
df$coef[df$coef < -2.5] <- -2.5

# plot heatmap using ggplot

p1 <- ggplot(df, aes(x=phenotype, y=trait, fill=coef)) +
      geom_tile() +
      scale_fill_gradient2(low="blue", mid="white", high="red") +
      labs(x="", y="", title="") + 
      main_theme +
      theme(legend.title=element_blank(),
            legend.position="left",
            axis.text.x=element_text(angle=45, hjust=1, size=10),
            axis.line=element_blank(),
            plot.margin=unit(c(0, 0, 0, 0), "mm"))

ggsave(file=paste(figures.dir, "heatmap.pdf", sep=""), p1, height=5, width=15)

# plot taxonomic assignment

df_x <- aggregate(df$coef, by=list(df$phenotype), FUN=sum)
colnames(df_x) <- c("phenotype", "sum")

df_x$tax <- taxonomy$Phylum[match(df_x$phenotype, taxonomy$isolate_ID)]
df_x$tax <- droplevels(df_x$tax)
colors <- rep(c_black, length(levels(df_x$tax)))
colors[levels(df_x$tax)=="Actinobacteria"] <- c_red
colors[levels(df_x$tax)=="Proteobacteria"] <- c_dark_green
colors[levels(df_x$tax)=="Bacteroidetes"] <- c_blue
colors[levels(df_x$tax)=="Firmicutes"] <- c_yellow

p2 <- ggplot(df_x, aes(x=phenotype, y=sum, fill=sum)) +
      geom_bar(stat="identity") +
      scale_fill_gradient2(low="blue", mid="white", high="red") +
      labs(x="", y="", title="") + 
      main_theme +
      theme(legend.title=element_blank(),
            legend.position="left",
            axis.text.x=element_text(angle=45, hjust=1, size=10),
            axis.line=element_blank(),
            plot.margin=unit(c(0, 0, 0, 0), "mm"))

ggsave(file=paste(figures.dir, "barplot_x.pdf", sep=""), p2, height=2, width=15)

# plot accumulated coefficients

df_y <- aggregate(df$coef, by=list(df$trait), FUN=sum)
colnames(df_y) <- c("trait", "sum")

p2 <- ggplot(df_y, aes(x=trait, y=sum, fill=sum)) +
      geom_bar(stat="identity") +
      scale_fill_gradient2(low="blue", mid="white", high="red") +
      labs(x="", y="", title="") + 
      main_theme +
      coord_flip() +
      theme(legend.title=element_blank(),
            legend.position="left",
            axis.text.x=element_text(angle=45, hjust=1, size=10),
            axis.line=element_blank(),
            plot.margin=unit(c(0, 0, 0, 0), "mm"))

ggsave(file=paste(figures.dir, "barplot_y.pdf", sep=""), p2, height=5, width=3)

# plot accumulated sensitivity for each target strain

sensitivity <- colSums(phenotypes)                  
names(sensitivity) <- mapping$Strain[match(names(sensitivity), mapping$ID)]
sensitivity <- sensitivity[names(sensitivity) %in% df$phenotype]
df_x2 <- data.frame(phenotype=names(sensitivity), sensitivity=sensitivity)
df_x2$phenotype <- factor(df_x2$phenotype, levels=levels(df$phenotype))

p3 <- ggplot(df_x2, aes(x=phenotype, y=sensitivity)) +
      geom_bar(stat="identity") +
      labs(x="", y="", title="") + 
      main_theme +
      theme(legend.title=element_blank(),
            legend.position="left",
            axis.text.x=element_text(angle=45, hjust=1, size=10),
            axis.line=element_blank(),
            plot.margin=unit(c(0, 0, 0, 0), "mm"))

ggsave(file=paste(figures.dir, "barplot_x2.pdf", sep=""), p3, height=2, width=15)

