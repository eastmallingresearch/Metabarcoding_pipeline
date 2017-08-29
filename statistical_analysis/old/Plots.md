##### Differential dispersion graphs plots
```{r}
with(res.merge,plot(log2FoldChange,log10(baseMean),pch=20, main="Volano like plot", xlim=c(-2.5,2)))
with(subset(res.merge, padj<.01 ), points(log2FoldChange, log10(baseMean), pch=20, col="red"))
with(subset(res.merge, abs(log2FoldChange)>1), points(log2FoldChange, log10(baseMean), pch=20, col="orange"))
with(subset(res.merge, padj<.01 & abs(log2FoldChange)>1), points(log2FoldChange, log10(baseMean), pch=20, col="green"))

rld <- varianceStabilizingTransformation(dds)
plot1 <- prune_taxa(rownames(sig.res[(sig.res$log2FoldChange>0),]),myfiltbiom)
plot1@otu_table@.Data <- assay(rld[rownames(otu_table(plot1))])
plotTaxa(plot1,"family","condition",type=2, others=F,fitType="local",ordered=T,trans=F,proportional=F)
plotTaxa(plot1,"family","condition",type=2, others=T,fitType="local",ordered=T,trans=F,proportional=T)
dev.off()
plot2 <- prune_taxa(rownames(sig.res[(sig.res$log2FoldChange<0),]),myfiltbiom)
```


### plots

#### alpha diversity
Phyloseq has a builtin method for plotting alpha diversity graphs (plot_richness). Estimate_richness will return the data points used to plot the graphs (I may edit this function as the graphs are pretty basic - ugly anyway)

```{r}
# res <- estimate_richness(mybiom) ## data used for plot_richness graphs
pdf("16S.alpha_bysex.pdf", height=8,width=8)
plot_richness(mybiom,x="condition",color="Sex",measures=c("Chao1", "ACE", "Shannon", "Simpson"))
dev.off()
```
#### beta diversity
NMDS etc. plots 

```{R}
obj <- mybiom
obj@otu_table@.Data <- assay(varianceStabilizingTransformation(phylo_to_des(obj)))
mynmds <- ordinate(obj,method = "NMDS",distance="bray",autotransform=F,try=100)
#plotOrd(mynmds$points,sample_data(obj),design="condition",shapes="location")
plot_ordination(obj,mynmds,color="condition",shape="location")
```

#### plotOrd
plotOrd is a ggplot wrapper that does something similar to the phyloseq plotting method, but without the background and with axes on the same scale.  

Example - PCA (PC1 vs PC2) with most abundant genera used to colour samples.
```{R}
myfiltbiom <- mybiom
myfiltbiom <- prune_taxa(rowSums(otu_table(myfiltbiom))>5,myfiltbiom)
mypca <- plotPCA(myfiltbiom,design="1",ntop= nrow(myfiltbiom@otu_table),returnData=T,fitType="local",blind=T)
sample_data(myfiltbiom)$sample <- row.names(sample_data(myfiltbiom))
myubiom <- phylo_to_ubiom(myfiltbiom)
myubiom$countData <- as.data.frame(assay(varianceStabilizingTransformation(ubiom_to_des(myubiom,design=~1),blind=T)))
myubiom$colData$sample <- row.names(myubiom$colData)
mytaxa <- sumTaxa(myubiom,"genus","sample")
mytaxa <- mytaxa[order(rowSums(mytaxa[,-1]),decreasing=T),]
myubiom$colData[mytaxa[1,1]] <- as.numeric(t(mytaxa[1,-1]))
df <- data.frame(mypca$x[,1]*mypca$percentVar[1],mypca$x[,2]*mypca$percentVar[2])
plotOrd(df,myubiom$colData,design=mytaxa[1,1],continuous=T)

row.names(mytaxa) <- mytaxa[,1]
lapply(seq(1:20),function(x) plotOrd(pc.res.var,as.data.frame(t(mytaxa[,-1])),design=as.character(mytaxa[x,1]),continuous=T,xlabel="PC1",ylabel="PC2")))

```
plotOrd has several other features - need to write about them at sometime.

Grid plot example (using combined Heineken and Goatham bioms)
This requires 2 legends one horizonatal and one vertical:
make legened free
get legend 1
get legenend 2
stick them together
```R
library(gridExtra)
library(grid)
g <-plotOrd(df,sample_data(myfiltbiom),shapes=c("Orchard","Sample"),design="Distance" ,xlabel="PC1",ylabel="PC2",continuous=T,ylims=c(-3,3),xlims=c(-6,10)) + theme(legend.position="None")
t <- plotOrd(df,sample_data(myfiltbiom),design="Distance") + theme(legend.direction="horizontal")
l1 <- ggplot_legend(t)
t <- plotOrd(df,sample_data(myfiltbiom),shapes=c("Orchard","Sample")) + theme(legend.direction="vertical")
l2 <- ggplot_legend(t)
#g <- g + geom_text(aes(label = LETTERS[i], x = 8, y = 11), hjust = -1, size=7)+theme(text = element_text(size=14))
g <- ggplot_gtable(ggplot_build(g))
g$layout$clip[g$layout$name == "panel"] <- "off"
ml <- list(g,l1,l2)
grid.arrange(grobs=ml,nrow=3)
```

Filled clusters at 95% with continuous distance scale in greyscale
```R
plotOrd(df,sample_data(myfiltbiom),shapes=c("Orchard","Sample"),cluster=0.95,design="Distance",xlabel="PC1",ylabel="PC2",continuous=T,pallete="greyscale",ylims=c(-3,3),xlims=c(-8,10),labels=F,centers=1)
```


#### plotPCA

plotPCA is a modified version of the DESeq2 version. 
It take the following options:

1. object (DESeq2 - required) a DESeq object 
2. intgroup (string - optional, default="condition") a column of colData used to describe (colour) the samples (e.g. infected/control)
3. labelby (string - optional) a 2nd column of colData used to descibe (shape) the samples (e.g. male/female)
4. ntop (int - optional, default=500) number of OTUs in descending count order to use in the PCA calculation
5. pcx (int - optional, default=1) X-axis pricipal component
6. pcy (int - optional, default=2) Y-axis pricipal component
7. returnData (bool - optional, default=F) returns the PCA results and exits
8. cofix (bool - optional, default=F) produces a graph with axes on the same scale
9. transform(fun - optional) a user supplied function to replace DESeq2 variance stabilising transform to transform the count matrix. A DESeq2 object will be passed to this function. 

plotPCAWithLabels will produce the same graph but with the addition of sample labels - useful for getting the name of outliers 
```{r}
pdf("16S.beta-diversity.pdf",height=8,width=8)
plotPCAWithLabels(dds)
dev.off()
```
NOTE - I  prefer to calculate and keep the PCA scores (using this function), then use plotOrd for plotting.

#### plotTaxa
plotTaxa produces a ggplot2 bar chart of taxa counts
It takes the following options:

1. obj (phyloseq - required) must is a phyloseq object which must include taxonomy and sample data
2. taxon (str - required) is the taxonomic level of interest
3. condition (str - required) describes how the samples should be grouped (must be column of sample data)
4. proportional (bool - optional, def=T) whether the graph should use proportional or absolute values
5. cutoff (double - optional, def =1) for proportional graphs. 
6. topn (int - optional)taxons to display (by total reads) for non-prortional graphs. 
7. others (bool - optional, def=T), combine values less than cutoff/topn into group "other"
8. ordered (bool - optional, =F) order by value (max to min)
9. type (int (1/2) - required, def=1) Type 1 produces a stacked barchart by sample, type 2 a barchart by taxa 
10. fixed (bool - optional, def=F) fixed is a ggplot parameter to apply coord_fixed(ratio = 0.1)
11. ncol (int - optional, def=1) ncol is a ggplot paramter to use n columns for the legend
12. transform (fun - optional) a user supplied function to replace DESeq2 variance stabilising transform to transform the count matrix. An S3 list of countData, Taxonomy and colData is passed to the transform function . 
13. ... arguments to pass to transform function (o.k. these could just be set in the function, but this is a neater solution if the default function is used)

The default transform firts creates a DESeq2 object using ubiom_to_des. If all OTUs contain at least one zero value DDS can't calculate sizeFactors. The transform function could be modified to calculate VST on a data matrix, or a function "calcFactors" can be passed to the ubiom_to_des method which can calculate the sizeFatcors (e.g, edgeR calcNormFactors) 

```{r}
pdf("16S.phylum.pdf",height=8,width=8)
plotTaxa(mybiom,"phylum","condition")
plotTaxa(mybiom,"phylum","condition",blind=F) # passes  blind=F to transform function

# all OTUs contain at least one zer - pass calcFactors to ubiom_to_des
t1 <- aggregate(t(otu_table(mybiom)),by=list(sample_data(mybiom)[[1]]),FUN=sum)[-1]
myfiltbiom <- prune_taxa(apply(t1,2,prod)>==0,mybiom) # Keep OTUs with at least one zero value 
plotTaxa(myfiltbiom,design="condition",calcFactors=function(d){library(edgeR);sizeFactors(d) <- calcNormFactors(counts(d))}) 

# precalculate transform
dds <- phylo_to_des(mybiom,design=~1)
rld <- varianceStabilizingTransformation(dds)
mybiom <- mytransbiom
mytransbiom@otu_table@.Data <- assay(rld) # couldn't get method otu_table(mybiom) to accept rld data..
plotTaxa(mytransbiom,taxon="genus",design="condition",trans=F)
dev.off()
```
