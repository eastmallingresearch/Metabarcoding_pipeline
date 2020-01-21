# Statistical analyses

There are many ways to analyse metabarcoding data. I've included here some of the more common options. All methods are written in R.

## Metafuncs R package

This package contains a lots of small functions to enable manipulation of metabarcoding data, and some plotting functions

Install the R scripts
```R
library(devtools)
install_github("eastmallingresearch/Metabarcoding_pipeline/scripts")
library(metafuncs)
```

TODO:
Write a package help file (some of the functions already have there own R help)

## Required libraries
```R
library(DESeq2)
library(data.table)
library(tidyverse)
library(Biostrings)
library(vegan)
library(lmPerm)
library(phyloseq)
library(metacoder)
library(ape)
library(metafuncs)
library(viridis)
library(BiocParallel)
```

## Load data
```R
# load the data into R using loadData (metafuncs) function
ubiome_BAC <- loadData("BAC.otu_table.txt","colData","BAC.utax.taxa",RHB="BAC")
ubiome_FUN <- loadData("FUN.otu_table.txt","colData","FUN.utax.taxa",RHB="FUN")

# It can be useful to combine the data for OTUs which are probably from the same species
# combine species at 0.95 (default) confidence (if they are species) 
invisible(mapply(assign, names(ubiome_FUN), ubiome_FUN, MoreArgs=list(envir = globalenv())))
combinedTaxa <- combineTaxa("FUN.utax.taxa")
countData <- combCounts(combinedTaxa,countData)
taxData <- combTaxa(combinedTaxa,taxData)
ubiome_FUN$countData <- countData
ubiome_FUN$taxData <- taxData

```

## Rarefaction curves

This is an updated verion of metafuncs plotRarefaction (need to update in the package...)  

```R
plotRarefaction <- function(countData) {        

  # descending order each sample 
  DT <- data.table(apply(countData,2,sort,decreasing=T))
  
  # get cummulative sum of each sample
  DT <- cumsum(DT)    
  
  # log the count values                            
  DT <- log10(DT)
  
  # set values larger than maximum for each column to NA
  DT <- data.table(apply(DT,2,function(x) {x[(which.max(x)+1):length(x)]<- NA;x}))
  
  # remove rows with all NA
  DT <- DT[rowSums(is.na(DT)) != ncol(DT), ]
  
  # add a count column to the data table
  DT$x <- seq(1,nrow(DT))
  
  # melt the data table for easy plotting 
  MDT <- melt(DT,id.vars="x")
  
  # create an empty ggplot object from the data table
  g <- ggplot(data=MDT,aes(x=x,y=value,colour=variable))
  
  # remove plot background and etc.
  g <- g + theme_classic_thin() %+replace% theme(legend.position="none")
  
  # plot cumulative reads
  g <- g + geom_line(size=1.5) + scale_colour_viridis(discrete=T)
  
  # add axis lables
  g <- g + ylab(expression("Log"[10]*" aligned sequenecs"))+xlab("OTU count")
  
  # Return the plot
  g
}
```

then run

```R
plotRarefaction(ubiome_BAC$countData)
plotRarefaction(ubiome_FUN$countData)
```

## Library size normalisation

```r
# Make DESeq objects for library size normalisation - filter based on rarefaction results
ubiome_FUN$dds <- ubiom_to_des(ubiome_FUN,filter=expression(colSums(countData)>=1000))
ubiome_BAC$dds <- ubiom_to_des(ubiome_BAC,filter=expression(colSums(countData)>=1000))

# Assign ubiome objects to global environment (either/or)
# for fungi
invisible(mapply(assign, names(ubiome_FUN), ubiome_FUN, MoreArgs=list(envir = globalenv())))
# for bacteria
invisible(mapply(assign, names(ubiome_BAC), ubiome_BAC, MoreArgs=list(envir = globalenv())))

```

## Alpha diversity

### Box plots of some measures
```R
plot_alpha(counts(dds,normalize=T),colData(dds),design="Status",colour=NULL,measures=c("Chao1", "Shannon", "Simpson","Observed"),type="box")		   
```

### Permutation ANOVA
```R
# get the diversity index data
all_alpha_ord <- plot_alpha(counts(dds,normalize=T),colData(dds),design="Treatment",returnData=T)

# join diversity indices and metadata
all_alpha_ord <- as.data.table(left_join(all_alpha_ord,colData,by=c("Samples"="Samples")))

# Chao1
setkey(all_alpha_ord,S.chao1)
summary(aovp(as.numeric(as.factor(all_alpha_ord$S.chao1))~Block + Treatment + Genotype + Treatment * Genotype,all_alpha_ord))

# Shannon
setkey(all_alpha_ord,shannon)
summary(aovp(as.numeric(as.factor(all_alpha_ord$shannon))~Block + Treatment + Genotype + Treatment * Genotype,all_alpha_ord))

# Simpson
setkey(all_alpha_ord,simpson)
summary(aovp(as.numeric(as.factor(all_alpha_ord$simpson))~Block + Treatment + Genotype + Treatment * Genotype,all_alpha_ord))
```

## Beta diveristy

Probably best to filter low count OTUs before doing beta diversity analyses

This is a pretty conservative filter
```R
dds <- dds[rowSums(counts(dds, normalize=T))>4,]
```

Most of the analyses area going to require a model. 

### PCA analysis

The count data should be transformed to make it homoskadistic - e.g. DESeq VST transformation

#### PCA plot
```R

# PC decomposition of variance stabilised reads reads
mypca <- des_to_pca(dds)

# to get pca plot axis into the same scale create a dataframe of PC scores multiplied by their variance
d <-t(data.frame(t(mypca$x)*mypca$percentVar))

# plot the PCA (plotOrd has a number of features - and a help file (wow))
plotOrd(d,colData(dds), design="Treatment",alpha=0.75,cbPalette=T,axes=c(1,2))

```

#### PCA ANOVA and sum of squares
```R
# statistical analysis of the first n PCs
lapply(1:4,function(i){
  summary(aov(mypca$x[,i]~ Block + Treatment,data=colData(dds)))
})

# Combined model sum of squares for all PCs
sum_squares <- apply(mypca$x,2,function(x) 
  summary(aov(x ~ Block + Treatment,data=colData(dds)))[[1]][2]
)
sum_squares <- do.call(cbind,sum_squares)
x<-t(apply(sum_squares,2,prop.table))
perVar <- x * mypca$percentVar

# sum squares values
colSums(perVar)

# sum squares %
colSums(perVar)/sum(colSums(perVar))*100
```

### ADONIS and MRPP

```R
# calculate bray-crtis distance matrix
vg <- vegdist(t(counts(dds,normalize=T)),method="bray")

# for replicable results seed needs to be set  
set.seed(sum(utf8ToInt("Greg Deakin")))

# run ADONIS analysis
(fm1 <- adonis(vg ~ Block + Treatment,colData(dds),permutations = 1000))

# run MRPP analysis
(fm2 <- mrpp(vg, colData$Treatment,permutations = 1000))
```

### Make phyloseq object

Phyloseq contains several useful functions for ordination analyses (they're actually wrapper scripts aroung Vegan functions)

```R
myphylo <- ubiom_to_phylo(list(counts(dds,normalize=T),as.data.frame(colData(dds)),taxData))
```

### NMDS

```R
# nmds ordination
ord_rda <- phyloseq::ordinate(myphylo,method="NMDS",distance="bray",formula= ~ Block + Treatment)		

otus <- scores(ord_rda,"species")
nmds <- scores(ord_rda)

# make plot
g <- plotOrd(nmds,colData(dds),design=Treatment,alpha=0.75,cbPalette=T)

# plot
g

# add OTUs to plot
g + geom_point(data=as.data.frame(otus),inherit.aes = F,aes(x=NMDS1,y=NMDS2))

# add some taxonomy arrows to the plot
taxmerge <-data.table(inner_join(data.table(OTU=rownames(otus),as.data.frame(otus)),data.table(OTU=rownames(taxData),taxData)))
taxmerge$phy <- taxaConfVec(taxData[,-8],conf=0.9,level=which(colnames(taxData)=="phylum"))
taxmerge$cls <- taxaConfVec(taxData[,-8],conf=0.9,level=which(colnames(taxData)=="class"))

phylum <- taxmerge[,lapply(.SD,mean),by=phy,.SDcols=c("NMDS1","NMDS2")]
cls <- taxmerge[,lapply(.SD,mean),by=cls,.SDcols=c("NMDS1","NMDS2")]

g + geom_segment(inherit.aes = F,data=cls,aes(xend=NMDS1,yend=NMDS2,x=0,y=0),size=1.5,arrow=arrow()) +
  geom_text(inherit.aes = F,data=cls,aes(x=NMDS1,y=(NMDS2+sign(NMDS2)*0.05),label=cls))  

```

### CCA/RDA and etc.
```R
# transform data using vst
otu_table(myphylo) <-  otu_table(assay(varianceStabilizingTransformation(dds)),taxa_are_rows=T)

# RDA
 ord_rda <- ordinate(myphylo,method="RDA","samples",formula= ~ Block + Treatment)		
anova.cca(ord_rda,permuations=1000)
 anova.cca(ord_rda,permuations=1000,by="terms")

# partial CCA

 ord_cca_partial <- ordinate(myphylo,method="CCA","samples",formula= ~Condition(Block) + Treatment )
 anova.cca(ord_cca_partial,permuations=1000)
 anova.cca(ord_cca_partial,permuations=1000,by="terms")

# plot 
	
#scores scaled by variation in each axes
sscores <- function(ord,axes=c(1,2)) {
	d <- scores(ord,axes)$sites 
	eigvec = eigenvals(ord)
	fracvar = eigvec[axes]/sum(eigvec)
	percVar = fracvar * 100
	d <- t(t(d)*percVar)
	d
}

plotOrd(sscores(ord_rda),colData,design="Treatment",shape="Genotype",pointSize=1.5,alpha=0.75)
plotOrd(sscores(ord_ccaa_partial),colData,design="Treatment",shape="Genotype",pointSize=1.5,alpha=0.75)

```

## Differential OTUs

Currently DESeq2 is the method of choice for differential OTU analysis

```R
# p value for FDR cutoff
alpha <- 0.05

# the model
design <- ~Block + Treatment
 
# add design to dds object
design(dds) <- design

# run model
dds <- DESeq(dds)

# Extract results (this can get complex - see examples in ?results)
res <- results(dds,alpha=alpha,contrast=c("Treatment","Treated","Control"))

```

It is also posible to combine OTU counts at various taxonomic levels and run a differential analysis (there are some fairly good reasons why this mught not always be a good idea)

The metafuncs combineTaxa function will produce a dataframe of every entry at a given taxonomic rank, and a list of OTUs for each entry.  
Actually it doesn't (would do if you could specify level) - but the combineByTaxa function does

```R

combineByTaxa <- 
function (taxData,countData, rank = "species", confidence = 0.95, column_order = -8) 
{
  require(plyr)
  require(data.table)
  
  # reorder taxonomy table (remove rank column by default)
  taxData <- taxData[, column_order]
  
  # get rank level
  level <- which(c("kingdom", "phylum", "class", "order", 
                   "family", "genus", "species")==rank)
  
  # select taxonomy at given confidence
  taxData <- phyloTaxaTidy(taxData, confidence,level=level)

  # convert to data table
  TD <- data.table(taxData, keep.rownames = "OTU")
  
  # get list of OTUs in each entry at the given rank
  combinedOTUs <- ddply(TD, ~rank, summarize, OTUS = list(as.character(OTU)))
  
  # combine the OTUs into list format
  combinedOTUs <- combinedOTUs[lapply(TD[, 2], function(x) length(unlist(x))) > 
            1, ]
  
  list(countData = combCounts(combinedOTUs,ubiome_FUN$countData),
       taxData   = combTaxa(combinedOTUs,taxData))

}


# get the combined OTU list
mylist <- combineByTaxa(taxData,countData,rank="family",confidence=0.8)

# combine the counts for each OTU
familyCounts <- mylist[[1]]

# produce a new taxonomy table for the combined counts (this is not strictly necessary, but useful)
familyTaxa <- mylist[[2]]

# create a new DESeq object
dds2 <- DESeqDataSetFromMatrix(familyCounts,colData,~1)
# these should be similar to original; it's possibly better to copy them from the original dds object...
sizeFactors(dds2) <- sizeFactors(estimateSizeFactors(dds2)) 

# then run as above

```


## Functional annotation 
```R

```

## Network analyis
