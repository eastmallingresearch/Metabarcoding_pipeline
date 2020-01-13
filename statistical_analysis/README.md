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

# Make DESeq objects for library size normalisation
ubiome_FUN$dds <- ubiom_to_des(ubiome_FUN,filter=expression(colSums(countData)>=1000))
ubiome_BAC$dds <- ubiom_to_des(ubiome_BAC,filter=expression(colSums(countData)>=1000))

# Assign objects to global environment (either/or)
# for fungi
invisible(mapply(assign, names(ubiome_FUN), ubiome_FUN, MoreArgs=list(envir = globalenv())))
# for bacteria
invisible(mapply(assign, names(ubiome_BAC), ubiome_BAC, MoreArgs=list(envir = globalenv())))

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
plotRarefaction(countData)
```

## Alpha diversity

### Box plots of some measures
```R
plot_alpha(counts(dds,normalize=T),colData(dds),design="Status",colour=NULL,measures=c("Chao1", "Shannon", "Simpson","Observed"),type="box")		   
```

### Statistical analysis (permutation ANOVA)
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


```R

```

### PCA analysis
```R

```
### NMDS
```R

```

### CCA/RDA and etc.
```R

```

## Differential OTUs
```R

```
## Functional annotation 
```R

```

## Network analyis
