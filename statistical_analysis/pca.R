#===============================================================================
#       Load libraries
#===============================================================================

library(DESeq2)
library("BiocParallel")
register(MulticoreParam(12))
library(data.table)
library(dplyr)
library(plyr)
library(devtools)
load_all("~/pipelines/metabarcoding/scripts/myfunctions")

#===============================================================================
#       Load data 
#===============================================================================

# create a simulated OTU table
countData <- otuSIM(samples=10,treatments=2,otu=100)

ubiom_FUN <- list(
  # load denoised otu count table
  countData=read.table("FUN_otu_table.txt",header=T,sep="\t",row.names=1,comment.char = ""),
  # load sample metadata
	colData=read.table("colData",header=T,sep=",",row.names=1,colClasses=c("character","factor","factor")),
  # load taxonomy data
	taxData=phyloTaxaTidy(read.table("FUN.taxa",header=F,sep=",",row.names=1)[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)],0.8),
  # ubiom name
	RHB="FUN"
)
# add DESeq2 object to Ubiom list
ubiom_FUN$dds <- ubiom_to_des(ubiom_FUN)

# attach objects
invisible(mapply(assign, names(ubiom_FUN), ubiom_FUN, MoreArgs=list(envir = globalenv())))

#===============================================================================
#       Filter Data and perform PCA decompossion
#===============================================================================

# simple filter to select OTUs with >5 counts actross all samples
myfilter <- rowSums(counts(dds,normalize=T))>5

# filter out low abundance OTUs
dds <- dds[myfilter,]

# perform PC decompossion on DES object
mypca <- des_to_pca(dds)

#===============================================================================
#      Anova 
#===============================================================================

# write to file
sink(paste(RHB,"PCA.txt"))

cat("
# ANOVA of first four PC scores \n")
apply(mypca$x[,1:4],2,function(x) summary(aov(x~Treatment,as.data.frame(dds@colData))))

# get sum of squares for all PC scores
sum_squares <- t(apply(mypca$x,2,function(x)t(summary(aov(x~Treatment,as.data.frame(dds@colData)))[[1]][2])))

# name sum_squares columns
colnames(sum_squares) <- c("Treatment","residual")

# proportion of total sum of squares for PC
perVar <- t(apply(sum_squares,1,prop.table)) * mypca$percentVar

cat("
# % total variance explained by the aov model and residual \n")
colSums(perVar)/sum(colSums(perVar))*100

# end write to file
sink()

#===============================================================================
#      Plot PCA
#===============================================================================

# to get pca plot axis into the same scale create a dataframe of PC scores multiplied by their variance
df  <- t(data.frame(t(mypca$x)*mypca$percentVar))

# output pdf
pdf(paste(RHB,"PCA.pdf",sep="_"))

# plot PC1 vs PC2
plotOrd(df,dds@colData,design="Year",shape="Site",xlabel="PC1",ylabel="PC2",labels=T,textSize=8)

# plot PC2 vs PC3
plotOrd(df,dds@colData,design="Year",shape="Site",xlabel="PC2",ylabel="PC3",dimx="PC2",dimy="PC3")




# write to file
dev.off()
