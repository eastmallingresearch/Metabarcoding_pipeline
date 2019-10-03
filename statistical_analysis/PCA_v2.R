#===============================================================================
#       Load libraries
#===============================================================================
library(data.table)
library(tidyverse)
library(metafuncs)

#===============================================================================
#       Load data 
#===============================================================================

dat  <- fread("myData.txt")
meta <- freaa("myMetaData.txt")

# or create some random data
dat  <- data.table(var1 = c(rnorm(20,10,3),rnorm(20,12,3),rnorm(20,8,2.5)),
                   var2 = c(rnorm(20,32,4),rnorm(20,30,2),rnorm(20,21,5)),
                   var3 = c(rnorm(20,50,4),rnorm(20,55,4),rnorm(20,20,6)),
                   var4 = c(runif(20,1,40),runif(20,1,40),runif(20,1,40)))

meta <- data.table(Sample=1:60,
                   Treatment=c(rep("treat1",20),rep("treat2",20),rep("treat3",20)),
                   Site=c("site1", "site2"))

#===============================================================================
#       Filter Data and perform PCA decompossion
#===============================================================================

mypca <- prcomp(dat,scale=T)

mypca$percentVar <- mypca$sdev^2/sum(mypca$sdev^2)

#===============================================================================
#      Anova 
#===============================================================================

apply(mypca$x[,1:4],2,function(x) summary(aov(x~Treatment,as.data.frame(dds@colData))))

# get sum of squares for all PC scores
sum_squares <- t(apply(mypca$x,2,function(x)t(summary(aov(x~Treatment,as.data.frame(dds@colData)))[[1]][2])))

# name sum_squares columns
colnames(sum_squares) <- c("Treatment","residual")

# proportion of total sum of squares for PC
perVar <- t(apply(sum_squares,1,prop.table)) * mypca$percentVar

# % total variance explained by the aov model and residual \n")
colSums(perVar)/sum(colSums(perVar))*100

# end write to file

#===============================================================================
#      Plot PCA
#===============================================================================

# to get pca plot axis into the same scale create a dataframe of PC scores multiplied by their variance
df  <- t(data.frame(t(mypca$x)*mypca$percentVar))

# output pdf
pdf(paste(RHB,"PCA.pdf",sep="_"))

# plot PC1 vs PC2
plotOrd(df,dds@colData,design="Year",shape="Site",xlabel="PC1",ylabel="PC2")

# write to file
dev.off()
