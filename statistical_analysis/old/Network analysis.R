#===============================================================================
#       Load libraries
#===============================================================================

library(phyloseq)
library(devtools)
load_all("metabarcoding_pipeline/scripts/myfunctions") # this is a set of scipts to do various things (plotOrd, plotPCA, geoSet)
library(data.table)
library(dplyr)
library(future)
plan(multiprocess)
library(EcoSimR)
library(vegan)
library(cooccur)

#===============================================================================
#       Load data 
#===============================================================================

##### 16S #####

biom_file = "analysis/16S.taxa.biom"
colData = "analysis/colData"
mybiom <- import_biom(biom_file) 
sample_data(mybiom) <- read.table(colData,header=T,sep="\t",row.names=1)
tax_table(mybiom) <- phyloTaxaTidy(tax_table(mybiom),0.65)
biom16<-mybiom

##### ITS #####

biom_file = "analysis/ITS.taxa.biom"
colData = "analysis/colData"
mybiom <- import_biom(biom_file) 
sample_data(mybiom) <- read.table(colData,header=T,sep="\t",row.names=1)
tax_table(mybiom) <- phyloTaxaTidy(tax_table(mybiom),0.65)
biomITS<-mybiom

mybioms <- list(bacteria=biom16,fungi=biomITS)


#===============================================================================
#       co-oocurance analysis with cooccur2- I've hacked cooccur to run a bit faster
#===============================================================================

myfiltbiom <- prune_samples(sample_data(biomITS)$site=="Chestnuts",biomITS)
cotable <- as.data.frame(as.matrix(otu_table(myfiltbiom)))

cotable_h <- cotable[,row.names(sample_data(prune_samples(sample_data(myfiltbiom)$condition=="Healthy",myfiltbiom)))]
cotable_h <- cotable_h[rowSums(cotable_h)>5,colSums(cotable_h)>5]
cotable_h[cotable_h>0] <- 1  
CHcoHmodel <- cooccur2(cotable_h,type = "spp_site",spp_names = T,thresh = T)

cotable_s <- cotable[,row.names(sample_data(prune_samples(sample_data(myfiltbiom)$condition=="Symptom",myfiltbiom)))]
cotable_s <- cotable_s[rowSums(cotable_s)>5,colSums(cotable_s)>5]
cotable_s[cotable_s>0] <- 1 
CHcoSmodel <- cooccur2(cotable_s,type = "spp_site",spp_names = T,thresh = T)

cotable <- cotable[rowSums(cotable)>5,colSums(cotable)>5]
cotable[cotable>0] <- 1
CHcomodel <-  cooccur2(cotable,type = "spp_site",spp_names = T,thresh = T)

CHcoHmodel$results$padj <- p.adjust(apply(CHcoHmodel$results[,8:9],1, min),"BH")
CHcoSmodel$results$padj <- p.adjust(apply(CHcoSmodel$results[,8:9],1, min),"BH")
CHcomodel$results$padj <- p.adjust(apply(CHcomodel$results[,8:9],1, min),"BH")



myfiltbiom <- prune_samples(sample_data(biomITS)$site=="Bigwood",biomITS)
cotable <- as.data.frame(as.matrix(otu_table(myfiltbiom)))

cotable_h <- cotable[,row.names(sample_data(prune_samples(sample_data(myfiltbiom)$condition=="Control",myfiltbiom)))]
cotable_h <- cotable_h[rowSums(cotable_h)>5,colSums(cotable_h)>5]
cotable_h[cotable_h>0] <- 1  
BWcoHmodel <- cooccur2(cotable_h,type = "spp_site",spp_names = T,thresh = T)

cotable_s <- cotable[,row.names(sample_data(prune_samples(sample_data(myfiltbiom)$condition=="Symptom",myfiltbiom)))]
cotable_s <- cotable_s[rowSums(cotable_s)>5,colSums(cotable_s)>5]
cotable_s[cotable_s>0] <- 1 
BWcoSmodel <- cooccur2(cotable_s,type = "spp_site",spp_names = T,thresh = T)

cotable <- cotable[rowSums(cotable)>5,colSums(cotable)>5]
cotable[cotable>0] <- 1
BWcomodel <-  cooccur2(cotable,type = "spp_site",spp_names = T,thresh = T)

BWcoHmodel$results$padj <- p.adjust(apply(BWcoHmodel$results[,8:9],1, min),"BH")
BWcoSmodel$results$padj <- p.adjust(apply(BWcoSmodel$results[,8:9],1, min),"BH")
BWcomodel$results$padj <- p.adjust(apply(BWcomodel$results[,8:9],1, min),"BH")




# cotable_h <- cotable_h[rowSums(cotable_h)>5,colSums(cotable_h)>5]
# cotable_s <- cotable_s[rowSums(cotable_s)>5,colSums(cotable_s)>5]

#fth1 <- invisible(future({cooccur(cotable_h,type = "spp_site",spp_names = T,thresh = T)}))
#fts1 <- invisible(future({cooccur(cotable_s,type = "spp_site",spp_names = T,thresh = T)}))
#ftsh2 <- future({cooccur(cotable_h,type = "spp_site",spp_names = T,thresh = F)})



# multiple testing correction

#coHmodel$results$padj_lt <-  p.adjust(coHmodel$results$p_lt,"BH")
#coHmodel$results$padj_gt <-  p.adjust(coHmodel$results$p_gt,"BH")
#coSmodel$results$padj_lt <-  p.adjust(coSmodel$results$p_lt,"BH")
#coSmodel$results$padj_gt <-  p.adjust(coSmodel$results$p_gt,"BH")


cotable_h <- cotable_h[rowSums(cotable_h)>5,colSums(cotable_h)>5]
cotable_s <- cotable_s[rowSums(cotable_s)>5,colSums(cotable_s)>5]

cotable_all <- cotable[rowSums(cotable)>5,colSums(cotable)>5]
BWft <- <- future({cooccur(cotable_all,type = "spp_site",spp_names = T,thresh = T)})
BWcomodel <- value(BWft)

BWfth1 <- future({cooccur(cotable_h,type = "spp_site",spp_names = T,thresh = T)})
BWfts1 <- future({cooccur(cotable_s,type = "spp_site",spp_names = T,thresh = T)})


#===============================================================================
#       co-oocurance analysis (EcoSimR)
#===============================================================================



fth <- future({cooc_null_model(cotable_h,nReps=10000,burn_in = 1000,suppressProg=T)})
fts <- future({cooc_null_model(cotable_s,nReps=10000,burn_in = 1000,suppressProg=T)})

myHModel <- value(fth)
mySModel <- value(fts)

summary(myHModel)
summary(mySModel)

plot(myHModel,type="hist")
plot(mySModel,type="hist")

plot(myHModel,type="cooc")
plot(mySModel,type="cooc")

plot(myHModel,type="burn_in")
plot(mySModel,type="burn_in")

dev.off()

# weights<- rowMeans(otu_table(mybiom))/sum(otu_table(mybiom))
# weights<- rowSums(otu_table(mybiom))/sum(otu_table(mybiom))


 

#===============================================================================
#       installing gmp (for cooccur)
#===============================================================================
```shell
mkdir gmp_build
cd gmp_build
wget https://gmplib.org/download/gmp/gmp-6.1.2.tar.bz2
tar --bzlib2 -xvf gmp-6.1.2.tar.bz2
cd gmp-6.1.2.tar.bz2
./configure --prefix=~/usr/local/ --disable-shared --enable-fat
make
make check
make install
```
```R
### I haven't checked this - I set the CPPFLAGS and LDFLAGS, PKG_CPPFLAGS and PKG_LIBS flags manually by hacking the configure script in the package
install.packages("gmp",lib.loc="~/usr/local")
```
