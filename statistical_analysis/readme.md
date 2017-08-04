# Statistical analysis

## R analysis
R will import table with character columns as factors. This might be prefereable if your doing lots of linear modeling, but for any data manipulation it is a disaster and can lead to unexpected (or worse unnoticed) errors.

Set stringAsFactor to false - apparently this can be set in .Rprofile, but this doesn't work on my production cluster environment
```{r}
options(stringsAsFactors = FALSE)
```
Changed my mind on this - it's not so difficult to work with factors (they're still a bit of a pain though)

### Load libraries
```{r}
library(DESeq2)
library(phyloseq)
library(ape)
library(vegan)
library(ncf)
library(data.table)

#There are several R functions throughout this module which are in the myfunctions package
library(devtools)
load_all("../..//metabarcoding_pipeline/scripts/myfunctions")
```

### phyloseq
Import biom table, taxonomy and sample data into a phyloseq object
Create phyloseq object in R

The taxonomy imported from UPARSE doesn't include rank names and adds extra stuff to the taxa names (k__,p__ and etc.)
phyloTaxaTidy will fix this and change "unknown" to the lowest known rank (and append a character indicating rank)
```{r}
biom_file = "16S.taxa.biom"
otu_file = "16S.otus.fa" # might be useful at some stage
colData = "colData"
mybiom <- import_biom(biom_file) # ,refseqfilename=out_file
sample_data(mybiom) <- read.table(colData,header=T,sep="\t",row.names=1)
tax_table(mybiom) <- tax_table(phyloTaxaTidy(tax_table(mybiom),0.65))

# An example of removing certain OTUs from a phyloseq object - and using S4 methods to access the data.
# This will filter based on OTU present in the first column "condition" of colData 
t1 <- aggregate(t(otu_table(mybiom)),by=list(sample_data(mybiom)[[1]]),FUN=sum)[-1]
myfiltbiom <- prune_taxa(apply(t1,2,prod)>0,mybiom)
```
### Modify coldata factor levels
```R
levels(mybiom@sam_data$condition)[levels(mybiom@sam_data$condition)=="something"] <- "something_else"
```

##### Merge ITS1 and ITS2
```{r}
ITS1 = "ITS1.taxa.biom"
ITS2 = "ITS2.taxa.biom"
colData = "colData"
mybiom <- merge_phyloseq(import_biom(ITS1),import_biom(ITS2))
sample_data(mybiom) <- read.table(colData,header=T,sep="\t",row.names=1)
tax_table(mybiom) <- tax_table(phyloTaxaTidy(tax_table(mybiom),0.65))
```

##### Create and add phylogentic tree to mybiom - need to add step for creating phy object
```shell
# doesn't work in usearch9...
usearch8 -calc_distmx 16S.otus.fa -distmxout 16S.phy  -threads 12 -format phylip_square
```


```{r}
library(ape)
temp_connection = file("16S.phy", 'r')
len = readLines(temp_connection, n=1)
len = as.numeric(len)
close(temp_connection)
phylip_data = read.table("16S.phy", fill=T, row.names=1, skip=1, col.names=1:len)
nj.16S <- nj(as.dist(phylip_data))
write.tree(nj.16S,"16S.tree")
phy_tree(mybiom) <- nj.16S
```

#### Core biom
```{r}
library(DESeq2)
# normalise the reads
mynormbiom <- mybiom
mynormbiom@otu_table@.Data <- counts(phylo_to_des(mynormbiom,fitType="local"),normalized=T)
## one of my samples had a massively over inflated count for a single OTU. Replacing calcFactors with a new geoMean calculation to ignore 0 values fixes this..
#mynormbiom@otu_table@.Data <- counts(phylo_to_des(mynormbiom,fitType="local",calcFactors=function(d)geoMeans(d),normalized=T))

cond="Y"
myfiltbiom <- mynormbiom
myfiltbiom <- prune_samples((sample_data(mynormbiom)[[10]]=="experiment"),mynormbiom)
myfiltbiom <- prune_samples(sample_data(myfiltbiom)[[1]]==cond,myfiltbiom)

otu_prop_table <- t(t(otu_table(myfiltbiom))/colSums(otu_table(myfiltbiom)))

min_freq <- 0.002   # the minimum count frequency (per sample) for OTU to be considered present (using 0.001 for fungi)
min_samp <- 0.8  # the minimum proportion of samples for OTU to be present ot be include in core biom 
mycorebiom <- prune_taxa(apply(otu_prop_table,1,function(x) (sum(x>=min_freq))/ncol(otu_prop_table)>=min_samp),myfiltbiom)
summary(colSums(otu_table(otu_prop_table)[rownames(otu_table(mycorebiom)),]))

write.table(tax_table(mycorebiom),"core.taxa",quote=F,sep="\t")
write.table(round(otu_table(mycorebiom),0),"core.otu",quote=F,sep="\t")
write.table(round(otu_table(otu_prop_table)[rownames(otu_table(mycorebiom)),],4),"core.prop",sep="\t",quote=F)

## plotting with plotTaxa will require converting back to unnormalised reads or supplying a seperate transform function
otu_table(mycorebiom) <- otu_table(mybiom)[rownames(otu_table(mycorebiom)),]
# produces a graph of proportions within the core microbiome
plotTaxa(mycorebiom,"genus","condition",type=2, others=F,fitType="local",ordered=T)
```

### Spatial analysis

#### Anova of PCA eigenvectors

```{r}
# filter out OTUs with less than 6 counts (and remove "extra samples")
myfiltbiom <- prune_samples(sample_data(mybiom)[[10]]!="duplicate",mybiom)
myfiltbiom@sam_data$location <- as.factor(myfiltbiom@sam_data$meters)

# plotPCA will return a prcomp object if returnData is set to TRUE
myfiltbiom <- prune_taxa(rowSums(otu_table(myfiltbiom))>5,myfiltbiom)
mypca <- plotPCA(myfiltbiom,design="1",ntop= nrow(myfiltbiom@otu_table),returnData=T,fitType="local",blind=T)

# for data from multiple sequencer runs (control samples named "internal" in column 10 of colData)
mypca <- plotPCA(
   myfiltbiom,
   design="1",
   ntop= nrow(myfiltbiom@otu_table),
   returnData=T,
   fitType="local",
   blind=T,
   filterFun=function(o,f){
    m <- prune_samples((sample_data(o)[[1]]!="C"),o);
    m <- prune_taxa(rowSums(otu_table(m))>5,m)
    return(m)
   },
   calcFactors=function(d,o){
     obj<-des_to_phylo(d);
     control_samples<-rownames(colData(d)[colData(d)[[10]]=="internal",]);
     normHTS(obj,control_samples)
   }
)
myfiltbiom <- prune_samples(sample_data(myfiltbiom)[[1]]!="C",myfiltbiom)
myfiltbiom <- prune_taxa(rowSums(otu_table(myfiltbiom))>5,myfiltbiom)


# get the sum of squares for tree/aisle, location and residual
sum_squares <- t(apply(mypca$x,2,function(x) t(summary(aov(x~condition+location,sample_data(myfiltbiom)))[[1]][2])))
colnames(sum_squares) <- c("condition","location","residual")
perVar <- sum_squares * mypca$percentVar
colSums(perVar)
colSums(perVar)/sum(colSums(perVar))*100
```
Plot residual after removing spatial component for first couple of PCA vectors
```{r}
# uncorrected
df <- data.frame(mypca$x[,1]*mypca$percentVar[1],mypca$x[,2]*mypca$percentVar[2])
plotOrd(df,sample_data(myfiltbiom))
# spatial removed
pc.res <- resid(aov(mypca$x~sample_data(myfiltbiom)$location))
d <- data.frame(pc.res[,1]*mypca$percentVar[1],pc.res[,2]*mypca$percentVar[2])
plotOrd(d,sample_data(myfiltbiom))

mypca$y <- t(t(mypca$x)*mypca$percentVar)
plotOrd(mypca$y,sample_data(myfiltbiom),dimx=1,dimy=2)

```
Manova/Anova of first couple of pca with tree/aisle
```{r}
fit <- manova(mypca$x[,1:3]~condition,as.data.frame(as.matrix(sample_data(myfiltbiom))))
summary(fit, test="Pillai") # could just call summary directly 


summary(aov(mypca$x[,1]~(condition*field)+location,as.data.frame(as.matrix(sample_data(myfiltbiom)))))
summary(aov(mypca$x[,2]~(condition*field)+location,as.data.frame(as.matrix(sample_data(myfiltbiom)))))
summary(aov(mypca$x[,3]~(condition*field)+location,as.data.frame(as.matrix(sample_data(myfiltbiom)))))


```

#### Alpha diversity measures
```R
myfiltbiom <- prune_samples(sample_data(mybiom)[[10]]!="duplicate",mybiom)
myfiltbiom<-prune_samples(sample_data(myfiltbiom)[[1]]!="C",myfiltbiom)
all_alpha <- plot_richness(myfiltbiom,returnData=T)
summary(aov(Chao1~condition*orchard+location,all_alpha))[[1]][2]
summary(aov(Shannon~condition*orchard+location,all_alpha))[[1]][2]
summary(aov(Simpson~condition*orchard+location,all_alpha))[[1]][2]
```
#### Auto correlation

```{r}
library(ape)
library(vegan)
library(ncf)
library(data.table)

myfiltbiom@sam_data$gap <- 0

cond <- "Y"

pc.x <- scores(mypca)[sample_data(myfiltbiom)$condition==cond,]
col.x <- sample_data(myfiltbiom)[sample_data(myfiltbiom)$condition==cond,]
pc.dt<- data.table(merge(pc.x,col.x,by="row.names"))
pc.reshape <- dcast(pc.dt,distance~.,fun.aggregate=function(x) mean(x,na.rm=T),value.var=c(names(pc.dt)[grep("PC",names(pc.dt))]))
names(pc.reshape)[grep("PC",names(pc.reshape))] <- sub("_.*","",names(pc.reshape)[grep("PC",names(pc.reshape))])
col.reshape <- sample_data(col.x)[sample_data(col.x)$replicate=="a"]
col.reshape <- col.reshape[order(col.reshape$meters)]

# Moran I test

distmat <- as.matrix(dist(cbind(sample_data(col.reshape)$distance, rep(0,24))))
distmat.inv <- 1/distmat
distmat.inv[is.infinite(distmat.inv)] <- 0
moran <- apply(pc.reshape,2,function(x) t(Moran.I(x,distmat.inv)))
names(moran)->temp
moran <- do.call(rbind,moran)
rownames(moran) <- temp
# this is similar to PCNM method shown below
distmat.bin <- (distmat > 0 $ distmat <=7.2)
moran.bin <- apply(pc.x,2,function(x) t(Moran.I(x,distmat.bin)))
rownames(moran.bin) <- rownames(moran)

# Moran correlogram

moran.mv  <- lapply(seq(1,10),function(y) correlog(sample_data(col.reshape)$distance,sample_data(col.reshape)$gap,pc.reshape$PC1,increment=y,quiet=T))
sapply(seq(1,10),function(x) plot.correlog(moran.mv[[x]]))
# I've knocked up a ggplot2 alternative plotting function, gets rid of the box around the plots
# second argument is the (two-tail) sig figure to colour points black
lapply(seq(1,10),function(x) plot.corr(moran.mv[[x]][c(1:3,5)],0.025))
dev.off()

#  Pearson Correlogram
cutoff <- 17
pc<-"PC1"
plotCorrelog(mypca,myfiltbiom,pc,cutoff=cutoff,xlim=NULL,ylim=c(-1,1),na.add=c(9,17))
dev.off()

### For H samples - due to experimental design
t1 <- plotCorrelog(mypca,prune_samples(sample_data(myfiltbiom)$block!=3,myfiltbiom),pc,na.add=c(9),returnCD=T)
t2 <- plotCorrelog(mypca,prune_samples(sample_data(myfiltbiom)$block==3,myfiltbiom),pc,returnCD=T)
t3 <- sapply(1:nrow(t1),function(i) if(i<=nrow(t2)){cbind(rbind(t1[[i,1]],t2[[i,1]]),rbind(t1[[i,2]],t2[[i,2]]))}else{cbind(t1[[i,1]],t1[[i,2]])})
d <- as.data.frame(t(sapply(1:length(t3),function(i) diag(cor(t3[[i]],use="pairwise.complete.obs")[c(1,3),c(2,4)]))))
d$V3 <- as.numeric(t1[[3]])

plotCorrelog(data=d,cutoff=cutoff,pc=pc,ylim=c(-1,1))
dev.off()
```


### DESeq2 differential OTU abundance
Using DESeq2 it's possible to calculate the probability of OTUs having different abundances between condtions. The default will use the the condition column of the dds object's colData table, and take the first two conditions. To specify a different column or use different "condtions use the contrast=c("column_name","condition_1","condition_2") construct when calling the results method.

It's possible to convert a phyloseq object to a DESeq datamatrix with the wrapper function phylo_to_des.
phylo_to_des has the option to specify the size factor calculation using the option calcFactors (see plotTaxa for examples of how to use this option). Set fit=T to fit a GLM model to the data. Further arguments will be passed to DESeq.
```{r}
# Can be a bit slow if lots of samples (100+), best to use multiple cores 
library("BiocParallel")
register(MulticoreParam(8))

myfiltbiom <- prune_samples((sample_data(mybiom)[[10]]=="experiment"&(sample_data(mybiom)[[1]]!="C")),mybiom)
myfiltbiom <- prune_taxa(rowSums(otu_table(myfiltbiom))>5,myfiltbiom)
myfiltbiom@sam_data$location <- as.factor( myfiltbiom@sam_data$meters)
dds <- phylo_to_des(myfiltbiom,fit=T, fitType="local",..)
# specify a design formula to look at conditon excluding info on loation (i.e. independent of spatial data) 
dds <- phylo_to_des(myfiltbiom,fit=T, fitType="local",design=~location+condition,parallel=T)
# for dds error all genes contain at least one zero
dds <- phylo_to_des(myfiltbiom,fit=T, fitType="local",design=~location+condition,parallel=T,calcFactors=geoMeans)
# complex design (cider + dessert orchards)
# This will give the difference in tree aisle controling for blocks and orchard (I think)
dds <- phylo_to_des(myfiltbiom,design=~orchard + orchard:location + condition,parallel=T,fit=T)
# This will give the difference in tree aisle controling for blocks at each orchard (I think)
dds2 <- phylo_to_des(myfiltbiom,design=~orchard + orchard:location + orchard:condition)
```
This code calculates the result object and mergse the phyloseq taxonomy object with the results.

```{r}
alpha <- 0.05 # significance level
res = results(dds, alpha=alpha,cooksCutoff=T,parallel=T)
#res <- results(dds,contrast=list(c("conditiongrass", "conditiontree" )),alpha=alpha,cooksCutoff=T,parallel=T)
#res = results(dds, alpha=alpha,contrast=c("condition","N","K"))	## specify different contrasts to make
res.merge <- merge(as.data.frame(res),tax_table(myfiltbiom),by="row.names",all.x=TRUE)
rownames(res.merge) <- res.merge$Row.names
res.merge <- res.merge[-1]
sig.res <- subset(res.merge,padj<=alpha)
```

#### Beta-diversity statistical analysis
Using PERMANOVA (best adonis method to use???)
PCA score are generated from library size normalised and variance stabilised (DESeq2) OTU counts.  
```{r}
myfiltbiom <- prune_samples(sample_data(mybiom)[[1]]=="some_condition",mybiom)
myfiltbiom <- prune_taxa(rowSums(otu_table(myfiltbiom))>5,myfiltbiom)
vld<-assay(varianceStabilizingTransformation(phylo_to_des(myfiltbiom,~1),blind=F)) # blind t/f should give same results
vld <- vld+abs(min(vld)) # add constant (min) to bring all values to 0 or above
adonis(t(vld)~condition,as.data.frame(as.matrix(sample_data(myfiltbiom))),method='bray') # bray is non-parametric
```

###### UniFrac
phyloseq object needs a tree adding (code is above)


### Spatial analyses not implemented

#### CCA
Simple first step - correspondence analysis

1. All data
2. Tree station
3. Aisle

```{r}
myfiltbiom <- prune_samples(sample_data(mybiom)[[10]]=="experiment",mybiom)

f1 <- prune_samples(sample_data(myfiltbiom)[[1]]!="C",myfiltbiom)
f2 <- prune_samples(sample_data(myfiltbiom)[[1]]=="Y",myfiltbiom)
f3 <- prune_samples(sample_data(myfiltbiom)[[1]]=="N",myfiltbiom)

f1 <- prune_taxa(rowSums(otu_table(f1))>2,f1)
f2 <- prune_taxa(rowSums(otu_table(f2))>2,f2)
f3 <- prune_taxa(rowSums(otu_table(f3))>2,f3)

CCA1 <- ordinate(f1,method="CCA",formula=~distance)
CCA2 <- ordinate(f2,method="CCA",formula=~distance)
CCA3 <- ordinate(f3,method="CCA",formula=~distance)

plot(CCA1)
plot(CCA2)
plot(CCA3)
dev.off()

anova(CCA1) # permutation analysis
anova(CCA2)
anova(CCA3)
```

#### PCNM
More advanced technique; Principal Coordinates of Neighbour Matrices
but see Gilbert & Bennett 2010 for probelms with this (and other) spatial analysis
```{R}
library(vegan)
library(packfor)

# apply various filters to data
cond="Y"
myfiltbiom <- prune_samples((sample_data(mybiom)[[10]]=="experiment")&(sample_data(mybiom)[[1]]==cond),mybiom)
myfiltbiom <- prune_taxa(rowSums(otu_table(myfiltbiom))>5,myfiltbiom)

# the phyloseq object structure messes up some of the following scripts, therefore convert to an S3 object
myubiom <- phylo_to_ubiom(myfiltbiom) 
myubiom$colData$gap <- 0
myubiom$colData$distance <- as.numeric(myubiom$colData$distance)
myubiom$colData$genotype <- as.factor(myubiom$colData$genotype)
myubiom$colData$block <- as.factor(myubiom$colData$block)
myubiom$colData <- myubiom$colData[order(myubiom$colData$distance),]
myubiom$countData <- myubiom$countData[,rownames(myubiom$colData)]

#calculate euclidean distance between sample points
euclid <- dist(myubiom$colData[,6:7],method="euclidean")

#rda
#unweighted PCNM for RDA, theshold is set to cover ~ 5 samples which seems to describe the maximum variation 
pcnm1 <- pcnm(euclid,threshold=7.2)
ord <- rda(t(myubiom$countData)~scores(pcnm1))
# test model significance (is the data described by a spatial dimension)
anova(ord)
# if yes find which spatial dimensions are significant
dim(scores(pcnm1)) # 14 eigenvectors
ord <- rda(t(myubiom$countData)~scores(pcnm1)[,1],...,scores(pcnm1)[,14])
anova(ord,by="terms",parellel=12,model="direct",permu=2000)
forward.sel(t(myubiom$countData),scores(pcnm1),Xscale=F,nperm=2000) # this should give results equivelent to anova above without having to redefine the model. 


ord <- rda(t(myubiom$countData)~genotype+block+Condition(scores(pcnm1)[,c(1,2,x,y)]))
anova(ord) # model shouldn't be significant for year 0  data

plot(ord)
msoplot(mso(ord, sample_data(myfiltbiom)[,6:7]))
dev.off()

#cca
#calculate PCNM weights for CCA
cs <- colSums(myubiom$countData/sum(myubiom$countData))
#weighted PCNM (this is more useful for cca)
pcnm2 <- pcnm(euclid,w=cs,threshold=7.2)
#cca on OTU data with positive eigen vectors from PCNM as the independent variables
#the residual should have no distance trend (i.e. the CA component of the ord object)
ord <- cca(t(otu_table(myfiltbiom))~scores(pcnm2))
anova(ord)
plot(ord)
msoplot(mso(ord, sample_data(myfiltbiom)[,6:7]))
dev.off()

#spatial free data
newCounts <- t(ord$CA$Xbar)
## 
ord2 <- cca(t(otu_table(myfiltbiom))~scores(pcnm2)[,1]+scores(pcnm2)[,2]+scores(pcnm2)[,3..n])
anova(ord2)
anova(ord2,by="terms",parellel=12,model="direct",permu=2000)
```

With threshold of 7.2 using RDA, for Goatham both tree and aisle have significant vectors 2,3,8 and 10, with 5 close to significant for both.

Calculate the pearson correlation coefficient for each OTU against the forward selected spatial parameters
Two methods - produce slightly different results 
```{R}
pc.median <- aggregate(pc.x,by=as.list(col.x[,4]),median)
pc.median <- pc.median[order(as.numeric(as.character(pc.median[,1]))),]

#library(Hmisc)
#library(stats)
#test2 <- rcorr(t(myubiom$countData),scores(pcnm2)[,c(14,15,16,1,8)],type="pearson")
#test3$r <- test2$r[1:7084,7085:7089]
#test3$P <- test2$P[1:7084,7085:7089]
#test3$n <- test2$n[1:7084,7085:7089]
#test3$p.adj <- apply(test3$P,2,function(x) p.adjust(x,"BH"))
library("psych")
test4 <- corr.test(t(myubiom$countData),scores(pcnm2)[,c(14,15,16,1,8)],adjust="BH")
```


###[16S workflow](../master/16S%20%20workflow.md)
###[ITS workflow](../master//ITS%20workflow.md)



OLD STUFF (to delete)
```shell
# Moran correlogram
#distmat <- as.matrix(dist(cbind(sample_data(col.x)$distance, sample_data(col.x)$gap)))
#distmat.inv <- 1/distmat
#diag(distmat.inv) <- 0
#distmat.inv[is.infinite(distmat.inv)] <- 0
#moran <- apply(pc.x,2,function(x) t(Moran.I(x,distmat.inv)))
#names(moran)->temp
#moran <- do.call(rbind,moran)
#rownames(moran) <- temp


#  Pearson Correlelog
ct1 <- correr2(pc.reshape$PC1)
ct2 <- correr2(pc.reshape$PC2)

##recalulate and find
ca1 <- correr2(pc.reshape$PC1)
ca2 <- correr2(pc.reshape$PC2)
d<-as.data.frame(cbind(ct1,ct2,ca1,ca2,pc.reshape$distance[1:22]))
names(d) <- c("Tree_PC1","Tree_PC2","Aisle_PC1","Aisle_PC2","Distance")

d2 <- melt(d[1:12,],id="Distance")
colnames(d2)[3] <- c("Correlation")

xlims=NULL
ylims=NULL#c((min(d2$Correlation)-0.1),(max(d2$Correlation)+0.1))
g <- ggplot(d2)
g <- g + coord_fixed(ratio = 10, xlim = xlims, ylim = ylims, expand = TRUE)
g <- g + theme_bw()
g <- g + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
g <- g + theme(axis.line.x = element_line(size=0.5,colour = "black"),axis.line.y = element_line(size=0.5,colour = "black"),axis.text = element_text(colour = "black"))
g <- g + geom_line(na.rm=T,aes(x=Distance, y=Correlation, colour=variable))
g <- g + scale_colour_manual(values=c("red","green","blue","orange"))
#g <- g + geom_point(na.rm=T,size=2.5,mapping=aes())
g


#pdf("bac.tree.correlogs.pdf")
#moran.res <- lapply(seq(1,10),function(y) correlog(sample_data(col.x)$distance,sample_data(col.x)$gap,pc.x[,1],increment=y,quiet=T))
#sapply(seq(1,10),function(x) plot.correlog(moran.res[[x]]))
#plot(correlog(sample_data(col.x)$distance,sample_data(col.x)$gap,pc.x[,1],increment=7.2,quiet=T))
#dev.off()

#pc <- 1
#df <- merge(pc.x[,pc],col.x,by="row.names")
#n.df <- #reshape(df,idvar="location",timevar="rep",drop=c("Row.names","condition","time_point","block","genotype_code"),direction="wide")
#n.df <- reshape(df,idvar="location",timevar="replicate",drop=c("Row.names","condition","time","block","genotype","genotype_name","sample_id","type"),direction="wide")
#n.df <- n.df[order(n.df$meters.a),]
#n.df <- n.df[,c(2,6,10,1,3,4,5)]
#moran.mv  <- lapply(seq(1,10),function(y) correlog(n.df$distance.a,n.df$gap.a,rowMeans(n.df[,1:3],na.rm=T),increment=y,quiet=T))
#sapply(seq(1,10),function(x) plot.correlog(moran.mv[[x]]))
# I've knocked up a ggplot2 alternative plotting function, gets rid of the box around the plots
# second argument is the (two-tail) sig figure to colour points black
#lapply(seq(1,10),function(x) plot.corr(moran.mv[[x]][c(1:3,5)],0.025))
#dev.off()

#moran.mv  <- lapply(seq(1,10),function(y) correlog(n.df$meters.a,rep(0,24),rowMeans(n.df[,1:3]),increment=y,quiet=T,na.rm=T))
#sapply(seq(1,10),function(x) plot.correlog(moran.mv[[x]]))

# Mantel test
mydist <- dist(cbind(sample_data(col.x)$distance, sample_data(col.x)$gap))
mantel.out <- t(apply(pc.x,2,function(x) unlist(mantel(mydist,dist(x),permutations=9999)[3:4])))

# Mantel correlogram
# this was copied from http://www.ats.ucla.edu not certain of the point of removing the spatial data before looking for autocorrelation
pc.res <- resid(aov(pc.x~sample_data(col.x)$location))
pc.dist <- dist(pc.res)
pc.correlog <- mantel.correlog(pc.dist,cbind(sample_data(col.x)$distance,sample_data(col.x)$gap),cutoff=F)
plot(pc.correlog)
```
