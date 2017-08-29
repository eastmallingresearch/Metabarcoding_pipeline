libray(phyloseq)
library(DESeq2)
library(devtools)
load_all("../..//metabarcoding_pipeline/scripts/myfunctions")
library("BiocParallel")
register(MulticoreParam(8))

# filter samples and set location to factor levels
myfiltbiom <- prune_samples(sample_data(mybiom)[[10]]!="duplicate",mybiom)
myfiltbiom <- prune_samples(sample_data(myfiltbiom)[[1]]!="C",myfiltbiom)
myfiltbiom@sam_data$location <- as.factor(myfiltbiom@sam_data$meters)
colnames(sample_data(myfiltbiom))[c(1,6,11)] <- c("Sample","Distance","Orchard")
levels(sample_data(myfiltbiom)[[1]]) <- c("C","Aisle","Tree")

design=~Orchard + Sample + Orchard:location + Orchard:Sample
dds <- phylo_to_des(myfiltbiom,design=design,parallel=T,fit=T)
contrast=c("Sample","Aisle", "Tree" ) # main effect
# contrast=list(c("conditionN","orchardCider.conditionN","orchardDessert.conditionN"),
# c("conditionY","orchardCider.conditionY","orchardDessert.conditionY"))
contrast=c("Orchard","Cider","Dessert") # orchard effect
contrast=list( "OrchardCider.SampleAisle","OrchardCider.SampleTree") # Grass:Tree (cider)
contrast=list( "OrchardDessert.SampleAisle","OrchardDessert.SampleTree") # Grass:Tree (Dessert)

alpha <- 0.05
res <- results(dds,contrast=contrast,alpha=alpha,parallel=T)
res.merge <- merge(as.data.frame(res),tax_table(myfiltbiom),by="row.names",all.x=TRUE)
rownames(res.merge) <- res.merge$Row.names
res.merge <- res.merge[-1]
sig.res <- subset(res.merge,padj<=alpha)

write.table(merge(as.data.frame(res_orchard),tax_table(myfiltbiom),by="row.names",all.x=TRUE),
"orchard.txt", sep="\t", quote=F,na="",row.names=F)
write.table(merge(as.data.frame(res_main),tax_table(myfiltbiom),by="row.names",all.x=TRUE),
"main.txt", sep="\t", quote=F,na="",row.names=F)
write.table(merge(as.data.frame(res_cider),tax_table(myfiltbiom),by="row.names",all.x=TRUE),
"cider.txt", sep="\t", quote=F,na="",row.names=F)
write.table(merge(as.data.frame(res_dessert),tax_table(myfiltbiom),by="row.names",all.x=TRUE),
"dessert.txt", sep="\t", quote=F,na="",row.names=F)

# gets combined counts at 0.65 conf for given level (1=kingdom,2=phylum,3=class,4=order,5=family,6=genus,7=species)
sig_res <- function(res) {
  res.merge <- merge(as.data.frame(res),tax_table(myfiltbiom),by="row.names",all.x=TRUE)
  rownames(res.merge) <- res.merge$Row.names
  res.merge <- res.merge[-1]
  return(subset(res.merge,padj<=alpha))
}  

sig.res <- sig_res(res_orchard_fun)
lclass <- countTaxa(taxaConf(sig.res[,7:21],0.65,3),"rank")
colnames(lclass)[2] <- "orchard"
lorder <- countTaxa(taxaConf(sig.res[7:21],0.65,4),"rank")
colnames(lorder)[2] <- "orchard"
lfamily <- countTaxa(taxaConf(sig.res[7:21],0.65,5),"rank") 
colnames(lfamily)[2] <- "orchard"

sig.res <- sig_res(res_main_fun)
lclass <- full_join(lclass,countTaxa(taxaConf(sig.res[sig.res$log2FoldChange>0,7:21],0.65,3),"rank"),by="Taxa")
colnames(lclass)[3] <- "main_up"
lclass <- full_join(lclass,countTaxa(taxaConf(sig.res[sig.res$log2FoldChange<0,7:21],0.65,3),"rank"),by="Taxa") 
colnames(lclass)[4] <- "main_down"
lorder <- full_join(lorder,countTaxa(taxaConf(sig.res[sig.res$log2FoldChange>0,7:21],0.65,4),"rank"),by="Taxa")
colnames(lorder)[3] <- "main_up"
lorder <- full_join(lorder,countTaxa(taxaConf(sig.res[sig.res$log2FoldChange<0,7:21],0.65,4),"rank"),by="Taxa")
colnames(lorder)[4] <- "main_down"
lfamily <- full_join(lfamily,countTaxa(taxaConf(sig.res[sig.res$log2FoldChange>0,7:21],0.65,5),"rank"),by="Taxa")
colnames(lfamily)[3] <- "main_up"
lfamily <- full_join(lfamily,countTaxa(taxaConf(sig.res[sig.res$log2FoldChange<0,7:21],0.65,5),"rank"),by="Taxa")
colnames(lfamily)[4] <- "main_down"

sig.res <- sig_res(res_cider_fun)
lclass <- full_join(lclass,countTaxa(taxaConf(sig.res[7:21],0.65,3),"rank"))
colnames(lclass)[5] <- "cider"
lorder <- full_join(lorder,countTaxa(taxaConf(sig.res[7:21],0.65,4),"rank"))
colnames(lorder)[5] <- "cider"
lfamily <- full_join(lfamily,countTaxa(taxaConf(sig.res[7:21],0.65,5),"rank")) 
colnames(lfamily)[5] <- "cider"

sig.res <- sig_res(res_dessert_fun)
lclass <- full_join(lclass,countTaxa(taxaConf(sig.res[7:21],0.65,3),"rank"))
colnames(lclass)[6] <- "dessert"
lorder <- full_join(lorder,countTaxa(taxaConf(sig.res[7:21],0.65,4),"rank"))
colnames(lorder)[6] <- "dessert"
lfamily <- full_join(lfamily,countTaxa(taxaConf(sig.res[7:21],0.65,5),"rank")) 
colnames(lfamily)[6] <- "dessert"

write.table(lclass,"sig.class.count.txt", sep="\t", quote=F,na="0",row.names=F)
write.table(lorder,"sig.order.count.txt", sep="\t", quote=F,na="0",row.names=F)
write.table(lfamily,"sig.family.count.txt", sep="\t", quote=F,na="0",row.names=F)

#### PLOTS

# base means vs log fold change

