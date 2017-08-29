library(phyloseq)
library(DESeq2)
library(gtable)
library(gridExtra)
library(devtools)
load_all("../..//metabarcoding_pipeline/scripts/myfunctions")

# filter out duplicates and control samples
mybiom<-biomITS
mybiom<-biom16

myfiltbiom <- prune_samples(sample_data(mybiom)[[10]]!="duplicate",mybiom)
colnames(sample_data(myfiltbiom))[c(1,6,11)] <- c("Sample","Distance","Orchard")
levels(sample_data(myfiltbiom)[[1]]) <- c("C","Aisle","Tree")
myfiltbiom<-prune_samples(sample_data(myfiltbiom)[[1]]!="C",myfiltbiom)

# myfiltbiom@otu_table@.Data <- round(counts(phylo_to_des(myfiltbiom),normalize=T),0) ## this is wrong, v. bad for Chao1 index

myfiltbiom@otu_table@.Data <- counts(dds,normalize=T)
myfiltbiom@otu_table@.Data[myfiltbiom@otu_table@.Data==0] <- NA
myfiltbiom@otu_table@.Data[myfiltbiom@otu_table@.Data<1] <- 1
myfiltbiom@otu_table@.Data[is.na(myfiltbiom@otu_table@.Data)] <- 0
myfiltbiom@otu_table@.Data <- round(myfiltbiom@otu_table@.Data,0)

myfiltbiom<-prune_samples(colSums(otu_table(myfiltbiom))>999,myfiltbiom)
myfiltbiom<-prune_samples(colSums(otu_table(myfiltbiom))<200001,myfiltbiom)
sample_data(myfiltbiom)$Class <- paste(sample_data(myfiltbiom)$Orchard,sample_data(myfiltbiom)$Sample,sep=" ")
sample_data(myfiltbiom)$Class[sample_data(myfiltbiom)$Class=="Cider Aisle"] <- "C-G"
sample_data(myfiltbiom)$Class[sample_data(myfiltbiom)$Class=="Cider Tree"] <- "C-T"
sample_data(myfiltbiom)$Class[sample_data(myfiltbiom)$Class=="Dessert Aisle"] <- "D-G"
sample_data(myfiltbiom)$Class[sample_data(myfiltbiom)$Class=="Dessert Tree"] <- "D-T"


# Get alpha data
all_alpha <- plot_richness(myfiltbiom,returnData=T)

# ANOVA of alpha data
data.frame(prop.table(summary(aov(Chao1~location+(Sample*Orchard),all_alpha))[[1]][c(2)])*100,
           summary(aov(Chao1~location+(Sample*Orchard),all_alpha))[[1]][5])
data.frame(prop.table(summary(aov(Shannon~location+(Sample*Orchard),all_alpha))[[1]][c(2)])*100,
           summary(aov(Shannon~location+(Sample*Orchard),all_alpha))[[1]][5])
data.frame(prop.table(summary(aov(Simpson~location+(Sample*Orchard),all_alpha))[[1]][c(2)])*100,
           summary(aov(Simpson~location+(Sample*Orchard),all_alpha))[[1]][5])

summary(aov(Chao1~location+(Sample*Orchard),all_alpha))[[1]][2]
summary(aov(Shannon~location+(Sample*Orchard),all_alpha))[[1]][2]
summary(aov(Simpson~location+(Sample*Orchard),all_alpha))[[1]][2]
#~condition*orchard+location

g1 <- plot_richness(myfiltbiom,x="Class",color="Distance",measures=c("Chao1", "Shannon", "Simpson"))
g2 <- plot_richness(myfiltbiom,x="Class",color="Distance",measures=c("Chao1", "Shannon", "Simpson"))

title.A <- textGrob(label = "A",x = unit(0, "lines"),y = unit(0, "lines"),hjust = -0.5, vjust = 0,gp = gpar(fontsize = 16))
title.B <- textGrob(label = "B",x = unit(0, "lines"),y = unit(0, "lines"),hjust = -0.5, vjust = -1,gp = gpar(fontsize = 16))

#dev.off()
g1 <- g1 + theme(legend.position="none",
                 axis.title.x=element_blank(),
                 axis.title.y=element_blank(),
                 plot.margin=unit(c(-1,1,2.5,0.5), "lines")
                )
g2 <- g2 + theme(legend.direction="horizontal",
                 legend.position="bottom",
                 legend.justification=c(0,0),
                 legend.box="vertical",
                 legend.box.just="left",
                 axis.title.x=element_blank(),
                 axis.title.y=element_blank(),
                 plot.margin=unit(c(-2,1,0.5,0.5), "lines")
                )
g3 <- arrangeGrob(g1, top = title.A)
g4 <- arrangeGrob(g2, top = title.B)

pdf("Alpha.pdf", height=8,width=8)
lay=rbind(c(1,1),c(2,2))
grid.arrange(g3,g4,layout_matrix=lay)
dev.off()
