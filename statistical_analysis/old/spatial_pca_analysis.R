### load libraries ###
library(DESeq2)
library(phyloseq)
library(data.table)
library(gtable)
library(gridExtra)
library(devtools)
load_all("../..//metabarcoding_pipeline/scripts/myfunctions")

### "Tidy" samples and names (specific to this project - could be done outside R) ###
myfiltbiom <- prune_samples(sample_data(mybiom)[[10]]!="duplicate",mybiom)
myfiltbiom@sam_data$location <- as.factor(myfiltbiom@sam_data$meters)
colnames(sample_data(myfiltbiom))[c(1,6,11)] <- c("Sample","Distance","Orchard")
levels(sample_data(myfiltbiom)[[1]]) <- c("C","Grass","Tree")
myfiltbiom<-prune_samples(colSums(otu_table(myfiltbiom))>999,myfiltbiom)

### Combined data ####

# Filter function for plotPCA
filterFun=function(o,f){
  o<-prune_samples(sample_data(o)[[1]]!="C",o)
  prune_taxa(rowSums(otu_table(o))>5,o)	
}
# Uses DESeq2 to variance stabelise counts on whole dataset then calcultes PCA scores 
# on filtered data
mypca <- plotPCA(myfiltbiom,design="1",returnData=T,filterFun=filterFun)
# Remove control samples
myfiltbiom <- prune_samples(sample_data(myfiltbiom)[[1]]!="C",myfiltbiom)
# Remove low count OTUs
myfiltbiom <- prune_taxa(rowSums(otu_table(myfiltbiom))>5,myfiltbiom)
sample_data(myfiltbiom)$orchloc <- paste(sample_data(myfiltbiom)$orchard,sample_data(myfiltbiom)$location,sep="_")

# Calculate ANOVA for first 4 PCs (the design is unbalanced for this data as we are missing 2 samples)
lapply(seq(1:4),function(x) 
	summary(aov(mypca$x[,x]~location+(Orchard*Sample),data.frame(as.matrix(sample_data(myfiltbiom)))))
)

# Calcultae sum of squares
sum_squares <- t(apply(mypca$x,2,function(x) 
  t(summary(aov(x~location+(Orchard*Sample),data.frame(as.matrix(sample_data(myfiltbiom)))))[[1]][2]))
)
colnames(sum_squares) <- c("location","orchard","condition","orchard:sample","residual")
x<-t(apply(sum_squares,1,prop.table))
perVar <- x * mypca$percentVar
#perVar <- sum_squares * mypca$percentVar
colSums(perVar)
colSums(perVar)/sum(colSums(perVar))*100

# PCA plots
# uncorrected
df <- t(data.frame(t(mypca$x)*mypca$percentVar))	
# location effect removed...
pc.res <- resid(aov(mypca$x~sample_data(myfiltbiom)$location))
d <- t(data.frame(t(pc.res)*mypca$percentVar))
# all spatial effect removed...
pc.res <- resid(aov(mypca$x~sample_data(myfiltbiom)$Orchard+sample_data(myfiltbiom)$location))
d2 <- t(data.frame(t(pc.res)*mypca$percentVar))

# plotting function with clusters at 95% confidence	
myglist <- lapply(list(df,d,d2),function(x)
  plotOrd(x,sample_data(myfiltbiom),
	shapes=c("Orchard","Sample"),
	cluster=0.95,
	design="Distance",
	xlabel="PC1",
	ylabel="PC2",
	continuous=T,
	dimx=1,
	dimy=2,
	colourScale=c("black","lightblue"),
	centers=1,
	ylims=c(-4,4),
	legend=F
))	

myglist_its<- myglist
myglist_16S<- myglist		  
		  
g1 <- myglist_its[[1]]
g2 <- myglist_16S[[1]]	  		  
		  
#g2 <- plotOrd(df,...)
	
g2 <- g2 + theme(legend.direction="horizontal", 
		 legend.position="bottom",
		 legend.justification=c(0,0),
		 legend.box="vertical",
		 legend.box.just="left",
		 axis.line.x = element_line(size=0.3,colour = "black"),
		 axis.line.y = element_line(size=0.3,colour = "black"),
		 axis.text = element_text(colour = "black"),
		 plot.margin=unit(c(2.5,1,0.5,0.5), "lines")
)	

### Main figure		  
pdf("ITS_16S_orchards2.pdf",width=7,height=6.5)	
g3 <- ggplotGrob(g1+annotate("text",label=paste("A"), x=20, y=3,size=5))
g4<-  ggplotGrob(g2+annotate("text",label=paste("B"), x=60, y=3,size=5))
g <- rbind(g3, g4, size="first")
grid.newpage()
grid.draw(g)
dev.off()			  
		  
		  
# get the two legends
l1 <- ggplot_legend(plotOrd(df,sample_data(myfiltbiom),design="Distance",continuous=T,colourScale=c("black","lightblue"))+ theme(legend.direction="horizontal"))
l2 <- ggplot_legend(plotOrd(df,sample_data(myfiltbiom),shapes=c("Orchard","Sample")) + theme(legend.direction="horizontal"))

	
# print the 3 graphs using grid.arrange	
lay=rbind(c(1,1),c(1,1),c(2,2),c(2,2),c(3,3),c(3,3),c(4,5))  
pdf("ITS_all_PCA.pdf",width=8,height=6)	
grid.arrange(
	myglist[[1]]+annotate("text",label=paste("A"), x=50, y=3,size=5),
	myglist[[2]]+annotate("text",label=paste("B"), x=50, y=3,size=5),
	myglist[[3]]+annotate("text",label=paste("C"), x=20, y=3,size=5),
	layout_matrix=lay
)
	


	
myglist <- append(myglist,myglist16S)
	,c(7,7,7,7)
lay=rbind(c(1,1,2,2),c(1,1,2,2),c(1,1,2,2),
	  c(3,3,4,4),c(3,3,4,4),c(3,3,4,4),
	  c(5,5,6,6),c(5,5,6,6),c(5,5,6,6),
	  c(7,8,8,NA),c(7,8,8,NA),c(9,NA,NA,NA)) 
pdf("all_all_PCA.pdf",width=8,height=9)	
grid.arrange(
	myglist[[1]]+annotate("text",label=paste("A"), x=20, y=3,size=5),
	myglist[[4]]+annotate("text",label=paste("D"), x=50, y=3,size=5)+scale_y_continuous(breaks=c(-4, 0, 4)),
	myglist[[2]]+annotate("text",label=paste("B"), x=20, y=3,size=5),
	myglist[[5]]+annotate("text",label=paste("E"), x=25, y=3,size=5),
	myglist[[3]]+annotate("text",label=paste("C"), x=20, y=3,size=5),
	myglist[[6]]+annotate("text",label=paste("F"), x=25, y=3,size=5),
	l1,l2, gblank, layout_matrix=lay
)	
dev.off()
	
#### Orchard specific ###
tempiom <- myfiltbiom
tempiom@otu_table@.Data <-  assay(varianceStabilizingTransformation(phylo_to_des(tempiom)))

filterFun=function(o,f){
  prune_samples(sample_data(o)[[11]]==f&sample_data(o)[[1]]!="C",o)
}	

mypca <- list(
  dessert=plotPCA(tempiom,returnData=T,trans=F,filterFun=filterFun,filter="Dessert"),
  cider=plotPCA(tempiom,returnData=T,trans=F,filterFun=filterFun,filter="Cider")
)
rm(tempiom)

myfiltbiom <- prune_samples(sample_data(myfiltbiom)[[1]]!="C",myfiltbiom)
myfiltbiom <- prune_taxa(rowSums(otu_table(myfiltbiom))>5,myfiltbiom)

myfiltbiom <- list( 
  dessert=filterFun(myfiltbiom,"Dessert"),
  cider=filterFun(myfiltbiom,"Cider")
)

# get the sum of squares for tree/aisle, location and residual
sum_squares <- list(
  dessert=t(apply(mypca$dessert$x,2,function(x) 
    t(summary(aov(x~Sample+location,data.frame(as.matrix(sample_data(myfiltbiom$dessert)))))[[1]][2]))),
  cider=t(apply(mypca$cider$x,2,function(x) 
    t(summary(aov(x~Sample+location,data.frame(as.matrix(sample_data(myfiltbiom$cider)))))[[1]][2])))
)

sum_squares <- lapply(sum_squares,function(x) {colnames(x) <- c("condition","location","residual");x})    

x<-lapply(sum_squares,function(x) t(apply(x,1,prop.table)))
perVar <- list(
  dessert=x$dessert * mypca$dessert$percentVar,
  cider=x$cider * mypca$cider$percentVar
)
x1 <- lapply(perVar,colSums)
x2 <- lapply(lapply(perVar,colSums),sum)
x1[[1]] <- rbind(t(data.frame(sum_sqr=x1[[1]])),"%"=x1[[1]]/x2[[1]]*100)
x1[[2]] <- rbind(t(data.frame(sum_sqr=x1[[2]])),"%"=x1[[2]]/x2[[2]]*100)

#Plot residual after removing spatial component for first couple of PCA vectors
# uncorrected
df <- lapply(mypca,function(x) {t(data.frame(t(x$x)*x$percentVar))})
# spatial removed
pc.res <- lapply(seq(1,2),function(x) resid(aov(mypca[[x]]$x~sample_data(myfiltbiom[[x]])$location)))
d <- lapply(seq(1,2),function(x) {t(data.frame(t(pc.res[[x]])*mypca[[x]]$percentVar))})

### 16S	
g1<-plotOrd(df[[1]][,1:2],sample_data(myfiltbiom[[1]]),
	design="Distance",shapes="Sample",continuous=T,colourScale=c("black","lightblue"),
	cluster=0.95,centers=1,xlabel="PC1",ylabel="PC2",legend=F,xlim=c(-6,6),ylim=c(-6,8),textSize=16
)
g2<-plotOrd(d[[1]][,1:2],sample_data(myfiltbiom[[1]]),
	design="Distance",shapes="Sample",continuous=T,colourScale=c("black","lightblue"),
	cluster=0.95,centers=1,xlabel="PC1",ylabel="PC2",legend=F,xlim=c(-6,6),ylim=c(-6,8),textSize=16
)

g3<-plotOrd(df[[2]][,1:2],sample_data(myfiltbiom[[2]]),
	design="Distance",shapes="Sample",continuous=T,colourScale=c("black","lightblue"),
	cluster=0.95,centers=1,xlabel="PC1",ylabel="PC2",legend=F,xlim=c(-8,15),ylim=c(-8,6),textSize=16
)
g4<-plotOrd(d[[2]][,1:2],sample_data(myfiltbiom[[2]]),
	design="Distance",shapes="Sample",continuous=T,colourScale=c("black","lightblue"),
	cluster=0.95,centers=1,xlabel="PC1",ylabel="PC2",legend=F,xlim=c(-8,15),ylim=c(-8,6),textSize=16
)

### ITS	
g1<-plotOrd(df[[1]][,1:2],sample_data(myfiltbiom[[1]]),
	design="Distance",shapes="Sample",continuous=T,colourScale=c("black","lightblue"),
	cluster=0.95,centers=1,xlabel="PC1",ylabel="PC2",legend=F,xlim=c(-8,8),ylims=c(-2.5,2.5),textSize=16
)
g2<-plotOrd(d[[1]][,1:2],sample_data(myfiltbiom[[1]]),
	design="Distance",shapes="Sample",continuous=T,colourScale=c("black","lightblue"),
	cluster=0.95,centers=1,xlabel="PC1",ylabel="PC2",legend=F,xlim=c(-8,8),ylims=c(-2.5,2.5),textSize=16
)
g3<-plotOrd(df[[2]][,1:2],sample_data(myfiltbiom[[2]]),
	design="Distance",shapes="Sample",continuous=T,colourScale=c("black","lightblue"),
	cluster=0.95,centers=1,xlabel="PC1",ylabel="PC2",legend=F,xlim=c(-5,5),ylims=c(-5,5),textSize=16
)
g4<-plotOrd(d[[2]][,1:2],sample_data(myfiltbiom[[2]]),
	design="Distance",shapes="Sample",continuous=T,colourScale=c("black","lightblue"),
	cluster=0.95,centers=1,xlabel="PC1",ylabel="PC2",legend=F,xlim=c(-5,5),ylims=c(-5,5),textSize=16
)	
#mygplots <-list(list(g1,g2),list(g3,g4))
	
mylegend <- ggplot_legend(
  plotOrd(df[[2]][,1:2],sample_data(myfiltbiom[[2]]),
	  design="Distance",shapes="Sample",continuous=T,colourScale=c("black","lightblue"),
    )	+theme(legend.direction="horizontal",legend.box="horizontal")
)
mylegend$layout$clip[mylegend$layout$name == "panel"] <- "off"  
  
lay=rbind(c(1,2),c(1,2),c(3,4),c(3,4),c(5,5))
 
pdf("ITS_PCA.pdf",width=6,height=6)
grid.arrange(
	g1+geom_text(aes(label = "A", x = 8, y = 6.5), color="black",size=3)+coord_fixed(ylims=c(-4,6.5)),
	g3+geom_text(aes(label = "C", x = 3, y = 3.5),color="black",size=3),
	g2+geom_text(aes(label = "B", x = 6.2, y = 6),color="black",size=3)+coord_fixed(ylims=c(-4,6)),
	g4+geom_text(aes(label = "D", x = 3, y = 6),color="black",size=3)+coord_fixed(ylims=c(-2,6)),
	mylegend, layout_matrix=lay
)
dev.off()

pdf("16S_PCA.pdf",width=6,height=6)	
grid.arrange(
	g1+geom_text(aes(label = "A", x = 8, y = 6.5), color="black",size=3),
	g3+geom_text(aes(label = "C", x = 14, y = 8),color="black",size=3),
	g2+geom_text(aes(label = "B", x =6, y=6), color="black",size=3),
	g4+geom_text(aes(label = "D", x =6, y = 6),color="black",size=3),
 	mylegend, layout_matrix=lay
)
	
dev.off()	
    
