library(data.table)
load_all("../..//metabarcoding_pipeline/scripts/myfunctions")
library(scales)    


myfiltbiom[[1]]@sam_data$gap <- 0
myfiltbiom[[2]]@sam_data$gap <- 0
levels(sample_data(myfiltbiom[[1]])[[1]]) <- c("Grass","Tree")
levels(sample_data(myfiltbiom[[2]])[[1]]) <- c("Grass","Tree")

temp1 <- prune_samples(sample_data(myfiltbiom[[1]])$Orchard=="Dessert",myfiltbiom[[1]])
temp2 <- prune_samples(sample_data(myfiltbiom[[1]])$Orchard=="Cider",myfiltbiom[[1]])
temp3 <- prune_samples(sample_data(myfiltbiom[[2]])$Orchard=="Dessert",myfiltbiom[[2]])
temp4 <- prune_samples(sample_data(myfiltbiom[[2]])$Orchard=="Cider",myfiltbiom[[2]])

myfiltbiom <- list(
	Bacteria_Dessert=temp1,
	Bacteria_Cider=temp2,
	Fungi_Dessert=temp3,
	Bacteria_Cider=temp4
)	
			
mypca <- lapply(myfiltbiom,function(obj) plotPCA(obj,design="1",returnData=T))
		
# Pearson Correlogram
cutoff <- 17
### For cider samples - c_fix corrects a slight issue
c_fix <- function(p,o,pn){
  t1 <- plotCorrelog(p,prune_samples(sample_data(o)$block!=3,o),pn,na.add=c(9),returnCD=T)
  t2 <- plotCorrelog(p,prune_samples(sample_data(o)$block==3,o),pn,returnCD=T)
  t3 <- sapply(1:nrow(t1),function(i) 
    if(i<=nrow(t2)){cbind(rbind(t1[[i,1]],t2[[i,1]]),rbind(t1[[i,2]],t2[[i,2]]))}else{cbind(t1[[i,1]],t1[[i,2]])}
  )
  d <- as.data.frame(t(sapply(1:length(t3),function(i)
    diag(cor(t3[[i]],use="pairwise.complete.obs")[c(1,3),c(2,4)])))
  )
  d$V3 <- as.numeric(t1[[3]])
  return(d)  
}
		       
g1<-plotCorrelog(mypca[[1]],myfiltbiom[[1]],"PC1",cutoff,xlim=NULL,ylim=c(-1,1),na.add=c(9,17),cols=c("black","lightblue"),legend=F)
g2<-plotCorrelog(data=c_fix(mypca[[2]],myfiltbiom[[2]],"PC1"), cutoff,pc="PC1",ylim=c(-1,1),cols=c("black","lightblue"),legend=F)
g3<-plotCorrelog(mypca[[1]],myfiltbiom[[1]],"PC2",cutoff,xlim=NULL,ylim=c(-1,1),na.add=c(9,17),cols=c("black","lightblue"),legend=F)
g4<-plotCorrelog(data=c_fix(mypca[[2]],myfiltbiom[[2]],"PC2"),cutoff,pc="PC2",ylim=c(-1,1),cols=c("black","lightblue"),legend=F)

g11<-plotCorrelog(mypca[[3]],myfiltbiom[[3]],"PC1",cutoff,xlim=NULL,ylim=c(-1,1),na.add=c(9,17),cols=c("black","lightblue"),legend=F)
g12<-plotCorrelog(data=c_fix(mypca[[4]],myfiltbiom[[4]],"PC1"), cutoff,pc="PC1",ylim=c(-1,1),cols=c("black","lightblue"),legend=T,lpos=c(0.275,0.25))
g13<-plotCorrelog(mypca[[3]],myfiltbiom[[3]],"PC2",cutoff,xlim=NULL,ylim=c(-1,1),na.add=c(9,17),cols=c("black","lightblue"),legend=F)
g14<-plotCorrelog(data=c_fix(mypca[[4]],myfiltbiom[[4]],"PC2"),cutoff,pc="PC2",ylim=c(-1,1),cols=c("black","lightblue"),legend=F)			      
		      
lay=cbind(c(1,3),c(2,4),c(5,7),c(6,8))			      
pdf("all_correlog_bb2.pdf",width=8,height=5)
grid.arrange(
	g11+geom_text(aes(label = "A", x = 15, y = 1), color="black",size=3),
	g13+geom_text(aes(label = "C",  x = 15, y = 1),color="black",size=3),
	g12+geom_text(aes(label = "B",  x = 15, y = 1),color="black",size=3),
	g14+geom_text(aes(label = "D",  x = 15, y = 1),color="black",size=3),
	g1+geom_text(aes(label = "E", x = 15, y = 1), color="black",size=3),
	g3+geom_text(aes(label = "G",  x = 15, y = 1),color="black",size=3),
	g2+geom_text(aes(label = "F",  x = 15, y = 1),color="black",size=3),
	g4+geom_text(aes(label = "H",  x = 15, y = 1),color="black",size=3),
	layout_matrix=lay
)
dev.off()			      

			      
pdf("correlog.pdf",width=6,height=6)
grid.arrange(
	g1+geom_text(aes(label = "A", x = 15, y = 1), color="black",size=3),
	g3+geom_text(aes(label = "C",  x = 15, y = 1),color="black",size=3),
	g2+geom_text(aes(label = "B",  x = 15, y = 1),color="black",size=3),
	g4+geom_text(aes(label = "D",  x = 15, y = 1),color="black",size=3)
)
dev.off()
			      
#gt <- gtable(widths = unit(c(1, 2,3,4), "null"), heights = unit(c(1, 2, 3,4), "null"))
#gt <- gtable_add_grob(gt, ga, t = 1, b = 4, l = 4, r = 1)
#gt <- gtable_add_grob(gt, mylegend, t = 4, r = 4,l=1)
#grid.draw(gt)
