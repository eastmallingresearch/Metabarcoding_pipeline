plotCorrelog <- function(
	pca,
	obj,
	pc="PC1",
	cutoff=15,
	xlim=NULL,
	ylim=NULL,
	na.add,
	returnData=F,
	returnInp=F,
	returnCD=F,
	data,
	legend=T,
	lpos="right",
	cols=c("#edf8b1","#2c7fb8"),
	lineWidth=1.5,
	useMeans=F,
	textSize=12,
	design=c("Tree","Grass") 
) {

	if(returnInp|returnCD){
		returnData=T
	}
	
	if(!missing(data)){
		d<-data
	} else {
		d <- calcCorrelog(pca,obj,pc,na.add,design,1,returnInp,returnCD,useMeans)
	}
	
	if(returnData) {
		return(d)
	}
	d<-d[order(d$V3),]
	d<- d[1:(length(d$V3[d$V3<=cutoff])),]
	names(d) <- c("Tree","Grass","Distance")
	d2 <- melt(d,id="Distance")
	colnames(d2) <- c("Distance","Sample","Correlation")

	g <- ggplot(d2)
	g <- g + coord_fixed(ratio = 10, xlim = xlim, ylim = ylim, expand = TRUE)
	
	# g <- g + theme_classic()

#	g <- g + theme_bw()
#	g <- g + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#		legend.position=lpos,legend.text=element_text(size=10),legend.title=element_text(size=10), legend.background = element_rect(fill=alpha('white', 0)))
# 	g <- g + theme(axis.line.x = element_line(size=0.5,colour = "black"),axis.line.y = element_line(size=0.5,colour = "black"),axis.text = element_text(colour = "black"))

	if(!legend) {
		g <- g + theme_classic(base_size=textSize) %+replace% theme(legend.position="none")
	} else {
		g <- g + theme_classic(base_size=textSize) + 
			theme(
			legend.position=lpos,
			legend.text=element_text(size=10),
			legend.title=element_text(size=10),
			legend.background = element_rect(fill=alpha('white', 0))
		)
	}


	g <- g + geom_line(na.rm=T,aes(x=Distance, y=Correlation, colour=Sample),size=lineWidth)
	g <- g + scale_colour_manual(values=cols)
	g
}