plotHeatmap <- 
function (dist,
	coldata,
	axis.labels=F,
	textSize=11,
	textFont="Helvetica"

) {
	require(ggplot2)
	#require(reshape)
	
	if (!missing(coldata)) {
		dist<-merge(dist,coldata,by="row.names",sort=F)
	}

	h1 <- melt(dist)
	colnames(h1) <-c("X","Y","Unidist")
	g <- ggplot(h1,aes(x=X,y=Y,fill=Unidist))
	if(axis.labels) {
		g <- g + theme_classic(textSize,textFont) %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1,),axis.title=element_blank())
	} else 	{
		g <- g + theme_classic(textSize,textFont) %+replace% theme(axis.text=element_blank(),axis.ticks=element_blank())
	}	

	g <- g + geom_raster()

	g <- g + scale_x_discrete("Var1")
	g <- g + scale_y_discrete("Var2")

	#g <- g + geom_tile(aes(fill=Unidist),colour="black")
	g <- g + scale_fill_gradient(low="#000033", high = "#66CCFF", na.value = "black")
	return(g)

}



