plot_ma <- function
(
	fitObj,
	xlim=c(-6,6),
	textsize=16,
	textFont="Helvetica",
	pointSize=3,
	legend=F,
	crush=T
)
{
	cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

	d <- fitObj[,c(2,3,6)]
	colnames(d) <- c("log2FoldChange","baseMean","padj")
	d$group<-1
	d[d$padj<=0.05,4]<-2 
	d[abs(d$log2FoldChange)>1,4]<-3 
	d[(d$padj<=0.05)&(abs(d$log2FoldChange)>1),4]<-4
	d$group<-as.factor(d$group)
	d$shape<-16

	if(crush){
		d[d$log2FoldChange<xlim[1],5]<-25
		d[d$log2FoldChange<xlim[1],1]<-xlim[1]
		d[d$log2FoldChange>xlim[2],5]<-24
		d[d$log2FoldChange>xlim[2],1]<-xlim[2]
	}

	g <- ggplot(data=d,aes(x=log2FoldChange,y=baseMean,colour=group,shape=shape))
	
	if(!legend) {
		g <- g + theme_classic_thin(textSize,textFont) %+replace% theme(legend.position="none")
	} else {
		g <- g + theme_classic_thin(textSize,textFont)
	}

	g <- g + scale_shape_identity() 
	g <- g + geom_point(size=pointSize)
	g <- g + scale_colour_manual(values=cbbPalette)
	g <- g + xlab(expression("Log"[2]*" Fold Change"))
	g <- g + ylab(expression("Log"[2]*" Mean Expression"))
	g <- g + xlim(xlim)
	g <- g + expand_limits(x = xlim[1], y = 5)
	g <- g + coord_flip()
	return(g)
}
