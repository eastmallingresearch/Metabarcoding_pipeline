plot_ma <- function
(
	fitObj,
	xlim=c(-6,6),
	textSize=16,
	textFont="Helvetica",
	pointSize=3,
	legend=F,
	LOG=2,
	crush=T
)
{
	cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

	d <- as.data.table(fitObj)
	colnames(d) <- c("log2FoldChange","baseMean","padj")
	d$group<-"Not sig"
	d$group[d$padj<=0.05]<- "p <= 0.05" 
	d$group[abs(d$log2FoldChange)>1] <- "FC >= 2"	
	d$group[(d$padj<=0.05&(abs(d$log2FoldChange)>1))] <- "p <= 0.05 &\nFC >= 2" 
	d$shape<-16
	d$group<-as.factor(d$group)
	d$group<-factor(d$group,levels(d$group)[c(2,3,1,4)])

	if(crush){
		d$shape[d$log2FoldChange<xlim[1]] <- 25
		d$log2FoldChange[d$log2FoldChange<xlim[1]] <-xlim[1]
		d$shape[d$log2FoldChange>xlim[2]] <- 24
		d$log2FoldChange[d$log2FoldChange>xlim[2]] <-xlim[2]
	}

	g <- ggplot(data=d,aes(x=log2FoldChange,y=log(baseMean,LOG),colour=group,shape=shape))
	
	if(!legend) {
		g <- g + theme_classic_thin(textSize,textFont) %+replace% theme(legend.position="none")
	} else {
		g <- g + theme_classic_thin(textSize,textFont) %+replace% theme(legend.title=element_blank())
	}

	g <- g + scale_shape_identity() 
	g <- g + geom_point(size=pointSize)
	g <- g + scale_colour_manual(values=cbbPalette)
	g <- g + xlab(expression("Log"[2]*" fold change"))
	g <- g + ylab(expression("Log"[2]*" mean expression"))
	g <- g + xlim(xlim)
	g <- g + expand_limits(x = xlim[1], y = 5)
	g <- g + coord_flip()
	return(g)
}
