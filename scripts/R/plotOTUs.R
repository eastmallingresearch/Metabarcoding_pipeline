plotOTUs <-
function (
	countData,
	colData,
	design="time",
	colour="Treatment",
	plotsPerPage=5,
	facet=formula(~OTU),
	Ylab="log_counts",
	line="smooth",
	returnData=F
) {
	suppressPackageStartupMessages(require(viridis))

	d <- data.frame(t(countData),colData)

	d <- melt(d,
		measure.vars=colnames(d)[1:(ncol(t(countData)))],
		id.vars = colnames(d)[(ncol(t(countData))+1):ncol(d)],
		variable.name = "OTU", 
		value.name = Ylab)

	d[[Ylab]] <- d[[Ylab]]+abs(min(d[[Ylab]]))

	ymax <- max(d[[Ylab]])
	allVars <- unique(d$OTU)
	noVars <- length(allVars)
	plotSequence <- c(seq(0, noVars-1, by = plotsPerPage), noVars)
	
	if(returnData) { return(d)}
	
	sapply(seq(2,length(plotSequence)),function(i) {
		start <- plotSequence[i-1] + 1
		end <- plotSequence[i]
		tmp <- d[d$OTU %in% allVars[start:end],]
		cat(unique(tmp$OTU), "\n")
		g <- ggplot(data=tmp,aes_string(y=Ylab, x=design,colour=colour),ylim=c(0,ymax))
		g <- g + theme_classic_thin(base_size = 16) %+replace% theme(panel.border=element_rect(colour="black",size=0.25,fill=NA),legend.position="bottom")
		g <- g + scale_colour_viridis(discrete=TRUE)
		g <- g + facet_grid(facet,scales="free_x")
		g <- g + geom_point(size=2)
		if(line=="smooth"){
			g <- g + stat_smooth(method=locfit, formula=y~lp(x),se=F)
		}else {
			g <- g + geom_line()
		}
		print(g)
	})

}	