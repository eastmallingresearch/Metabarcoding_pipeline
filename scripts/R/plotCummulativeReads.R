plotCummulativeReads <- 
function(countData,cutoffs=c(0.8,0.9,0.99,0.999))
{
	
	# calculate row sums from normalised read counts
	df <- data.table(CD=rowSums(countData),keep.rownames=F)

	# calculate cumulative sum of ordered OTU counts 
	df <- cumsum(df[order("CD",decreasing=T)])
 
	suppressWarnings(if(cutoffs) {
		# get cutoff poisiotns
		mylines <- data.table(RHB=sapply(cutoffs,function(i) {nrow(df) - sum(df >= max(df,na.rm=T)*i)}),Cutoffs=cutoffs)	
	}) 

	# log the count values
	df <- log10(df)

	#if(returnData){return(df)}

	# create an empty ggplot object from the data table
	g <- ggplot(data=df,aes_string(x=seq(1,nrow(df)),y="CD"))

	# remove plot background and etc.
	g <- g + theme_classic_thin()

	# a colour pallete that can be descriminated by the colur blind 
	cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

	# plot cumulative reads
	g <- g + geom_line(size=1.5) + scale_colour_manual(values=cbbPalette) 

	if(exists("mylines")) {
		# add cutoff lines
		g <- g + geom_vline(data=mylines,aes(xintercept=RHB),colour=cbbPalette[2:(length(cutoffs)+1)])
	
		# label cutoff lines
		g <- g + geom_text(data = mylines, aes(label=Cutoffs,x=RHB, y = -Inf),angle = 90, inherit.aes = F, hjust = -.5, vjust = -.5)
	}

	# add axis lables
	g <- g + ylab(expression("Log"[10]*" aligned sequenecs"))+xlab("OTU count")

	# print the plot
	g
}
