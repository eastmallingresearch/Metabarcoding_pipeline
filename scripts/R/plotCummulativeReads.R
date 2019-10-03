plotCummulativeReads <- 
function(countData,cutoffs=c(0.8,0.9,0.99,0.999),returnData="dtt",plot=TRUE,bysample=F)
{
	
	require(ggplot2)
	require(data.table)
	
	# calculate row sums from normalised read counts
    if(bysample){
		DT <- data.table(apply(countData,2,sort,decreasing=T))
   	} else {
		DT <- data.table(CD=rowSums(countData),OTU=rownames(countData))
		setorder(DT,-CD)
	}


#	DT <- DT[order("CD",decreasing=T)]
	print(paste("Returning data table as",returnData))
 	assign(returnData, DT, envir=globalenv())
	if(!plot)return("Not plotting")

	# calculate cumulative sum of ordered OTU counts 
    #DT <- cumsum(DT)
    if(!bysample){DT<-DT[-2]}
	DT <- cumsum(DT)
	#DT$CD <- cumsum(DT[,"CD"])

	suppressWarnings(if(cutoffs&!bysample) {
		# get cutoff poisiotns
		mylines <- data.table(RHB=sapply(cutoffs,function(i) {nrow(DT) - sum(DT$CD >= max(DT$CD,na.rm=T)*i)}),Cutoffs=cutoffs)	
#mylines <- data.table(RHB=sapply(cutoffs,function(i) {nrow(DT) - sum(DT >= max(DT,na.rm=T)*i)}),Cutoffs=cutoffs)	
	}) 

	# log the count values
	DT$CD <- log10(DT[,"CD"])

	# create an empty ggplot object from the data table
	g <- ggplot(data=DT,aes_string(x=seq(1,nrow(DT)),y="CD"))

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


plotRarefaction <- function(X,cutOff) {
  
  cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")  
  DT <- as.data.table(X)
  DT <- DT[,lapply(.SD,sort,decreasing=T)]
  DT <- cumsum(DT)
  
  
  DT <- DT[,lapply(.SD,log10)]
  DT[,Count:=1:nrow(DT)]
  
  DT <- melt(DT,id.vars="Count")
  
  g <- ggplot(DT[Count<=cutOff,],aes(x=Count,y=value,colour=variable))    
  g + 
    geom_line(size = 1.5)  + 
    scale_colour_manual(values = cbPalette) + 
    theme_classic_thin() %+replace% theme(legend.position = "none") + 
    xlab("Number of OTUs") + ylab(expression(~Log[10]~" sequence count"))
}
