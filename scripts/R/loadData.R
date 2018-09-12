loadData <-
function(
	countData,
	colData,
	taxData,
	phylipData,
	tax_order=c(1,3,5,7,9,11,13,2,4,6,8,10,12,14),
	tax_conf=0.65,
	RHB=NULL
){
	if(missing(countData)|missing(colData)|missing(taxData)) {
		if(missing(countData)) missing_data("countData")
		if(missing(colData)) missing_data("colData")
		if(missing(taxData)) missing_data("taxData")
		stop("Missing values")
	}

	# load otu count table
	countData <- read.table(countData,header=T,sep="\t",row.names=1, comment.char = "")

	# load sample metadata
	colData <- read.table(colData,header=T,sep="\t",row.names=1)

	# load taxonomy data
	taxData <- read.table(taxData,header=F,sep=",",row.names=1)

	# reorder columns
	taxData<-taxData[,tax_order]

	# add best "rank" at 0.65 confidence and tidy-up the table
	taxData<-phyloTaxaTidy(taxData,tax_conf)

	# save data into a list
	ubiom <- list(
		countData=countData,
		colData=colData,
		taxData=taxData,
		RHB=RHB
	)

	# get unifrac dist
	if(!missing(phylipData)) ubiom$phylipData <- fread.phylip(phylipData)

	return(ubiom)

	missing_data<-function(x) {
		print(paste(x, "must be specified"))
	}
}