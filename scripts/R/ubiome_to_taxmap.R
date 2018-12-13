ubiome_to_taxmap <- function(
	obj
){
	suppressPackageStartupMessages(require(metacoder))
	# tibblise data
	otu_table <- as.tibble(obj[[1]],rownames="otu_id")
	sample_data  <- as.tibble(obj[[2]],rownames="sample_id")
	tax_data <- as.tibble(obj[[3]],rownames="otu_id")
	parsed_tax <- apply(tax_data,1,function(x) {
		x <- sub(".*\\(.*",NA,x[2:8])
		x <- x[!is.na(x)]
		x <- sub("_SH.*","",x)
		x <- gsub("_"," ",x)
		xx <- hierarchy(x)
		lapply(seq_along(xx$taxa),function(i) {
			xx$taxa[[i]]$rank$name <<- names(x)[i]
		})
		xx
	})
	output <- taxmap(.list = parsed_tax, named_by_rank = T)
	# set the taxon_id to the rank (rank is the lowest defined rank with a given confidence)
	t1 <- output$taxon_names()
	t2 <- sub("\\(.*","",tax_data$rank)
	t2 <- sub("_SH.*","",t2)
	t2 <- gsub("_"," ",t2)	
	t3 <- sapply(t2,function(x) names(t1[t1==x])[1])
	output$data <- list(
		otu_table = as.tibble(cbind(taxon_id=t3,otu_table,stringsAsFactors=F)),
		otu_counts = otu_table,
		sample_data = sample_data
	)
	output	
}
