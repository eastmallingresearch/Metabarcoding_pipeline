ubiome_to_taxmap <- function(
	obj
){
	suppressPackageStartupMessages(require(metacoder))
	# tibblise data
	otu_table <- as.tibble(obj[[1]],rownames="otu_id")
	#sample_data  <- as.tibble(obj[[2]],rownames="sample_id")
	tax_data  <- as.tibble(obj[[3]],rownames="otu_id")
	#phy_data  <- if(!is.null(obj[[4]])){obj[[4]]} else{NULL}
	parsed_tax <- lapply(seq_len(nrow(tax_data)), function(i) {
		class_source <- unlist(lapply(tax_data[i, 2:8], as.character))
		unlist(class_source)
	})
	
	#	parsed_tax <- lapply(parsed_tax, trimws) # removes white space (should't be any)
	output <- taxmap(.list = parsed_tax, named_by_rank = F)
	# set the taxon_id to the rank (rank is the lowest defined rank with a given confidence)
	t1 <- output$taxon_names()
	t2 <- sub("\\(.*","",tax_data$rank) 
	t3 <- sapply(t2,function(x) names(t1[t1==x])[1])
	output$data$otu_table<-as.tibble(cbind(taxon_id=t3,otu_table))
	output	
}

