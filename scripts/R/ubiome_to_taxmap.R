ubiome_to_taxmap <- function(
	obj
){
	suppressPackageStartupMessages(require(metacoder))
	# tibblise data
	otu_table <- as.tibble(obj[[1]],rownames="otu_id")
	sample_data  <- as.tibble(obj[[2]],rownames="sample_id")
	tax_data  <- as.tibble(obj[[3]],rownames="otu_id")
	#phy_data  <- if(!is.null(obj[[4]])){obj[[4]]} else{NULL}

	parsed_tax <- lapply(unlist(tax_data[,2:8]), trimws)
	output <- taxmap(.list = parsed_tax, named_by_rank = F)
	output$data <- list(otu_table=otu_table,tax_data=tax_data,sample_data=sample_data)
	output
}

