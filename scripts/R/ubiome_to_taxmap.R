ubiome_to_taxmap <- function(
	obj
){
	suppressPackageStartupMessages(require(metacoder))
	# tibblise data
	otu_table <- as.tibble(obj[[1]],rownames="otu_id")
	#sample_data  <- as.tibble(obj[[2]],rownames="sample_id")
	tax_data  <- as.tibble(obj[[3]],rownames="otu_id")
	parsed_tax <- apply(td[1:3,2:8],1, function(x) {
		x <- x[!is.na(x)]
		xx<-hierarchy(x)
		xx <- lapply(seq_along(xx),function(i) {
			xx$taxa[[i]]$rank$name <- names(x)[i]
			xx
		})
		xx
	})
#		xx$taxa[[1]]$rank$name <- "kingdom"
#		xx$taxa[[2]]$rank$name <- "phylum"
#		xx$taxa[[3]]$rank$name <- "class"
#		xx$taxa[[4]]$rank$name <- "order"
#		xx$taxa[[5]]$rank$name <- "family"
#		xx$taxa[[6]]$rank$name <- "genus"
#		xx$taxa[[7]]$rank$name <- "species"	
	#	parsed_tax <- lapply(parsed_tax, trimws) # removes white space (should't be any)
	output <- taxmap(.list = parsed_tax, named_by_rank = T)
	# set the taxon_id to the rank (rank is the lowest defined rank with a given confidence)
	t1 <- output$taxon_names()
	t2 <- sub("\\(.*","",tax_data$rank) 
	t3 <- sapply(t2,function(x) names(t1[t1==x])[1])
	output$data$otu_table<-as.tibble(cbind(taxon_id=t3,otu_table))
	output	
}

