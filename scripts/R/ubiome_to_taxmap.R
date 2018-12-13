ubiome_to_taxmap <- function(
	obj#,
	#filter="rowSums(otu_table[,c(-1,-2)])>=0"
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
	otu_table <- as.tibble(cbind(taxon_id=t3,otu_table,stringsAsFactors=F))
	#otu_counts <- otu_table[eval(filter),]
	output$data <- list(
		otu_table   = otu_table#,
		#otu_counts  = otu_counts,
		#sample_data = sample_data
	)
	output	
}

calc_taxon_abund <-
function (obj, data, cols = NULL, groups = NULL, out_names = NULL)
{
  do_it <- function(count_table, cols = cols, groups = groups) {
        my_print("Summing per-taxon counts from ", length(cols),
            " columns ", ifelse(length(unique(groups)) == length(unique(cols)),
                "", paste0("in ", length(unique(groups)), " groups ")),
            "for ", length(obj$taxon_ids()), " taxa")
        obs_indexes <- obj$obs(data)
        output <- lapply(split(cols, groups), function(col_index) {
            col_subset <- count_table[, col_index]
            vapply(obs_indexes, function(i) {
            	if(length(i)>0) {
                	sum(col_subset[i, ])
                } else numeric(1)
            }, numeric(1))
        })
        output <- as.data.frame(output, stringsAsFactors = FALSE)
        return(output)
    }
    output <- do_calc_on_num_cols(obj, data, cols = cols, groups = groups,
        other_cols = NULL, out_names = out_names, func = do_it)
    output <- cbind(data.frame(taxon_id = obj$taxon_ids(), stringsAsFactors = FALSE),
        output)
    dplyr::as.tbl(output)
}		     
