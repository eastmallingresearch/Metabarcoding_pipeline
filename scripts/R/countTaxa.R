countTaxa <- function(
	obj,
	taxon="phylum"
){
	# sums number of unique entries at given level 
	suppressPackageStartupMessages(require(dplyr))
	#suppressPackageStartupMessages(require(data.table))
	data.frame(Taxa=as.data.frame(as.matrix(obj))[[taxon]])

	dtx <-data.frame(Taxa=as.data.frame(as.matrix(obj))[[taxon]])
	dtx %>% group_by(Taxa) %>% summarise(count=length(Taxa)) #dply method
	
#	return(setDT(dtx)[, .N, keyby=Taxa]) # data table method - doesn't work in function??
	
}

countTaxa2 <- 
function(
        obj,
        taxon="phylum"
){
        # sums number of unique entries at given level
        suppressPackageStartupMessages(require(data.table))
        data.table(Taxa=as.data.frame(as.matrix(obj))[[taxon]])[, .N, keyby=Taxa]
}