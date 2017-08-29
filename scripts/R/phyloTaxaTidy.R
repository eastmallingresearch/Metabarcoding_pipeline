phyloTaxaTidy <- function(obj,...) {
	if(ncol(obj)==7){
		colnames(obj) <- c(
			"kingdom","phylum","class", "order", "family","genus", "species"
		)
	} else {
		colnames(obj) <- c(
			"kingdom","phylum","class", "order", "family","genus", "species",
			"k_conf", "p_conf","c_conf","o_conf","f_conf","g_conf","s_conf"
		)
		
		obj[,8:14] <- as.numeric(unlist(obj[,8:14]))
		obj[obj[,8:14]==100,8:14] <- 0
		#obj[obj[,8:14]==-1,8:14] <- 1
		obj <- taxaConf(obj,...)
		obj <- obj[,c(1:7,15,8:14)]
	}
	#obj <- sub("_+"," ",obj)
	
	obj[,1:7] <- t(apply(obj[,1:7],1,taxonomyTidy))
	return(obj)
}
