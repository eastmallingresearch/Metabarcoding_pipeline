normHTS <- function(
	obj,
	controlSamples)
{

	gm_mean = function(x, na.rm=TRUE){
		exp(sum(log(x[x > 0]), na.rm=na.rm) / sum(!is.na(x)))
	}

	uniqueRuns <- unique(sub("S[0-9]*","",row.names(sample_data(obj))))

	dds.all <- sapply(uniqueRuns,function(x) phylo_to_des(
		prune_samples(sub("S[0-9]*","",row.names(sample_data(obj)))==x,obj)))
	
	c.counts <- lapply(dds.all,function(x) counts(x,normalize=T)[,colnames(x)%in%controlSamples])

	cc1 <- as.data.frame(lapply(c.counts,function(x) return(tryCatch(x[,1],error=function(e)NA))))
	cc2 <- as.data.frame(lapply(c.counts,function(x) return(tryCatch(x[,2],error=function(e)NA))))
	cc3 <- as.data.frame(lapply(c.counts,function(x) return(tryCatch(x[,3],error=function(e)NA))))
	cc4 <- as.data.frame(lapply(c.counts,function(x) return(tryCatch(x[,4],error=function(e)NA))))
	
	x2 <- t(data.frame(
		cc1=colSums(cc1)/colSums(cc1),
		cc2=colSums(cc2)/colSums(cc2),
		cc3=colSums(cc3)/colSums(cc3),
		cc4=colSums(cc4)/colSums(cc4)	
	))

	cc1[is.na(cc1)] <- apply(cc1,1,gm_mean)
	cc2[is.na(cc2)] <- apply(cc2,1,gm_mean)
	cc3[is.na(cc3)] <- apply(cc3,1,gm_mean)
	cc4[is.na(cc4)] <- apply(cc4,1,gm_mean)

	x1 <- rbind(
		estimateSizeFactorsForMatrix(cc1),
		estimateSizeFactorsForMatrix(cc2),
		estimateSizeFactorsForMatrix(cc3),
		estimateSizeFactorsForMatrix(cc4)
	) 
	
	x3 <- x1 * x2
	x4 <- apply(x3,2,function(x) median(x,na.rm=T))	
	x4[is.na(x4)] <- 1	
	
	xx <- unlist(lapply(dds.all,sizeFactors))

	xy <- xx*sapply(names(xx),function(x) x4[sub("\\..*","",x)])
	
	names(xy) <- sub(".*\\.","",names(xy))

	return(xy)

	#k <- which(is.na(x3), arr.ind=TRUE)
	#x3[k] <- apply(x3,1,gm_mean)[k[,1]]
	
	#estimateSizeFactorsForMatrix(x3)
}



