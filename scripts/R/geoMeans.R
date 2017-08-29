# function sometimes useful for replacing calcfactors
geoMeans <- function(d,dummy) {
	gm_mean = function(x, na.rm=TRUE){
		exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
	}	
	gm = apply(counts(d), 1, gm_mean)
	sizeFactors(estimateSizeFactors(d, geoMeans =gm))
}

geoSet <- function(d,dummy) {
	tryCatch ({
		sizeFactors(estimateSizeFactors(d))
	}, error = function(e) {
		cat("Sizefactors: Using geoMeans","\n")
		return(geoMeans(d))
	})
	
}