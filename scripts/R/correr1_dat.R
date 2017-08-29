correr1_dat <- function(x) {
	y <- x
	mycorr <- NULL
	count<-1
	while (length(x) >2) {
		mycorr[[count]] <- cbind(y[1:(length(x)-1)],x[2:length(x)])
		count <- count +1
		x<-x[-1]
	}
	return(mycorr)
}