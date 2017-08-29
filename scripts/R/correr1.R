correr1 <- function(
	x,
	returnData=F
) {
	y <- x
	count<-1
	mycorr <- NULL
	while (length(x) >2) {
		if(returnData) {
			mycorr[[count]] <- cbind(y[1:(length(x)-1)],x[2:length(x)])
		} else {
			mycorr[count] <- cor(y[1:(length(x)-1)],x[2:length(x)],use="pairwise.complete.obs")
		}
		count <- count +1
		x<-x[-1]
	}
	return(mycorr)
}
