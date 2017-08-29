correr2 <- function(
	X,
	breaks
){
	y <- seq(0,max(X[,2])+1,breaks)
	vec <- X[,2] 
	test <- nearest.vec(y,vec)	
	test[9] <- NA
	x <- X[X[,2]==test,1]
	mycorr <- numeric(0)
	count<-1
	while (length(x) >2) {
		mycorr[count] <- cor(y[1:(length(x)-1)],x[2:length(x)],use="pairwise.complete.obs")
		count <- count +1
		x<-x[-1]
	}
	return(mycorr)
}