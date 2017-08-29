prcomp2 <-
function(obj) {
	X <- prcomp(t(obj))
	X$percentVar<-X$sdev^2/sum(X$sdev^2)
	X
}