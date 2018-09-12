as.number <-
function(f,convert=as.numeric) {
	convert(levels(f))[f]
}