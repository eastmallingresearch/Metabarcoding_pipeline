nearest.vec <- function(
	x, 
	vec,
	breaks=1
) {
    smallCandidate <- findInterval(x, vec, all.inside=TRUE)
    largeCandidate <- smallCandidate + 1
    nudge <- 2 * x > vec[smallCandidate] + vec[largeCandidate]
    empty <- abs(vec[smallCanditate]-break)
    return(smallCandidate + nudge)
}