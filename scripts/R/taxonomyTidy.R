taxonomyTidy <- function(x) {
	if (x[2]=="unknown") {x[2] <- paste(x[1],"(k)",sep="")}
	if (x[3]=="unknown") {if(any(grep('\\(',x[2]))) {x[3]<-x[2]}else{x[3]<-paste(x[2],"(p)",sep="")}}
	if (x[4]=="unknown") {if(any(grep('\\(',x[3]))) {x[4]<-x[3]}else{x[4]<-paste(x[3],"(c)",sep="")}}
	if (x[5]=="unknown") {if(any(grep('\\(',x[4]))) {x[5]<-x[4]}else{x[5]<-paste(x[4],"(o)",sep="")}}
	if (x[6]=="unknown") {if(any(grep('\\(',x[5]))) {x[6]<-x[5]}else{x[6]<-paste(x[5],"(f)",sep="")}}
	if (x[7]=="unknown") {if(any(grep('\\(',x[6]))) {x[7]<-x[6]}else{x[7]<-paste(x[6],"(g)",sep="")}}
	return(x)
}