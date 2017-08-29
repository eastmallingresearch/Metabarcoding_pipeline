calcCorrelog<- function(
	pca,
	obj,
	pc,
	na.add,
	condition=c("Y","N"),
	condColumn=1,
	returnInp,
	returnCD,
	useMeans=F
){
	cond<-condition[1]

	pc.x <- scores(pca)[rownames(scores(pca))%in%rownames(sample_data(obj)[sample_data(obj)[[condColumn]]==cond]),]
	col.x <- sample_data(obj)[sample_data(obj)[[condColumn]]==cond,]

	pc.dt<- data.table(merge(pc.x,col.x,by="row.names"))
	if(useMeans) {
		pc.reshape <- dcast(pc.dt,Distance~.,fun.aggregate=function(x) mean(x,na.rm=T),value.var=c(names(pc.dt)[grep("PC",names(pc.dt))]))
	} else {
		pc.reshape <-pc.dt
	}
	names(pc.reshape)[grep("PC",names(pc.reshape))] <- sub("_.*","",names(pc.reshape)[grep("PC",names(pc.reshape))])
	pc.reshape<-pc.reshape[order(pc.reshape$Distance),]
	dvec <- pc.reshape$Distance

	if (!missing(na.add)) {
		inp <- pc.reshape[[pc]]
		dvec<- unlist(sapply(1:length(inp),function(i) if(i%in%na.add){return(c(mean(c(dvec[i],dvec[(i-1)])),dvec[i]))}else{return(dvec[i])}))
		inp1 <- sapply(1:length(inp),function(i) if(i%in%na.add){return(c(NA,inp[i]))}else{return(inp[i])})
	}else {
		inp1<-pc.reshape[[pc]]
	}
	
	cond<-condition[2]

	pc.x <- scores(pca)[rownames(scores(pca))%in%rownames(sample_data(obj)[sample_data(obj)[[condColumn]]==cond]),]
	col.x <-  sample_data(obj)[sample_data(obj)[[condColumn]]==cond,]
	pc.dt<- data.table(merge(pc.x,col.x,by="row.names"))
	if(useMeans) {
		pc.reshape <- dcast(pc.dt,Distance~.,fun.aggregate=function(x) mean(x,na.rm=T),value.var=c(names(pc.dt)[grep("PC",names(pc.dt))]))
	} else {
		pc.reshape <-pc.dt
	}	

	names(pc.reshape)[grep("PC",names(pc.reshape))] <- sub("_.*","",names(pc.reshape)[grep("PC",names(pc.reshape))])	

	if (!missing(na.add)) {
		inp <- pc.reshape[[pc]]
		inp2 <- sapply(1:length(inp),function(i) if(i%in%na.add){return(c(NA,inp[i]))}else{return(inp[i])})
	}else {
		inp2<-pc.reshape[[pc]]
	}
	
	if(returnInp) {
		return(cbind(unlist(inp1),unlist(inp2),dvec))
	}

	ct1 <- correr1(unlist(inp1),returnCD)
	ca1 <- correr1(unlist(inp2),returnCD)

	d<-as.data.frame(cbind(ct1,ca1,dvec[1:(length(dvec)-2)]))
	d
}