plotOrd <- function (
	obj,
	colData,
	design,
	shapes,
	labels=F,
	sublabels=F,
	cluster=NULL,
	continuous=F,
	colourScale=c(low="red", high="yellow"), #greyscale low="#000000", high="#DCDCDC"
	cbPalette=F,
	pointSize=2,
	textSize=11,
	textFont="Helvetica",
	xlims=NULL,
	ylims=NULL,
	legend=T,
	title=NULL,
	xlabel,
	ylabel,
	axes=c(1,2),
	dimx=1,
	dimy=2,
	alpha=NULL,
	exclude=F, # sometimes it can be useful to include a vector of points to exclude 
	...
) {

	obj     <- obj[,!exclude]
	colData <- colData[,!exclude]

	
	if(!is.null(title)){if(title=="debugging"){invisible(mapply(assign, ls(),mget(ls()), MoreArgs=list(envir = globalenv())));return(title)}}

	suppressPackageStartupMessages(require(ggplot2))
	suppressPackageStartupMessages(require(viridis))

	d <- as.data.frame(obj[,axes])
	#as.data.frame(obj[, axes[1]], obj[, axes[2]])

	if(!is.null(cluster)) {
	#	km <- kmeans(obj,...)
	#	colData$Cluster<-as.factor(km$cluster)
	}

	x = colnames(d)[1]
	y = colnames(d)[2]

	aes_map <-aes_string(x=x,y=y)

	if (!missing(design)) {
		colour <- if (length(design) > 1) {
			factor(apply(as.data.frame(colData[, design,drop = FALSE]), 1, paste, collapse = " : "))
		}
		else {
			as.factor(colData[[design]])
		}
		d <- cbind(d,colour=colour)
		if(continuous) {
			aes_map <- modifyList(aes_map,aes_string(colour=as.number(colour)))
		} else {
			aes_map <- modifyList(aes_map,aes_string(colour=colour))
		}
	}

	if (!missing(shapes)) {
		shape <- if (length(shapes) > 1) {
			factor(apply(as.data.frame(colData[, shapes,drop = FALSE]), 1, paste, collapse = " : "))
		} else {
			as.factor(colData[[shapes]])
		}
		d <- cbind(d,shapes = shape)
		aes_map <- modifyList(aes_map,aes_string(shape="shapes"))
	}



	if(labels) {
		label <- rownames(d)
		label[sublabels] <- NA
		aes_map <- modifyList(aes_map,aes_string(label="label"))
	}
 #return(aes_map)
	cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

	### ggplot ###
	g <- ggplot(data=d,aes_map) + geom_point(size=pointSize,na.rm = TRUE)
	g <- g + coord_fixed(ratio = 1, xlim = xlims, ylim = ylims, expand = TRUE)

	g <- g + theme_classic_thin(textSize,textFont)

	if(!legend) {
		g <- g + theme_classic_thin(textSize,textFont) %+replace% theme(legend.position="none")
	}

	if(!missing(design)) {
#		g <- g + aes_string(colour=colour)
		if(continuous) {
#			g <- g + aes_string(colour=as.number(colour))
			g <- g + scale_colour_gradient(low=colourScale[1], high=colourScale[2],name=design)
		} else {
			if(cbPalette) {
				g<-g+scale_colour_manual(values=cbbPalette)	+ guides(colour=guide_legend(title=design))
			} else {
				g<-g+scale_colour_viridis(discrete=TRUE)+ guides(colour=guide_legend(title=design))
			}
		}
	}

	if (!missing(shapes)) {
		#g <- g + aes(shape=shapes)
		g <- g + scale_shape_discrete(name=shapes)
	}

#	g <- g + geom_point(size=pointSize)

	if(!is.null(alpha)) g <-g+ geom_point(size=(pointSize+(pointSize*1.5)),alpha=alpha)

	if(labels) {
		#g <- g + aes(label=row.names(obj))
		g <- g + geom_text(size=(pointSize+2), vjust=-1.5, hjust=0.5)
	}

	if(!is.null(cluster)) {
			g<-g+stat_ellipse(geom="polygon", level=cluster, alpha=0.2)
	}

	if (!missing(xlabel)) {g <- g + xlab(xlabel)}
	if (!missing(ylabel)) {g <- g + ylab(ylabel)}
	return(g)
}
