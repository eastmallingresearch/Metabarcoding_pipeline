##phyloseq function modified to allow return of data

plot_richness <- 
function (physeq, x = "samples", color = NULL, shape = NULL, size=2, textSize=14,
    title = NULL, scales = "free_y", nrow = 1, shsi = NULL, measures = NULL,
    sortby = NULL,
    limits,
    returnData=F	
    )
{
	library(phyloseq)
    erDF = estimate_richness(physeq, split = TRUE, measures = measures)
    measures = colnames(erDF)
    ses = colnames(erDF)[grep("^se\\.", colnames(erDF))]
    measures = measures[!measures %in% ses]
    if (!is.null(sample_data(physeq, errorIfNULL = FALSE))) {
        DF <- data.frame(erDF, sample_data(physeq))
    }
    else {
        DF <- data.frame(erDF)
    }
    if (!"samples" %in% colnames(DF)) {
        DF$samples <- sample_names(physeq)
    }
    if (!is.null(x)) {
        if (x %in% c("sample", "samples", "sample_names", "sample.names")) {
            x <- "samples"
        }
    }
    else {
        x <- "samples"
    }
    if(returnData) {
        return(DF) 
    }


    if(!missing(limits)) {
	DF<-DF[!rowSums(DF[,1:(ncol(DF)-ncol(sample_data(physeq))-1)]>limits[2]),]
    }

    mdf = melt(DF, measure.vars = measures)
    mdf$se <- NA_integer_
    if (length(ses) > 0) {
        selabs = ses
        names(selabs) <- substr(selabs, 4, 100)
        substr(names(selabs), 1, 1) <- toupper(substr(names(selabs),
            1, 1))
        mdf$wse <- sapply(as.character(mdf$variable), function(i,
            selabs) {
            selabs[i]
        }, selabs)
        for (i in 1:nrow(mdf)) {
            if (!is.na(mdf[i, "wse"])) {
                mdf[i, "se"] <- mdf[i, (mdf[i, "wse"])]
            }
        }
        mdf <- mdf[, -which(colnames(mdf) %in% c(selabs, "wse"))]
    }
    if (!is.null(measures)) {
        if (any(measures %in% as.character(mdf$variable))) {
            mdf <- mdf[as.character(mdf$variable) %in% measures,
                ]
        }
        else {
            warning("Argument to `measures` not supported. All alpha-diversity measures (should be) included in plot.")
        }
    }
    if (!is.null(shsi)) {
        warning("shsi no longer supported option in plot_richness. Please use `measures` instead")
    }
    if (!is.null(sortby)) {
        if (!all(sortby %in% levels(mdf$variable))) {
            warning("`sortby` argument not among `measures`. Ignored.")
        }
        if (!is.discrete(mdf[, x])) {
            warning("`sortby` argument provided, but `x` not a discrete variable. `sortby` is ignored.")
        }
        if (all(sortby %in% levels(mdf$variable)) & is.discrete(mdf[,
            x])) {
            wh.sortby = which(mdf$variable %in% sortby)
            mdf[, x] <- factor(mdf[, x], levels = names(sort(tapply(X = mdf[wh.sortby,
                "value"], INDEX = mdf[wh.sortby, x], mean, na.rm = TRUE,
                simplify = TRUE))))
        }
    } 

    richness_map = aes_string(x = x, y = "value", colour = color,
        shape = shape)
    p = ggplot(mdf, richness_map) + geom_point(na.rm = TRUE,position = position_dodge(width = 0.5),size=size)
    if (any(!is.na(mdf[, "se"]))) {
        p = p + geom_errorbar(aes(ymax = value + se, ymin = value -
            se), width = 0.5,position = position_dodge(width = 0.5))
    }
    p = p +theme_bw(base_size = textSize, base_family = "Helvetica")

    p = p + theme(axis.text.x = element_text(angle = -90, vjust = 0.5,
        hjust = 0),panel.grid.major.x = element_line(size = .5, color = "grey"))
    p = p + ylab("Alpha Diversity Measure")
    p = p + facet_wrap(~variable, nrow = nrow, scales = scales)
    if (!is.null(title)) {
        p <- p + ggtitle(title)
    }
    return(p)
}

# cut down version of above - useful for projects without phyloseq objects

plot_alpha <- 
function (countdata,colData,design="Site",colour="Samples",shape=NULL,returnData=F)
{	
	# function to convet input to integer values with ceiling for values < 1
	alpha_counts <- function(X) {
		X[X==0]<-NA
		X[X<1] <- 1
		X[is.na(X)] <- 0
		return(round(X,0))
	}

	# convert counts to integers and transpose 
	OTU <- t(alpha_counts(countData))

	# calculate diversity measures
	all_alpha <- data.table(
		t(estimateR(OTU)),
		shannon = diversity(OTU, index = "shannon"),
		simpson = diversity(OTU, index = "simpson"),
		keep.rownames="Samples"
	)	
	
	if(returnData) {return(all_alpha)}

	# set colData to a data table	
	colData <- as.data.table(as.data.frame(colData),keep.rownames="Samples")


	# get column names of all_alpha
	measures = colnames(all_alpha[,-1])

	# get standard error columns (Chao1 and ACE)
	ses = colnames(all_alpha)[grep("^se\\.", colnames(all_alpha))]

	# temp holder for se labels
	selabs = ses

	# rename se labels to the same as the corresponding measures column
	names(selabs) <- sub("se","S",selabs)

	# set measures to (measures - se) columns
	measures = measures[!measures %in% ses]

	all_alpha <- as.data.table(inner_join(all_alpha,colData))

	id.vars <- colnames(colData)

	# melt the data table by Sample
	mdf <- melt(all_alpha,meaures.vars=meaures,id.vars=id.vars,variable.factor = TRUE)

	X<-all_alpha[,c(id.vars,selabs),with=F]
	names(X) <- c(id.vars,names(selabs))
	X <- melt(X,id.vars=id.vars,value="se",variable.factor = TRUE)

	mdf <- left_join(mdf,X)

	# remove se rows
	mdf <- mdf[as.character(mdf$variable) %in% measures,]

	# refactor variable
	mdf$variable <- as.factor(mdf$variable)
	mdf$value <- as.numeric(mdf$value)

	# create ggplot object
	g <- ggplot(data=mdf,aes_string(x=design,y="value",colour=colour))

	# add a theme
	g <- g + theme_classic_thin() %+replace% theme(	
		panel.border = element_rect(colour = "black", fill=NA, size=0.5),
		axis.text.x = element_text(angle = -90, vjust = 0.5,hjust = 0)
	)

	# add points
	g <- g + geom_point(na.rm = TRUE,position = position_dodge(width = 0.5),size=1)

	# add error bars
	g <- g + geom_errorbar(aes(ymax = value + se, ymin = value -  se), width = 0.5,position = position_dodge(width = 0.5))

	# add y label
	g <- g + ylab("Alpha Diversity Measure")

	# Add heading to each graph
	g <- g + facet_wrap(~variable, nrow = 1,scales="free_y")

	g

}

