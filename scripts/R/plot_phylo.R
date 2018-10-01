##phyloseq alpha diversity function modified to allow return of data

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

#	DF[,3:4] <- round(DF[,3:4],2)

    if(!missing(limits)) {
	  DF <- DF[DF[,1]>=limits[1],]
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
    p = ggplot(mdf, richness_map) + geom_point(na.rm = TRUE,position = position_dodge(width = 0.75),size=size)
    if (any(!is.na(mdf[, "se"]))) {
        p = p + geom_errorbar(aes(ymax = value + se, ymin = value -
            se), width = 0.2,position = position_dodge(width = 0.75))
    }
    p = p +theme_bw(base_size = textSize, base_family = "Helvetica")
#axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0),
    p = p + theme(
      panel.grid.major.x = element_line(size = .5, color = "lightgrey"),
      panel.grid.major.y =element_blank())
    p = p + ylab("Alpha Diversity Measure")
    p = p + facet_wrap(~variable, nrow = nrow, scales = scales)
    if (!is.null(title)) {
        p <- p + ggtitle(title)
    }
    return(p)
}


# cut down version of above - useful for projects without phyloseq objects
plot_alpha <-
function (	
  countData,
  colData,
  design=NULL,
  colour=NULL,
  shape=NULL,
  returnData=F,
  limits,
  measures=NULL,
  discrete=T,
  legend="right",
  pointSize=2,
  cbPalette=F,
  ...
){
	suppressPackageStartupMessages(require(viridis))

	simpleCap <- function(x) {
        	s <- strsplit(x, " ")[[1]]
        	paste(toupper(substring(s, 1,1)), substring(s, 2),
        	sep="", collapse=" ")
    	}

	retain<-measures
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

	# set colData to a data table (no ides why this is to data frame first...)
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
	
	# capitalise indices
	mdf$variable <- sub("S\\.","",mdf$variable)
	mdf$variable <- sapply(mdf$variable,simpleCap)

	# refactor variable
	mdf$variable <- sub("Obs","Observed",mdf$variable)
	mdf$variable <- as.factor(mdf$variable)
	mdf$value <- as.numeric(mdf$value)

	# retain indices were interested in (o.k. this would be better done before calculating them)
	if (!is.null(retain)) {
	    retain<-sapply(tolower(retain),simpleCap)
    	mdf <- mdf[as.character(mdf$variable) %in% retain,]
    	mdf <- droplevels(mdf)
    	mdf$variable <- factor(mdf$variable, levels = retain)
    }
	
	arguments <- list(...)
	if("debugging"%in%names(arguments))if(arguments$debugging) return(mdf)

	if(!missing(limits)) {
		mdf <- mdf[!(mdf$variable==limits[3]&(mdf$value<as.numeric(limits[1])|mdf$value>as.numeric(limits[2]))),]
	}

	cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

	aes_map <-aes_string(y="value",variable="variable")
	if(length(design)) aes_map <- modifyList(aes_map,aes_string(x=design))
	if(length(colour)) aes_map <- modifyList(aes_map,aes_string(colour=colour))
	if(length(shape)) aes_map <- modifyList(aes_map,aes_string(shape=shape))

	# create ggplot object
	g <- ggplot(data=mdf,aes_map)

	# add a theme
	g <- g + theme_classic_thin() %+replace% theme(
		panel.border = element_rect(colour = "black", fill=NA, size=0.5),
		axis.text.x = element_text(angle = -90, vjust = 0.5,hjust = 0),
		legend.position=legend
	)

	# add points
	g <- g + geom_point(na.rm = TRUE,position = position_dodge(width = 0.5),size=pointSize)

	# add error bars
	g <- g + geom_errorbar(aes(ymax = value + se, ymin = value -  se), width = 0.5,position = position_dodge(width = 0.5))

	# add y label
	g <- g + ylab("Alpha Diversity Measure")

	# change colours to viridis
	if(cbPalette) {
		g<-g+scale_colour_manual(values=cbbPalette) #+ guides(colour=guide_legend(title=design))
	} else {
		g<-g+scale_colour_viridis(discrete=TRUE) #+ guides(colour=guide_legend(title=design))
	}

	# Add heading to each graph
	g <- g + facet_wrap(~variable, nrow = 1,scales="free_y")

	g

}

# modified phyloseq ordination plot (allows rescaling by percentage variation in the axes)
plot_ordination <-
function (
  physeq,
  ordination,
  type = "samples",
  axes = 1:2,
  color = NULL,
  shape = NULL,
  label = NULL,
  title = NULL,
  justDF = FALSE,
  continuous=F,
  colourScale=c(low="red", high="yellow"),
  cbPalette=F,
  xlims=NULL,
  ylims=NULL,
  alpha=NULL,
  rescale=T,
  ...
)
{
	suppressPackageStartupMessages(require(ggplot2))
	suppressPackageStartupMessages(require(viridis))

    if (length(type) > 1) {
        warning("`type` can only be a single option,\\n            but more than one provided. Using only the first.")
        type <- type[[1]]
    }
    if (length(color) > 1) {
        warning("The `color` variable argument should have length equal to 1.",
            "Taking first value.")
        color = color[[1]][1]
    }
    if (length(shape) > 1) {
        warning("The `shape` variable argument should have length equal to 1.",
            "Taking first value.")
        shape = shape[[1]][1]
    }
    if (length(label) > 1) {
        warning("The `label` variable argument should have length equal to 1.",
            "Taking first value.")
        label = label[[1]][1]
    }
    official_types = c("sites", "species", "biplot", "split",
        "scree")
    if (!inherits(physeq, "phyloseq")) {
        if (inherits(physeq, "character")) {
            if (physeq == "list") {
                return(official_types)
            }
        }
        warning("Full functionality requires `physeq` be phyloseq-class ",
            "with multiple components.")
    }
    type = gsub("^.*site[s]*.*$", "sites", type, ignore.case = TRUE)
    type = gsub("^.*sample[s]*.*$", "sites", type, ignore.case = TRUE)
    type = gsub("^.*species.*$", "species", type, ignore.case = TRUE)
    type = gsub("^.*taxa.*$", "species", type, ignore.case = TRUE)
    type = gsub("^.*OTU[s]*.*$", "species", type, ignore.case = TRUE)
    type = gsub("^.*biplot[s]*.*$", "biplot", type, ignore.case = TRUE)
    type = gsub("^.*split[s]*.*$", "split", type, ignore.case = TRUE)
    type = gsub("^.*scree[s]*.*$", "scree", type, ignore.case = TRUE)
    if (!type %in% official_types) {
        warning("type argument not supported. `type` set to 'samples'.\\n",
            "See `plot_ordination('list')`")
        type <- "sites"
    }
    if (type %in% c("scree")) {
        return(plot_scree(ordination, title = title))
    }
    is_empty = function(x) {
        length(x) < 2 | suppressWarnings(all(is.na(x)))
    }
    specDF = siteDF = NULL
    trash1 = try({
        siteDF <- scores(ordination, choices = axes, display = "sites",
            physeq = physeq)
    }, silent = TRUE)
    trash2 = try({
        specDF <- scores(ordination, choices = axes, display = "species",
            physeq = physeq)
    }, silent = TRUE)
    siteSampIntx = length(intersect(rownames(siteDF), sample_names(physeq)))
    siteTaxaIntx = length(intersect(rownames(siteDF), taxa_names(physeq)))
    specSampIntx = length(intersect(rownames(specDF), sample_names(physeq)))
    specTaxaIntx = length(intersect(rownames(specDF), taxa_names(physeq)))
    if (siteSampIntx < specSampIntx & specTaxaIntx < siteTaxaIntx) {
        co = specDF
        specDF <- siteDF
        siteDF <- co
        rm(co)
    }
    else {
        if (siteSampIntx < specSampIntx) {
            siteDF <- specDF
            specDF <- NULL
        }
        if (specTaxaIntx < siteTaxaIntx) {
            specDF <- siteDF
            siteDF <- NULL
        }
    }
    if (is_empty(siteDF) & is_empty(specDF)) {
        warning("Could not obtain coordinates from the provided `ordination`. \\n",
            "Please check your ordination method, and whether it is supported by `scores` or listed by phyloseq-package.")
        return(NULL)
    }
    if (is_empty(specDF) & type != "sites") {
        message("Species coordinates not found directly in ordination object. Attempting weighted average (`vegan::wascores`)")
        specDF <- data.frame(wascores(siteDF, w = veganifyOTU(physeq)),
            stringsAsFactors = FALSE)
    }
    if (is_empty(siteDF) & type != "species") {
        message("Species coordinates not found directly in ordination object. Attempting weighted average (`vegan::wascores`)")
        siteDF <- data.frame(wascores(specDF, w = t(veganifyOTU(physeq))),
            stringsAsFactors = FALSE)
    }
    specTaxaIntx <- siteSampIntx <- NULL
    siteSampIntx <- length(intersect(rownames(siteDF), sample_names(physeq)))
    specTaxaIntx <- length(intersect(rownames(specDF), taxa_names(physeq)))
    if (siteSampIntx < 1L & !is_empty(siteDF)) {
        warning("`Ordination site/sample coordinate indices did not match `physeq` index names. Setting corresponding coordinates to NULL.")
        siteDF <- NULL
    }
    if (specTaxaIntx < 1L & !is_empty(specDF)) {
        warning("`Ordination species/OTU/taxa coordinate indices did not match `physeq` index names. Setting corresponding coordinates to NULL.")
        specDF <- NULL
    }
    if (is_empty(siteDF) & is_empty(specDF)) {
        warning("Could not obtain coordinates from the provided `ordination`. \\n",
            "Please check your ordination method, and whether it is supported by `scores` or listed by phyloseq-package.")
        return(NULL)
    }
    if (type %in% c("biplot", "split") & (is_empty(siteDF) |
        is_empty(specDF))) {
        if (is_empty(siteDF)) {
            warning("Could not access/evaluate site/sample coordinates. Switching type to 'species'")
            type <- "species"
        }
        if (is_empty(specDF)) {
            warning("Could not access/evaluate species/taxa/OTU coordinates. Switching type to 'sites'")
            type <- "sites"
        }
    }

	if (length(extract_eigenvalue(ordination)[axes]) > 0) {
        eigvec = extract_eigenvalue(ordination)
        fracvar = eigvec[axes]/sum(eigvec)
        percvar = round(100 * fracvar, 1)
    } else {
		percvar = 1
	}

    if (type != "species") {
		if(rescale) siteDF <- t(t(siteDF)*percvar)
        sdf = NULL
        sdf = data.frame(access(physeq, slot = "sam_data"), stringsAsFactors = FALSE)
        if (!is_empty(sdf) & !is_empty(siteDF)) {
            siteDF <- cbind(siteDF, sdf[rownames(siteDF), ])
        }
    }
    if (type != "sites") {
        tdf = NULL
        tdf = data.frame(access(physeq, slot = "tax_table"),
            stringsAsFactors = FALSE)
        if (!is_empty(tdf) & !is_empty(specDF)) {
            specDF = cbind(specDF, tdf[rownames(specDF), ])
        }
    }
    if (!inherits(siteDF, "data.frame")) {
        siteDF <- as.data.frame(siteDF, stringsAsFactors = FALSE)
    }
    if (!inherits(specDF, "data.frame")) {
        specDF <- as.data.frame(specDF, stringsAsFactors = FALSE)
    }
    DF = NULL
    DF <- switch(EXPR = type, sites = siteDF, species = specDF,
        {
            specDF$id.type <- "Taxa"
            siteDF$id.type <- "Samples"
            colnames(specDF)[1:2] <- colnames(siteDF)[1:2]
            DF = merge(specDF, siteDF, all = TRUE)
            if (!is.null(shape)) {
                DF <- rp.joint.fill(DF, shape, "Samples")
            }
            if (!is.null(shape)) {
                DF <- rp.joint.fill(DF, shape, "Taxa")
            }
            if (!is.null(color)) {
                DF <- rp.joint.fill(DF, color, "Samples")
            }
            if (!is.null(color)) {
                DF <- rp.joint.fill(DF, color, "Taxa")
            }
            DF
        })
    if (justDF) {
        return(DF)
    }
    if (!is.null(color)) {
        if (!color %in% names(DF)) {
            warning("Color variable was not found in the available data you provided.",
                "No color mapped.")
            color <- NULL
        }
    }



    if (!is.null(shape)) {
        if (!shape %in% names(DF)) {
            warning("Shape variable was not found in the available data you provided.",
                "No shape mapped.")
            shape <- NULL
        }
    }
    if (!is.null(label)) {
        if (!label %in% names(DF)) {
            warning("Label variable was not found in the available data you provided.",
                "No label mapped.")
            label <- NULL
        }
    }
    x = colnames(DF)[1]
    y = colnames(DF)[2]
    if (ncol(DF) <= 2) {
        message("No available covariate data to map on the points for this plot `type`")
        ord_map = aes_string(x = x, y = y)
    }
    else if (type %in% c("sites", "species", "split")) {
        ord_map = aes_string(x = x, y = y, color = color, shape = shape,
            na.rm = TRUE)
    }
    else if (type == "biplot") {
        if (is.null(color)) {
            ord_map = aes_string(x = x, y = y, size = "id.type",
                color = "id.type", shape = shape, na.rm = TRUE)
        }
        else {
            ord_map = aes_string(x = x, y = y, size = "id.type",
                color = color, shape = shape, na.rm = TRUE)
        }
    }
# return(ord_map)
    p <- ggplot(DF, ord_map) + geom_point(na.rm = TRUE)
	p <- p + coord_fixed(ratio = 1, xlim = xlims, ylim = ylims, expand = TRUE)

	if(continuous) {
		p <- p + scale_colour_gradient(low=colourScale[1], high=colourScale[2])
	} else {
		if(cbPalette) {
			p<-p+scale_colour_manual(values=cbbPalette)#	+ guides(colour=guide_legend(title=design))
		} else {
			p<-p+scale_colour_viridis(discrete=TRUE) # + guides(colour=guide_legend(title=design))
		}
	}

    if (type == "split") {
        p <- p + facet_wrap(~id.type, nrow = 1)
    }
    if (type == "biplot") {
        if (is.null(color)) {
            p <- update_labels(p, list(colour = "Ordination Type"))
        }
        p <- p + scale_size_manual("type", values = c(Samples = 5,
            Taxa = 2))
    }
    if (!is.null(label)) {
        label_map <- aes_string(x = x, y = y, label = label,
            na.rm = TRUE)
        p = p + geom_text(label_map, data = rm.na.phyloseq(DF,
            label), size = 2, vjust = 1.5, na.rm = TRUE)
    }
    if (!is.null(title)) {
        p = p + ggtitle(title)
    }
    if (length(extract_eigenvalue(ordination)[axes]) > 0 & !rescale) {
        strivar = as(c(p$label$x, p$label$y), "character")
        strivar = paste0(strivar, "   [", percvar, "%]")
        p = p + xlab(strivar[1]) + ylab(strivar[2])
    }
    return(p)
}

