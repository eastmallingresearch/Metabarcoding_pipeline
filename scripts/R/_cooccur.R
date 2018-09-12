########################################################################
#
# This is a modified (for speed) version of the CRAN package cooccur
# only hyper is currently available for prob calculation
#
########################################################################

_coprob2 <-
function(max_inc,j,min_inc,nsite){
    require(gmp)
    as.matrix(round(chooseZ(max_inc,j) * chooseZ(nsite - max_inc, min_inc - j),0) / round(chooseZ(nsite,min_inc),0))
}

_cooccur2 <- 
function (mat, type = "spp_site", thresh = TRUE, spp_names = FALSE,
    true_rand_classifier = 0.1, prob = "hyper", site_mask = NULL,
    only_effects = FALSE, eff_standard = TRUE, eff_matrix = FALSE)
{
    require(cooccur)
    if (type == "spp_site") {
        spp_site_mat <- mat
    }
    if (type == "site_spp") {
        spp_site_mat <- t(mat)
    }
    if (spp_names == TRUE) {
        spp_key <- data.frame(num = 1:nrow(spp_site_mat), spp = row.names(spp_site_mat))
    }
    if (!is.null(site_mask)) {
        if (nrow(site_mask) == nrow(spp_site_mat) & ncol(site_mask) ==
            ncol(spp_site_mat)) {
            N_matrix <- create.N.matrix(site_mask)
        }
        else {
            stop("Incorrect dimensions for site_mask, aborting.")
        }
    }
    else {
        site_mask <- matrix(data = 1, nrow = nrow(spp_site_mat),
            ncol = ncol(spp_site_mat))
        N_matrix <- matrix(data = ncol(spp_site_mat), nrow = nrow(spp_site_mat),
            ncol = nrow(spp_site_mat))
    }
    spp_site_mat[spp_site_mat > 0] <- 1
    tsites <- ncol(spp_site_mat)
    nspp <- nrow(spp_site_mat)
    spp_pairs <- choose(nspp, 2)
    incidence <- prob_occur <- obs_cooccur <- prob_cooccur <- exp_cooccur <- matrix(nrow = spp_pairs,
        ncol = 3)
    incidence <- prob_occur <- matrix(nrow = nrow(N_matrix),
        ncol = ncol(N_matrix))

### hacked bit - need updating to add back only_effects=T (easy) and combinations probablity (not so easy, my current method requires huge amounts of memory) 

	mat_mask <- as.matrix(mat*site_mask)
	incidence <- t(apply(mat_mask,1,function(m) site_mask%*%m))
	pairs <- t(apply(mat_mask,1,function(m) mat_mask%*%m))
	diag(incidence) <- NA
	prob_occur <- incidence/N_matrix

	obs_cooccur <- pairs
	prob_cooccur <- prob_occur*t(prob_occur)
	exp_cooccur <- prob_cooccur*N_matrix

	if (thresh) {
		n_pairs <- sum(prob_cooccur>=0,na.rm=T)/2
		t_table <- exp_cooccur>=1
		prob_cooccur <- prob_cooccur*t_table
		obs_cooccur <- obs_cooccur*t_table
		exp_cooccur <- exp_cooccur*t_table
		n_omitted <- n_pairs - sum(t_table,na.rm=T)
	}

	sp1_inc=incidence
	sp2_inc=t(incidence)
	max_inc <- pmax(sp1_inc, sp2_inc)
	min_inc <- pmin(sp1_inc, sp2_inc)
	nsite <- N_matrix
	psite <- nsite + 1

	arr <- array(c(min_inc,nsite,max_inc,(sp1_inc + sp2_inc)),c(nrow(nsite),ncol(nsite),4))

	effect_func <- function(a=arr){
		X <- apply(a,1:2, function(x) {
			x[is.na(x)]<-0
			i <- x[1]-((x[2]-x[4])*(-1)+abs((x[2]-x[4])*(-1)))/2
			ii <- (i+abs(i))/2
			#ii[is.na(ii)]<-0
			y<-rep(1,ii)
			return(y)
		})
	}
	
	hyper_func <- function(a=arr){
		X <- apply(a, 1:2 , function(x) {
			x[is.na(x)]<-0
			y<-phyper(0:x[1],x[1],x[2]-x[1],x[3])
			y<-c(y[1],(y[-1]-y[-length(y)]))
			return(y)
		})
		return(X)
	}

	comb_func2  <- function(a=arr){
		y <- sapply(0:max(a[,,2]),function(i) coprob2(arr[,,3],i,arr[,,1],arr[,,2]))
		a[is.na(a)]<-0
		start <- a[,,4]-a[,,2]			
		i <- a[,,1]-((a[,,2]-a[,,4])*(-1)+abs((a[,,2]-a[,,4])*(-1)))/2
		len <- (i+abs(i))/2
		return(y[start:(start+len)])		
	}

	comb_func  <- function(a=arr){
		X <- apply(a,1:2, function(x) {
			x[is.na(x)]<-0
			y<-coprob(x[3],0:x[2],x[1],x[2])
			start <- ((x[4]-x[2])+abs(x[4]-x[2]))/2			
			#i <- x[1]-((x[2]-x[4])*(-1)+abs((x[2]-x[4])*(-1)))/2
			#len <- x[1]-start#(i+abs(i))/2
			#ii[is.na(ii)]<-0
			return(y[start:(x[1]+1)])
		})
		return(X)
	}

	if (only_effects) {
		prob_share_site <- effect_func()
	} else {
		if (prob == "hyper") {
			prob_share_site <- hyper_func()
		} else if (prob == "comb") {		
			prob_share_site <- comb_func()			
		} else {
			print("Unsupported probability model specified\n. Returning effect sizes only")
			prob_share_site<-effect_func()
			only_effect=T
		}
    	}
#return(prob_share_site)
	prob_share_site<- prob_share_site[which(lower.tri(prob_share_site))]
	obs_cooccur<- obs_cooccur[which(lower.tri(obs_cooccur))]
	prob_cooccur<- prob_cooccur[which(lower.tri(prob_cooccur))]
	exp_cooccur<- exp_cooccur[which(lower.tri(exp_cooccur))]
	t_table <- t_table[which(lower.tri(t_table))]

	sp <- matrix(rep(seq(1,nspp),nspp),nrow=nspp,ncol=nspp)
	sp1 <- t(sp)[which(lower.tri(sp[-nrow(sp),-ncol(sp)],diag=T))]
	sp2 <- sp[which(lower.tri(sp,diag=F))]

	sp <- matrix(rep(rownames(sp1_inc),ncol(sp1_inc)),nrow=ncol(sp1_inc),ncol=ncol(sp1_inc))
	sp1_name <- t(sp)[which(lower.tri(sp[-nrow(sp),-ncol(sp)],diag=T))]
	sp2_name <- sp[which(lower.tri(sp,diag=F))]

	sp1_inc<- sp1_inc[which(lower.tri(sp1_inc))]
	sp2_inc<- sp2_inc[which(lower.tri(sp2_inc))]

	p_lt <- sapply(seq(1,length(prob_share_site)),function(i) sum(unlist(prob_share_site[i])[1:(obs_cooccur[i]+1)],na.rm=T))
	p_gt <- sapply(seq(1,length(prob_share_site)),function(i) sum(unlist(prob_share_site[i])[(obs_cooccur[i]+1):length(unlist(prob_share_site[i]))],na.rm=T))
	p_exactly_obs <- sapply(seq(1,length(prob_share_site)),function(i) sum(unlist(prob_share_site[i])[(obs_cooccur[i]+1)],na.rm=T))

	p_lt <- round(p_lt, 5)
	p_gt <- round(p_gt, 5)
	p_exactly_obs <- round(p_exactly_obs, 5)

	prob_cooccur <- round(prob_cooccur, 3)
	exp_cooccur <- round(exp_cooccur, 1)

	output<-data.frame(
		sp1=sp1,
		sp2=sp2,
		sp1_inc=sp2_inc,
		sp2_inc=sp1_inc,
		obs_cooccur=obs_cooccur,
		prob_cooccur=prob_cooccur,
		exp_cooccur=exp_cooccur,
		p_lt=p_lt,
		p_gt=p_gt,
		sp1_name=sp1_name,
		sp2_name=sp2_name
	)
	if (thresh) {
		# could do this earlier and save a couple of minutes execution time
		output <- output[t_table,]
	}
####

    true_rand <- (nrow(output[(output$p_gt >= 0.05 & output$p_lt >=
        0.05) & (abs(output$obs_cooccur - output$exp_cooccur) <=
        (tsites * true_rand_classifier)), ]))
    output_list <- list(call = match.call(), results = output,
        positive = nrow(output[output$p_gt < 0.05, ]), negative = nrow(output[output$p_lt <
            0.05, ]), co_occurrences = (nrow(output[output$p_gt <
            0.05 | output$p_lt < 0.05, ])), pairs = nrow(output),
        random = true_rand, unclassifiable = nrow(output) - (true_rand +
            nrow(output[output$p_gt < 0.05, ]) + nrow(output[output$p_lt <
            0.05, ])), sites = N_matrix, species = nspp, percent_sig = (((nrow(output[output$p_gt <
            0.05 | output$p_lt < 0.05, ])))/(nrow(output))) *
            100, true_rand_classifier = true_rand_classifier)
    if (spp_names == TRUE) {
        output_list$spp_key <- spp_key
        output_list$spp.names = row.names(spp_site_mat)
    }
    else {
        output_list$spp.names = c(1:nrow(spp_site_mat))
    }
    if (thresh == TRUE) {
        output_list$omitted <- n_omitted
        output_list$pot_pairs <- n_pairs
    }
    class(output_list) <- "cooccur"
    if (only_effects == F) {
        output_list
    }
    else {
        effect.sizes(mod = output_list, standardized = eff_standard,
            matrix = eff_matrix)
    }
}
