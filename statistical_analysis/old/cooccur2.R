############################################################################################
#
# This is a modified version of the CRAN package cooccur for metabarcoding projects
# i.e. large numbers of OTUs (1000s) and sites (100s)
# 5000 OTUs and 60 sites takes about 30 mins (single processor, uses <32G memory)
#
# There are afew bugs and limitations at the moment
# 1. BUG: If using thershold the numbers of removed species is not reported correctly
# 2. Combinations model is not implemented (needs some work to implement with site mask)
# 3. Only effects isn't implemented at the moment (this should be easy to add back in)
#############################################################################################

coprob2 <-
function(max_inc,j,min_inc,nsite){
    as.matrix(round(choose(max_inc,j) * choose(nsite - max_inc, min_inc - j),0) / round(choose(nsite,min_inc),0))
}

cooccur2 <- 
function (mat, type = "spp_site", thresh = TRUE, spp_names = FALSE,
    true_rand_classifier = 0.1, prob = "hyper", site_mask = NULL,
    only_effects = FALSE, eff_standard = TRUE, eff_matrix = FALSE)
{
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

	if (thresh == TRUE) {
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
	only_effects <- FALSE
	prob <- "hyper"	

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
		
	if (only_effects) {
		arr <- array(c(min_inc,(sp1_inc + sp2_inc)),c(nrow(nsite),ncol(nsite),3))
		prob_share_site <- apply(arr,1:2, function(x) {
			i <- x[1]-((x[2]-x[3])*(-1)+abs((x[2]-x[3])*(-1)))/2
			ii <- (i+abs(i))/2
			ii[is.na(ii)]<-0
			y<-rep(1,ii)
			return(y)
		})
	} else {
		if (prob == "hyper") {
			arr <- array(c(min_inc,nsite,max_inc),c(nrow(nsite),ncol(nsite),3))
			prob_share_site <- apply( arr , 1:2 , function(x) {
				x[is.na(x)]<-0
				y<-phyper(0:x[1],x[1],x[2]-x[1],x[3])
				y<-c(y[1],(y[-1]-y[-length(y)]))
				return(y)
			})
		} else if (prob == "comb") {		
			arr <- array(c(min_inc,nsite,(sp1_inc + sp2_inc),max_inc),c(nrow(nsite),ncol(nsite),4))
			prob_share_site <- apply(arr,1:2, function(x) {
				i <- x[1]-((x[2]-x[3])*(-1)+abs((x[2]-x[3])*(-1)))/2
				ii <- (i+abs(i))/2
				ii[is.na(ii)]<-0
				js <- sapply(ii,function(i)seq(1,i))
				y<-coprob(x[4],js,x[1],x[2])
				return(y)
			})			
		}

    	}		

	prob_share_site<- prob_share_site[which(lower.tri(prob_share_site))]

	p_lt <- sapply(seq(1,length(prob_share_site)),function(i) sum(unlist(prob_share_site[i])[1:(obs_cooccur[i]+1)]))
	p_gt <- sapply(seq(1,length(prob_share_site)),function(i) sum(unlist(prob_share_site[i])[(obs_cooccur[i]+1):length(unlist(prob_share_site[i]))]))
	p_exactly_obs <- sapply(seq(1,length(prob_share_site)),function(i) sum(unlist(prob_share_site[i])[(obs_cooccur[i]+1)]))

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
	if (thresh == TRUE) {
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
