
#' Get pi projection of normalized data
#'
#' This function calculates pi_low and pi_high, defined as the inner products of 
#' normalized gene expression and log fold changes for DE genes enriched in the low (background-enriched) 
#' and high (target-enriched) count droplets, respectively.
#'
#' @param x SCE. SCE object.
#'
#' @return An matrix of pi values
#' @importFrom Matrix Matrix t
#' @export
get_pi <- function(x, iteration=NULL){
	if (is.null(iteration)){
		iteration <- length(x@diem)
	} else {
		if (iteration > length(x@diem)) stop("iteration cannot be larger than the number of DIEM iterations run.")
	}

	if (length(x@norm) == 0){
		stop("Normalize counts data before calculating pi.")
	}

	de_obj <- x@diem[[iteration]]@de

	if (length(de_obj@log2fc) == 0){
		stop("No genes with log fold changes available, they have either not been calculated or no genes meet the criteria.")
	}

	log2fc_l <- de_obj@log2fc[de_obj@deg_low]
	log2fc_h <- de_obj@log2fc[de_obj@deg_high]

	if (any(is.na(log2fc_l)) | any(is.na(log2fc_h))){
		stop("Fold changes cannot be NA.")
	}
	if (any(is.infinite(log2fc_l)) | any(is.infinite(log2fc_h))){
		stop("Fold changes cannot be infinite.")
	}

	genes_low <- intersect(rownames(x@norm), de_obj@deg_low)
	genes_high <- intersect(rownames(x@norm), de_obj@deg_high)

	if ( (length(genes_low) == 0) | (length(genes_high) == 0) ){
		stop("No differential genes between background and signal.")
	}

	# log2fc_l is positive. log2fc_h is negative (so take negative)
	pi_l <- t(x@norm[genes_low,]) %*% log2fc_l[genes_low]
	pi_h <- t(x@norm[genes_high,]) %*% -log2fc_h[genes_high]

	pi_df <- as.data.frame(cbind(pi_l, pi_h))
	rownames(pi_df) <- colnames(x@norm)
	colnames(pi_df) <- c("pi_l", "pi_h")

	x@diem[[iteration]]@pi <- as.matrix(pi_df)

	return(x)
}

#' Initialize EM parameters for k=2 low and high groups
#'
#' Initialize the means, variances of the 2 Gaussians, as well as tau, 
#' the mixing parameter. Means and variances are intialized to the sample medians 
#' and variances of top 5% and bottom 5% (or by given \code{pct}).
#' 
#' @param x SCE. SCE object with pi calculated
#' @param pct Numeric. Percent of droplets at each tail to fix. A number from 0 to 1.
#'
#' @return A list with 3 components
#' \itemize{
#' \item{"mu0"}{"A list of size 2, with each element being a size 2 vector containing the means of pi_l and pi_h."}
#' \item{"sgma0"}{"A list of size 2, with each element being a size 2 vector containing the variances of pi_l and pi_h."}
#' \item{"tau0"}{"A numeric of size 2, containing the tau mixing coefficients."}
#'}
#' The first and second element of the lists of mu0 and sgma0 are for the high and low groups respectively
#' @importFrom Matrix colSums
#' @export
pi_init_tails <- function(x, pim, pct=0.05){
	if (pct <= 0 | pct > 0.5) pct <- 0.05
	
	if (length(pim) == 0){
		stop("Pi matrix is empty.")
	}

	n <- floor(nrow(pim) * pct)
	dc <- Matrix::colSums(x@counts[,rownames(pim)])
	topn <- names(dc)[order(dc, decreasing=TRUE)[1:n]]
	botn <- names(dc)[order(dc, decreasing=FALSE)[1:n]]

	mu_top <- apply(pim[topn,,drop=FALSE], 2, median)
	mu_bot <- apply(pim[botn,,drop=FALSE], 2, median)
	sgma_top <- cov(pim[topn,,drop=FALSE])
	sgma_bot <- cov(pim[botn,,drop=FALSE])
	# sgma_top <- apply(dat[topn,,drop=FALSE], 2, var)
	# sgma_bot <- apply(dat[botn,,drop=FALSE], 2, var)
	# mu0 <- list(mu_top, mu_bot)
	# sgma0 <- list(sgma_top, sgma_bot)
	mu0 <- rbind(mu_top, mu_bot)
	sgma0 <- rbind(sgma_top[lower.tri(sgma_top, diag=TRUE)],
				   sgma_bot[lower.tri(sgma_bot, diag=TRUE)])
	tau0 <- c(0.5, 0.5)
	ret <- list(Mu=mu0, LTSigma=sgma0, pi=tau0)
	return(ret)
}

#' x is a vector
#' @export
fraction_log <- function(x){
	x_c = x - max(x)
	x_c = exp(x_c)
	frac = x_c / sum(x_c);
	return(frac);
}

# mu1, mu2 - vectors of length 2
# Sigma1, Sigma2 - symmetric matrices of size 2 by 2
# Returns the Overlapping Coefficient (OVL) of the two distributions defined by (mu1, Sigma1) and (mu2, Sigma2)
#' @importFrom EMCluster dmvn
#' @importFrom pracma integral2
#' @export  
calc_ovl <- function(mu1, mu2, Sigma1, Sigma2){
	# define the integration boundaries based on 10 standard deviations of the larget element in Sigma1, Sigma2 and the largest component in mu1, mu2
	boundary <- max(max(abs(mu1)),max(abs(mu2))) + max(max(Sigma1),max(Sigma2))*10
	min.mvns <- function(x, mu1, mu2, Sigma1, Sigma2){
		f1.x <- EMCluster::dmvn(x=x, mu1, Sigma1, log=FALSE)
		f2.x <- EMCluster::dmvn(x=x, mu2, Sigma2, log=FALSE)
		return(min(f1.x, f2.x))
	}
	f <- function(x,y) (min.mvns(c(x,y), mu1 = mu1, mu2 = mu2, Sigma1 = Sigma1, Sigma2 = Sigma2))
	ovl <- pracma::integral2(fun = f, -boundary, boundary, -boundary, boundary, vectorized = FALSE)
	return(ovl$Q)
}

lower_to_full <- function(x){
	lx <- length(x)
	p <- (-1+sqrt(1+(8*lx)))/2
	y <- matrix(nrow=p, ncol=p)
	y[lower.tri(y, diag=TRUE)] <- x
	y[upper.tri(y)] <- t(y)[upper.tri(y)]
	return(y)
}

#' Check initial separation of tails
#'
#'
#' @export
check_init_sep <- function(x, pct=0.05){
	if (pct <= 0 | pct > 0.5) pct <- 0.05

	pim <- x@diem[[1]]@pi
	
	if (length(pim) == 0){
		stop("Pi matrix is empty.")
	}

	n <- floor(nrow(pim) * pct)
	dc <- Matrix::colSums(x@counts[,rownames(pim)])
	topn <- names(dc)[order(dc, decreasing=TRUE)[1:n]]
	botn <- names(dc)[order(dc, decreasing=FALSE)[1:n]]

	mu_top <- apply(pim[topn,,drop=FALSE], 2, median)
	mu_bot <- apply(pim[botn,,drop=FALSE], 2, median)
	sgma_top <- cov(pim[topn,,drop=FALSE])
	sgma_bot <- cov(pim[botn,,drop=FALSE])
	
	mu <- rbind(mu_top, mu_bot)
	sgma <- rbind(sgma_top[lower.tri(sgma_top, diag=TRUE)],
				  sgma_bot[lower.tri(sgma_bot, diag=TRUE)])

	ovl <- calc_ovl(mu1=mu[1,], mu2=mu[2,], Sigma1=sgma[1,], Sigma2=sgma[2,])

	cat(paste0("Overlap of +/- 5% tails is ", as.character(round(ovl, 2)*100), "%.\n"))
}

#' Run EM on pi features
#' 
#' Run expectation maximization (EM) using diagonal covariance on pi matrix in an SCE object with k=2 groups. 
#' EM is run \code{n_runs} times with random initializations, and parameters with 
#' the best log likelihood are output. Within EM, iterate a minimum of \code{min_iter} times and 
#' a maximum of \code{max_iter} times. The \code{eps} value gives the threshold of percentage change in 
#' log likelihood until convergence is reached.
#'
#' @param x SCE. SCE object with PCs.
#' @param verbose Boolean. Verbosity.
#'
#' @return SCE object with EM output in the \code{emo} slot. See \code{\link{run_mv_em_diag}} for details
#' @importFrom EMCluster emcluster
#' @importFrom mvtnorm dmvnorm
#' @export
run_em <- function(x, init_em_tails=0.05, iteration=NULL, verbose=TRUE){
	if (is.null(iteration)){
		iteration <- length(x@diem)
	} else {
		if (iteration > length(x@diem)) stop("iteration cannot be larger than the number of DIEM iterations run.\n")
	}

	pim <- x@diem[[iteration]]@pi

	if (length(pim) == 0){
		stop("Calculate pi before running EM")
	}

	# Run EM
	eminit <- pi_init_tails(x, pim, pct=init_em_tails)
	require(EMCluster)
	if (sum(x@labels) > 0){
		if (verbose) cat("Running semi-supervised EM\n")
		emo <- emcluster(pim, eminit, lab=x@labels)
	} else {
		if (verbose) cat("Running unsupervised EM\n")
		emo <- emcluster(pim, eminit)
	}
	emo <- list(Mu=emo$Mu,
				LTSigma=emo$LTSigma, 
				pi=emo$pi, 
				llhdval=emo$llhdval, 
				conv.iter=emo$conv.iter)

	# Get likelihoods
	Signal <- apply(pim, 1, function(a){
					k=1
					mvtnorm::dmvnorm(x=a, 
					mean=emo$Mu[k,], 
					sigma=lower_to_full(emo$LTSigma[k,]), 
					log=TRUE) + log(emo$pi[k])
					})
	names(Signal) <- rownames(pim)
	Background <- apply(pim, 1, function(a){
					k=2
					mvtnorm::dmvnorm(x=a, 
					mean=emo$Mu[k,], 
					sigma=lower_to_full(emo$LTSigma[k,]), 
					log=TRUE) + log(emo$pi[k])
					})
	names(Background) <- rownames(pim)

	z <- cbind(Signal, Background)
	zp <- t(apply(z, 1, fraction_log)) # Get the responsibilities (posterior probabilities)
	emo$Z <- as.matrix(zp)

	# Overlap of estimated MVNs
	ovl <- calc_ovl(mu1=emo$Mu[1,], mu2=emo$Mu[2,], Sigma1=emo$LTSigma[1,], Sigma2=emo$LTSigma[2,])

	if (ovl > 0.25){
		cat("========================================================================\n")
		cat(paste0("Warning: Overlap of estimated distributions is ", 
					as.character(round(ovl, 2)*100), "%.\n"))
		cat("         Check if data separates into clusters\n")
		cat("========================================================================\n")
	}

	emo$ovl <- ovl

	emo$assign <- c("Signal"=1, "Background"=2)

	if (emo$Mu[emo$assign["Signal"],1] > emo$Mu[emo$assign["Background"],1] & 
		emo$Mu[emo$assign["Signal"],2] < emo$Mu[emo$assign["Background"],2]){
		cat("========================================================================\n")
		cat("Warning: Centers of Signal and Background have flipped.\n")
		cat("         Flipping labels\n")
		cat("         Check if data separates into clusters, and/or try running as semi-supervised.\n")
		cat("========================================================================\n")
		emo$assign <- emo$assign[c(2,1)]
	} else if (emo$Mu[emo$assign["Signal"],1] > emo$Mu[emo$assign["Background"],1]){
		cat("========================================================================\n")
		cat("Warning: pi_low mu of signal is greater than that of background\n")
		cat("         Check if data separates into clusters\n")
		cat("========================================================================\n")
	} else if (emo$Mu[emo$assign["Signal"],2] < emo$Mu[emo$assign["Background"],2]){
		cat("========================================================================\n")
		cat("Warning: pi_high mu of signal is less than that of background\n")
		cat("         Check if data separates into clusters\n")
		cat("========================================================================\n")
	}

	# Check if any signal calls are below the background center
	signal_cells <- rownames(zp[zp[,emo$assign["Signal"]] > x@p_thresh,])
	pi_signal <- pim[signal_cells,]
	pi_lb <- pi_signal[,"pi_l"] > emo$Mu[emo$assign["Background"],1]
	pi_hb <- pi_signal[,"pi_h"] < emo$Mu[emo$assign["Background"],2]
	pi_crossd <- pi_lb & pi_hb
	if (any(pi_crossd)){
		cat("========================================================================\n")
		cat("Warning: Tail density of signal crossed over and above density of background.\n")
		cat("         Misclassification of background as signal!!!\n")
		cat("         Check if data separates into clusters\n")
		cat("========================================================================\n")
	}

	x@diem[[iteration]]@emo <- emo
	if (verbose){
		cat("Finished EM\n")
	}
	return(x)
}

#' Call targets after EM
#'
#' Call targets from droplets if the log-likelihood membership 
#' probability is higher than \code{lk_fraction} given during 
#' initialization.
#'
#' @param x SCE. SCE object.
#' @param interation Integer. DIEM iteration number to use for calliing targets
#'
#' @return SCE object
#' @export
call_targets <- function(x, iteration=NULL){
	if (is.null(iteration)){
		iteration <- length(x@diem)
	} else {
		if (iteration > length(x@diem)) stop("iteration cannot be larger than the number of DIEM iterations run.\n")
	}

	Z <- x@diem[[iteration]]@emo[["Z"]]
	asgn <- x@diem[[iteration]]@emo[["assign"]]

	calls <- rep("Debris", nrow(Z))
	names(calls) <- rownames(Z)

	tb <- Z[,asgn["Signal"]] > x@p_thresh
	target_names <- rownames(Z)[tb]
	calls[tb] <- "Signal"

	x@diem[[iteration]]@calls <- calls

	return(x)
}

#' Return target IDs
#'
#' @param x SCE. SCE object
#'
#' @return Character vector with droplet IDs of called targets
#' @export
targets_ids <- function(x, iteration=NULL){
	if (is.null(iteration)){
		iteration <- length(x@diem)
	} else {
		if (iteration > length(x@diem)) stop("iteration cannot be larger than the number of DIEM iterations run.\n")
	}

	calls <- x@diem[[iteration]]@calls

	return(names(calls)[calls == "Signal"])
}

#' Insert pi and calls to dropl_info
#'
#' @param x SCE. SCE object
#'
#' @return Character vector with droplet IDs of called targets
#' @export
fill_dropl_info <- function(x, iteration=NULL){
	if (is.null(iteration)){
		iteration <- length(x@diem)
	} else {
		if (iteration > length(x@diem)) stop("iteration cannot be larger than the number of DIEM iterations run.\n")
	}

	# Fill calls
	x@dropl_info[,"Call"] <- rep(NA, nrow(x@dropl_info))
	calls <- x@diem[[iteration]]@calls
	x@dropl_info[names(calls),"Call"] <- calls

	# Fill pi
	pi_df <- as.data.frame(x@diem[[iteration]]@pi)
	pi_df[,"pi"] <- pi_df[,"pi_l"] - pi_df[,"pi_h"]
	x@dropl_info[,"pi_l"] <- rep(NA, nrow(x@dropl_info))
	x@dropl_info[,"pi_h"] <- rep(NA, nrow(x@dropl_info))
	x@dropl_info[rownames(pi_df),"pi_l"] <- pi_df[,"pi_l"]
	x@dropl_info[rownames(pi_df),"pi_h"] <- pi_df[,"pi_h"]
	x@dropl_info[rownames(pi_df),"pi"] <- pi_df[,"pi"]

	return(x)

}

#' Get percent of reads aligning to MT genome and MALAT1 gene
#'
#' Places percent of reads aligning to MT genome and MALAT1 gene in 
#' \code{x@dropl_info}, under columns MALAT1 and MT_PCT. The gene 
#' names in the rows of \code{x@counts} should be MALAT1 and only mitochondrial 
#' genes should start with MT-.
#'
#' @param x SCE.
#'
#' @return SCE object
#' @export
get_mt_malat1 <- function(x){
	n_drop <- nrow(x@dropl_info)
	
	malat_gene <- grep("^malat1", rownames(x@counts), ignore.case=T, value=TRUE)
	if ( length(malat_gene) == 0 ) malat1 <- NA
	else malat1 <- x@counts[malat_gene,] / x@dropl_info[,"total_counts"]

	mt_genes <- grep("^mt-", rownames(x@counts), ignore.case=T, value=TRUE)
	mt_pct <- Matrix::colSums( x@counts[mt_genes,,drop=FALSE] ) / x@dropl_info[,"total_counts"]

	x@dropl_info[,"MALAT1"] <- malat1
	x@dropl_info[,"MT_PCT"] <- mt_pct
	return(x)
}

#' Get percent of reads aligning to given genes
#'
#' @param x SCE.
#' @param genes Character. Genes to calculate percentage of in counts
#' @param name Character. Column name to place in dropl_info
#'
#' @return SCE object
#' @export
get_gene_pct <- function(x, genes, name){
	n_drop <- nrow(x@dropl_info)
	expr <- x@counts[genes,]
	if (length(expr) == 0){
		stop("None of genes found in counts.")
	}
	gene_pct <- Matrix::colSums(x@counts[genes,]) / x@dropl_info[,"total_counts"]
	x@dropl_info[names(gene_pct),] <- gene_pct
	return(x)
}
