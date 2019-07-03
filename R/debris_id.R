
#' Get pi
#'
#' This function calculates pi_low and pi_high, defined as the inner products of 
#' normalized gene expression and log fold changes for DE genes enriched in the low (background-enriched) 
#' and high (target-enriched) count droplets, respectively.
#'
#' @param x SCE. SCE object.
#'
#' @return An SCE object.
#' @importFrom Matrix Matrix t
#' @export
get_pi <- function(x){
	log2fc_l <- x@de@log2fc[x@de@deg_low]
	log2fc_h <- x@de@log2fc[x@de@deg_high]

	if (any(is.na(log2fc_l)) | any(is.na(log2fc_h))){
		stop("Fold changes cannot be NA")
	}
	if (any(is.infinite(log2fc_l)) | any(is.infinite(log2fc_h))){
		stop("Fold changes cannot be infinite")
	}

	genes_low <- intersect(rownames(x@diem@norm), x@de@deg_low)
	genes_high <- intersect(rownames(x@diem@norm), x@de@deg_high)

	pi_l <- t(x@diem@norm[genes_low,]) %*% log2fc_l[genes_low]
	pi_h <- t(x@diem@norm[genes_high,]) %*% -log2fc_h[genes_high]

	pi_df <- as.data.frame(cbind(pi_l, pi_h))
	rownames(pi_df) <- colnames(x@diem@norm)
	colnames(pi_df) <- c("pi_l", "pi_h")
	x@diem@pi <- pi_df

	x@dropl_info[,"pi_l"] <- rep(NA, nrow(x@dropl_info))
	x@dropl_info[,"pi_h"] <- rep(NA, nrow(x@dropl_info))

	x@dropl_info[rownames(pi_df),"pi_l"] <- pi_df[,"pi_l"]
	x@dropl_info[rownames(pi_df),"pi_h"] <- pi_df[,"pi_h"]
	return(x)
}

#' Initialize EM parameters for k=2 low and high groups
#'
#' Initialize the means and variances of the 2 Gaussians, as well as tau, 
#' the mixing parameter.
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
pi_init_tails <- function(x, pct=0.05){
	n <- floor(ncol(x@diem@counts) * pct)
	tc <- Matrix::colSums(x@diem@counts)[rownames(x@diem@pi)]
	topn <- rownames(x@diem@pi)[order(tc, decreasing=TRUE)[1:n]]
	botn <- rownames(x@diem@pi)[order(tc, decreasing=FALSE)[1:n]]
	mu_top <- apply(x@diem@pi[topn,], 2, mean)
	mu_bot <- apply(x@diem@pi[botn,], 2, mean)
	sgma_top <- apply(x@diem@pi[topn,], 2, var)
	sgma_bot <- apply(x@diem@pi[botn,], 2, var)
	mu0 <- list(mu_top, mu_bot)
	sgma0 <- list(sgma_top, sgma_bot)
	tau0 <- c(0.5, 0.5)
	ret <- list(mu0=mu0, sgma0=sgma0, tau0=tau0)
	return(ret)
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
#' score and the absolute value of the R coefficient must be greater than this threshold.
#' @param min_iter Integer. Minimum number of iterations.
#' @param max_iter Integer. Maximum number of iterations.
#' @param eps Numeric. Epsilon of log-likilhood change to call convergence.
#' @param verbose Boolean. Verbosity.
#'
#' @return SCE object with EM output in the \code{emo} slot. See \code{\link{run_mv_em_diag}} for details
#' @export
run_em <- function(x, 
				   pct_tails=0.05, 
				   min_iter=5, 
				   max_iter=1000, 
				   eps=1e-10, 
				   verbose=TRUE){
	runs <- list()

	if (length(x@diem@pi) == 0){
		stop("Calculate pi before running EM")
	}

	# Run EM
	if (verbose){
		cat("Running EM\n")
	}

	# Get starting parameters
	init <- pi_init_tails(x, pct=pct_tails)

	# Run EM
	x@diem@emo <- run_mv_em_diag(x@diem@pi, 
								 k=2, 
								 min_iter=min_iter, 
								 labels=x@diem@labels, 
								 max_iter=max_iter, 
								 eps=eps, 
								 verbose=FALSE,
								 mu0=init[["mu0"]], 
								 sgma0=init[["sgma0"]], 
								 tau0=init[["tau0"]])

	x@diem@emo@assign <- vector()
	x@diem@emo@assign["Target"] <- 1
	x@diem@emo@assign["Background"] <- 2

	# Flip if incorrectly assigned
	mean_tg <- colMeans(x@diem@pi[x@diem@emo@Z[,1] > 0.5,])
	mean_bg <- colMeans(x@diem@pi[x@diem@emo@Z[,2] > 0.5,])

	if ((mean_tg["pi_l"] > mean_bg["pi_l"]) & (mean_tg["pi_h"] < mean_bg["pi_h"])){
		x@diem@emo@assign["Target"] <- 2
		x@diem@emo@assign["Background"] <- 1
		cat("Warning: centers flipped after initialization\n")
	}
	
	if (verbose){
		cat("Finished EM\n")
	}

	return(x)
}


#' Call targets after EM
#'
#' Call targets from droplets if the log-likelihood membership probability is higher than \code{lk_fraction}.
#'
#' @param x SCE. SCE object.
#' @param lk_fraction Numeric. Select droplets that have a fraction of target likelihood greater than this number.
#'
#' @return SCE object
#' @export
call_targets <- function(x, lk_fraction=0.99){
	x@dropl_info[,"Call"] <- rep(NA, nrow(x@dropl_info))
	x@dropl_info[rownames(x@diem@emo@Z),"Call"] <- "Debris"
	target_names <- rownames(x@diem@emo@Z)[ x@diem@emo@Z[,x@diem@emo@assign["Target"]] > lk_fraction ]
	x@dropl_info[target_names,"Call"] <- "Nucleus"
	return(x)
}

#' Return target IDs
#'
#' @param x SCE. SCE object
#'
#' @return Character vector with droplet IDs of called targets
#' @export
targets_ids <- function(x){
	return(rownames(x@dropl_info)[x@dropl_info[,"Call"] %in% "Nucleus"])
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
