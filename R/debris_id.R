
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

	pi_l <- t(x@diem@norm[x@de@deg_low,]) %*% log2fc_l
	pi_h <- t(x@diem@norm[x@de@deg_high,]) %*% -log2fc_h

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

#' Run EM on pi features
#' 
#' Run expectation maximization (EM) using diagonal covariance on pi matrix in an SCE object with k=2 groups. 
#' EM is run \code{n_runs} times with random initializations, and parameters with 
#' the best log likelihood are output. Within EM, iterate a minimum of \code{min_iter} times and 
#' a maximum of \code{max_iter} times. The \code{eps} value gives the threshold of percentage change in 
#' log likelihood until convergence is reached. \code{seedn} gives the seed for random number generation, 
#' which can be used to reproduce results.
#'
#' @param x SCE. SCE object with PCs and labels.
#' score and the absolute value of the R coefficient must be greater than this threshold.
#' @param min_iter Integer. Minimum number of iterations.
#' @param max_iter Integer. Maximum number of iterations.
#' @param eps Numeric. Epsilon of log-likilhood change to call convergence.
#' @param n_runs Integer. Number of EM initializations to run, with best result returned.
#' @param seedn Numeric. Seed for random number generation.
#' @param verbose Boolean. Verbosity.
#'
#' @return SCE object with EM output in the \code{emo} slot. See \code{\link{run_mv_em_diag}} for details
#' @export
run_em <- function(x, 
				   min_iter=5, 
				   max_iter=1000, 
				   eps=1e-10, 
				   n_runs=10, 
				   seedn=NULL, 
				   verbose=TRUE){
	runs <- list()

	if (length(x@diem@pi) == 0){
		stop("Calculate pi before running EM")
	}

	# Run EM
	if (verbose){
		cat("Running semi-supervised EM\n")
	}
	pb <- txtProgressBar(min = 0, max = n_runs, initial = 0, style = 3)
	for (i in 1:n_runs){
		runs[[i]] <- run_mv_em_diag(x@diem@pi, k=2, min_iter=min_iter, max_iter=max_iter,
									labels = x@diem@labels, eps=eps, seedn=seedn, verbose=FALSE)
		if (!is.null(seedn)) seedn <- seedn + 1
		setTxtProgressBar(pb,i)
	}
	close(pb)
	if (verbose){
		cat("Finished EM\n")
	}

	# Get run with max likelihood of parameters
	llks <- sapply(runs, function(x) x@llks[length(x@llks)])
	if (sum(!is.na(llks)) == 0){
		stop("No EM runs converged successfully.\n")
	}
	max_i <- which.max(llks)
	x@diem@emo <- runs[[max_i]]
	colnames(x@diem@emo@Z) <- c("Target", "Background")
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
call_targets <- function(x, lk_fraction=0.95){
	# if (is.null(lk_fraction)){
	# 	lk_fraction <- 1 - (0.05/nrow(x@diem@emo@Z))
	# }
	x@dropl_info[,"Call"] <- rep(NA, nrow(x@dropl_info))
	x@dropl_info[rownames(x@diem@emo@Z),"Call"] <- "Debris"
	target_names <- rownames(x@diem@emo@Z)[ x@diem@emo@Z[,"Target"] > lk_fraction ]
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
	return(rownames(x@dropl_info)[x@dropl_info[,"Call"]])
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
