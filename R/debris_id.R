

#' Get PCs from expression data
#'
#' @param x SCE. An SCE object
#' @param n_pcs Integer. Number of PCs to calculate
#' @param center Boolean. Whether to center variables
#' @param scale. Boolean. Whether to scale variables to have unit variance
#'
#' @return An SCE object with \code{\link[irlba]{prcomp_irlba}} output in 
#' the \code{pcs} slot
#' @importFrom irlba prcomp_irlba
#' @export
get_pcs <- function(x, n_pcs=15, center = TRUE, scale. = TRUE){
	pcs <- prcomp_irlba(Matrix::t(x@norm[x@vg,]), 
						n=n_pcs, 
						retx=TRUE, 
						center=center, 
						scale.=scale.)
	rownames(pcs$x) <- colnames(x@norm)
	rownames(pcs$rotation) <- x@vg
	colnames(pcs$x) <- paste0("PC", as.character(1:n_pcs))
	x@pcs <- pcs
	return(x)
}

#' Specify background droplets and candidate droplets
#' 
#' This function returns the barcodes that are considered originating from background, and those that are
#' considered candidate targets (cells/nuclei). Barcodes are considered background if either they fall below
#' the top \code{top_n_cand} barcodes, or if the number of reads is less than bg_max.
#' 
#' total read/UMI count is greater than or equal to \code{bg_min} and less than or equal to \code{bg_max}. 
#' Candidate cells/nuclei are the top \code{top_n_cand} droplets ranked by read counts. This should be 
#' at least the expected number of targets recovered from the experiment. If the 
#' top \code{top_n_cand} droplets includes those with read/UMI counts < \code{bg_max}, the \code{bg_max} 
#' parameter is lowered so that \code{top_n_cand} candidates are always output.
#'
#' @param x SCE. An SCE object
#' @param bg_max Numeric. Maximum number of counts for a droplet to be considered background. Ignored if top_n_cand is set.
#' @param top_n_cand Numeric. Specifies the top number of droplets by expression to considered candidates
#'
#' @return SCE object with bg_info filled, a list with components
#' \item{Background}{Character vector of column names of x that are considered background droplets}
#' \item{Candidate}{Character vector of column names of x that are considered candidate droplets}
#' \item{Labels}{Numeric vector for droplets specifying 0 for candidate and 1 for background}
#' @export
get_bg_cand <- function(x, bg_max=NULL, top_n_cand=NULL){
	drplt_ids <- colnames(x@norm)
	labels <- rep(1, length(drplt_ids))
	names(labels) <- drplt_ids

	if (!is.null(top_n_cand)){
		top_n_ix <- order(x@dropl_counts, na.last = TRUE, decreasing = TRUE)[1:top_n_cand]
		labels[top_n_ix] <- 0
	} else if (!is.null(bg_max)){
		labels[x@dropl_counts > bg_max] <- 0
	}
	else {
		stop("Either top_n_cand or bg_max must be set")
	}
	cand_ids <- names(labels)[labels == 0]
	bg_ids <- names(labels)[labels == 1]

	x@bg_info <- list(Background = bg_ids, Candidate = cand_ids, Labels = labels)
	return(x)
}

#' Run EM on PCs
#' 
#' Run expectation maximization (EM) using diagonal covariance on expression PCs in an SCE object with k=2 groups. 
#' EM is run \code{n_runs} times with random initializations, and the EM output parameters with 
#' the best log likelihood is output. Within EM, iterate a minimum of \code{min_iter} times and 
#' a maximum of \code{max_iter} times. The \code{eps} value gives the threshold of percentage change in 
#' log likelihood until covnergence is called. \code{seedn} gives the seed for random number generation, 
#' which can be used to reproduce results.
#'
#' @param x SCE. SCE object with PCs and labels
#' @param n_pcs Integer. Number of expression PCs to use for running EM
#' @param min_iter Integer. Minimum number of iterations
#' @param max_iter Integer. Maximum number of iterations
#' @param eps Numeric. Epsilon of log-likilhood change to call convergence
#' @param n_runs Integer. Number of EM runs, with best result returned
#' @param seedn Numeric. Seed for random number generation
#' @param verbose Boolean. Print out logging information
#'
#' @return SCE object with EM output in the \code{emo} slot. See \code{\link{run_mv_em_diag}} for details
#' @export
run_em_pcs <- function(x, 
					   n_pcs=NULL, 
					   min_iter=5, 
					   max_iter=1000, 
					   eps=1e-10, 
					   n_runs=5, 
					   seedn=NULL, 
					   verbose=TRUE){
	runs <- list()
	for (i in 1:n_runs){
		if (!is.null(n_pcs)){
			runs[[i]] <- run_mv_em_diag(x@pcs$x[,1:n_pcs], k=2, min_iter=min_iter, max_iter=max_iter,
										labels = x@bg_info$Labels, eps=eps, seedn=seedn, verbose=verbose)
		} else {
			runs[[i]] <- run_mv_em_diag(x@pcs$x, k=2, min_iter=min_iter, max_iter=max_iter, 
										labels = x@bg_info$Labels, eps=eps, seedn=seedn, verbose=verbose)
		}
		if (!is.null(seedn)) seedn <- seedn + 1
	}
	llks <- sapply(runs, function(x) x$llks[length(x$llks)])
	max_i <- which.max(llks)
	x@emo <- runs[[max_i]]
	# Add column names
	colnames(x@emo$Z) <- c("Background", "Target")
	return(x)
}

#' Call targets after EM
#'
#' Call targets, typically cells or nuclei, from droplets if the 
#' log-likelihood membership probability is higher than \code{p}.
#'
#' @param x SCE. SCE object with EM data
#' @param min_count Integer. Minimum number of read/UMI counts to call a target
#' @param min_gene Integer. Minimum number of genes detected to call a target
#' @param p Numeric. Proportion of likelihood to call a target
#'
#' @return SCE object with \code{targets} identified
#' @export
call_targets <- function(x, min_count=200, min_gene=200, p=0.95){
	keep <- (x@emo$Z[,2] > p) & (x@dropl_counts >= min_count) & (x@n_genes >= min_gene)
	x@targets <- rownames(x@emo$Z)[keep]
	return(x)
} 

#' Summarize results
#' 
#' Create a data frame summarizing summary stats for target and background calls.
#'
#' @param x SCE. SCE object with PCs and labels
#'
#' @return SCE object with summary stats in \code{other_info}
#' @export
summary_results <- function(x){
	n_drop <- length(x@dropl_counts)
	target_call <- integer(n_drop); names(target_call) <- colnames(x@raw); target_call[x@targets] <- 1
	malat_gene <- grep("malat1", rownames(x@raw), ignore.case=T, value=TRUE)
	if ( length(malat_gene) == 0 ) malat1 <- NA
	else malat1 <- x@raw[malat_gene,] / x@dropl_counts
	mt_genes <- grep("MT-", rownames(x@raw), ignore.case=T, value=TRUE)
	mt_pct <- Matrix::colSums( x@raw[mt_genes,] ) / x@dropl_counts
	x@other_data <- data.frame(Call=target_call, n_counts=x@dropl_counts, n_genes=x@n_genes, MALAT1=malat1, MT_PCT=mt_pct, x@pcs$x)
	return(x)
}

#' Run diem pipeline directly from 10X output
#'
#' @param path_10X Character. Path to 10X output, containing genes, barcodes, and expression files
#' @param min_bg_count Numeric. Minimum number of counts for a droplet to be considered background
#' @param max_bg_count Numeric. Maximum number of counts for a droplet to be considered background
#' @param n_pcs Integer. Number of expression PCs to use for running EM
#' @param n_runs Integer. Number of EM runs, with best result returned
#' @param expected_targets Integer. Number of expected targets
#' @param min_count Integer. Minimum number of counts to call a target
#' @param p Numeric. Minimum membership probability to call a target
#' @param seedn Numeric. Seed for random number generation
#' @param verbose Boolean. Print out logging information
#' 
#' @useDynLib diem
#' @export
DIEM_10X <- function(path_10X, min_bg_count=30, max_bg_count=150, n_pcs=5, n_runs=5, expected_targets=NULL, min_count=200, p=0.95, seedn=NULL, verbose=TRUE){
	x <- read_10x(path_10X)
	sce <- create_SCE(x)
	sce <- subset_dropls(sce, min_c=min_bg_count)
	sce <- normalize(sce)
	sce <- get_var_genes(sce)
	sce <- get_pcs(sce, n_pcs)
	sce <- get_bg_cand(sce, bg_max=max_bg_count)
	sce <- run_em_pcs(sce, n_pcs=n_pcs, n_runs=n_runs, seedn=seedn, verbose=verbose)
	sce <- call_targets(sce, min_count=min_count, p=p)
	sce <- summary_results(sce)
	return(sce)
}


#' Run diem pipeline for 10X samples
#'
#' Run diem pipeline on 10X samples, returning \code{SCE} objects.
#'
#' @param path_10X Character. Path to 10X output, containing genes, barcodes, and expression files
#' @param min_bg_count Numeric. Minimum number of counts for a droplet to be considered background
#' @param max_bg_count Numeric. Maximum number of counts for a droplet to be considered background
#' @param n_pcs Integer. Number of expression PCs to use for running EM
#' @param n_runs Integer. Number of EM runs, with best result returned
#' @param expected_targets Integer. Number of expected targets
#' @param min_count Integer. Minimum number of counts to call a target
#' @param p Numeric. Minimum membership probability to call a target
#' @param seedn Numeric. Seed for random number generation
#' @param verbose Boolean. Print out logging information
#' 
#' @export
run_diem_10x <- function(path_10X, min_bg_count=30, max_bg_count=150, n_pcs=5, n_runs=5, expected_targets=NULL, min_count=200, p=0.95, seedn=NULL, verbose=TRUE){
	sce_list <- list()
	for (path in path_10X){
		if (verbose) cat(paste0("Running ", path, "\n"))
		sce_list[[path]] <- DIEM_10X(path, min_bg_count=min_bg_count, max_bg_count=max_bg_count, n_pcs=n_pcs, n_runs=n_runs, expected_targets=expected_targets, min_count=min_count, p=p, seedn=seedn, verbose=verbose)
	}
	return(sce_list)
}
