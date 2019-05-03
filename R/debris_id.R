

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
#' considered candidate targets (cells/nuclei). Barcodes are considered background if their 
#' total read/UMI count is greater than or equal to \code{bg_min} and less than or equal to \code{bg_max}. 
#' Candidate cells/nuclei are the top \code{top_n_cand} droplets ranked by read counts. If the 
#' top \code{top_n_cand} droplets includes those with read/UMI counts < \code{bg_max}, the \code{bg_max} 
#' parameter is lowered so that \code{top_n_cand} candidates are always output.
#'
#' @param x SCE. An SCE object
#' @param bg_min Numeric. Minimum number of counts for a droplet to be considered background
#' @param bg_max Numeric. Maximum number of counts for a droplet to be considered background
#' @param top_n_cand Numeric. Specifies the top number of droplets by expression to considered candidates
#'
#' @return SCE object with bg_info filled, a list with components
#' \item{Background}{Character vector of column names of x that are considered background droplets}
#' \item{Candidate}{Character vector of column names of x that are considered candidate droplets}
#' \item{Labels}{Numeric vector for droplets specifying 0 for candidate and 1 for background}
#' @export
get_bg_cand <- function(x, bg_min, bg_max, top_n_cand){
	if (top_n_cand < 1) stop(paste0("top_n_cand must be a positive integer, not ", str(top_n_cand)))
	if (bg_min < 1) stop(paste0("bg_min must be a positive integer, not ", str(bg_min)))
	if (bg_max < 1) stop(paste0("bg_max must be a positive integer, not ", str(bg_max)))
	if (bg_min > bg_max) stop("bg_min must be less than bg_max")
	
	drplt_ids <- colnames(x@norm)
	bg_ids <- drplt_ids[ x@dropl_counts >= bg_min & x@dropl_counts <= bg_max ]
	top_n_ix <- order(x@dropl_counts, na.last = TRUE, decreasing = TRUE)[1:top_n_cand]
	cand_ids <- drplt_ids[top_n_ix]
	n_overlap <- length(intersect(bg_ids, cand_ids))
	bg_max_new <- x@dropl_counts[top_n_cand]
	if (n_overlap > 0){
		cat(paste0("Warning: ", as.character(n_overlap), " candidates have counts <= bg_max.\n"))
		cat(paste0("Lowering bg_max threshold to ", as.character(bg_max_new), "\n"))
	}

	labels <- integer(length(drplt_ids))
	names(labels) <- drplt_ids
	labels[bg_ids] <- 1
	labels[cand_ids] <- 0
	bg_ids <- names(labels)[labels == 1] # reset if candidates overlapped background IDs.
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
					   n_runs=10, 
					   seedn=NULL, 
					   verbose=TRUE){
	runs <- list()
	for (i in 1:n_runs){
		if (!is.null(n_pcs)){
			runs[[i]] <- run_mv_em_diag(sce@pcs$x[,1:n_pcs], k=2, min_iter=min_iter, max_iter=max_iter,
										labels = x@bg_info$Labels, eps=eps, seedn=seedn, verbose=verbose)
		} else {
			runs[[i]] <- run_mv_em_diag(sce@pcs$x, k=2, min_iter=min_iter, max_iter=max_iter, 
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
#' @param p Numeric. Proportion of likelihood to call a target
#'
#' @return SCE object with \code{targets} identified
#' @export
call_targets <- function(x, min_count=200, p=0.5){
	keep <- (x@emo$Z[,2] > p) & (x@dropl_counts >= min_count)
	x@targets <- rownames(x@emo$Z)[keep]
	return(x)
} 

#' Run DIEM pipeline directly from 10X output
#'
DIEM_10X <- function(path_10X, min_bg_count=25, max_bg_count=150, n_pcs=5, expected_targets=10000, min_count=200, p=0.5, seedn=NULL){
	x <- read_10x(path_10X)
	sce <- create_SCE(x)
	sce <- subset_dropls(sce, min_c=min_bg_count)
	sce <- normalize(sce)
	sce <- get_var_genes(sce)
	sce <- get_pcs(sce, n_pcs)
	sce <- subset_n_dropls(sce, 2*expected_targets, slot_name=c("raw", "norm"))
	sce <- get_bg_cand(sce, min_bg_count, max_bg_count, expected_targets)
	print(head(sce@bg_info$Labels))
	print(dim(sce@norm))
	sce <- run_em_pcs(sce, n_pcs=n_pcs, seedn=seedn)
	sce <- call_targets(sce, min_count=min_count, p=p)
	return(sce)
}

run_diem_10x <- function(path_10X, min_bg_count=25, max_bg_count=150, n_pcs=5, expected_targets=10000, min_count=200, p=0.5, seedn=NULL){
	sce_list <- lapply(path_10X, FUN=DIEM_10X,
					   min_bg_count=min_bg_count, max_bg_count=max_bg_count, n_pcs=n_pcs, expected_targets=expected_targets, min_count=min_count, p=p, seedn=seedn)
	return(sce_list)
}
