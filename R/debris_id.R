

#' Get PCs from expression data
#'
#' Runs irlba implementation of prcomp. Takes normalized expression data in 
#' \code{x@norm}, subsets using genes in the \code{deg} or \code{vg} slot, 
#' specified with parameter \code{genes}. The first \code{n_pcs} are calculated.
#' The return object for \code{prcomp_irlba} is placed in \code{x@pcs}. The 
#' principal component scores for the droplets are thus stored in \code{x@pcs$x}.
#'
#' @param x SCE. An SCE object
#' @param n_pcs Integer. Number of PCs to calculate
#' @param center Boolean. Whether to center variables
#' @param scale Boolean. Whether to scale variables to have unit variance
#' @param genes Character. Use these genes stored in SCE object. Either "vg" or "de"
#' @param verbose Boolean.
#'
#' @return An SCE object with \code{\link[irlba]{prcomp_irlba}} output in 
#' the \code{pcs} slot
#' @importFrom irlba prcomp_irlba
#' @export
get_pcs <- function(x, n_pcs=30, center = TRUE, scale. = TRUE, genes="deg", verbose=FALSE){
	if (verbose) cat("Running PCA\n")
	gn <- slot(x, genes)
	if (n_pcs < 30) n_pcs <- 30
	if (length(gn) <= n_pcs) n_pcs <- length(gn)-1
	if (length(gn) == 0){
		stop(paste0("No genes listed in slot ", genes))
	}
	if (length(gn) > 100){
	pcs <- prcomp_irlba(Matrix::t(x@norm[gn,]), 
						n=n_pcs, 
						retx=TRUE, 
						center=TRUE, 
						scale.=scale.)
	} else {
		pcs <- prcomp(Matrix::t(x@norm[gn,]), retx=TRUE, center=TRUE, scale.=scale.)
	}
	rownames(pcs$x) <- colnames(x@norm)
	rownames(pcs$rotation) <- gn
	colnames(pcs$x) <- paste0("PC", as.character(1:ncol(pcs$x)))
	x@pcs <- pcs
	if (verbose) cat("Found PCs\n")
	return(x)
}

#' @export
sweep_total <- function(x, vec){
	d <- Matrix::Diagonal(x = 1/vec)
	x <- x %*% d
}

#' Get background score for each cell
#'
#' Background score is calculated as the percentage of reads mapping to a 
#' DE gene enriched in the background droplets. Differentially expressed genes 
#' must be calculated beforehand. Background score is stored in \code{x@dropl_info} 
#' in the bg_score column.
#'
#' @param x SCE. SCE object.
#'
#' @return SCE object.
#' @export
get_bgscore <- function(x){
    expr <- x@raw
    expr <- sweep_total(expr, x@dropl_info[,"total_counts"])
    fc <- x@gene_info[x@deg, "fc", drop=FALSE]
    bg_genes <- rownames(fc)
    x@dropl_info[,"bg_score"] <- Matrix::colSums(expr[bg_genes,])
    return(x)
}

##' Specify background droplets and candidate droplets
##' 
##' This function returns the barcodes that are considered originating from background, and those that are
##' considered candidate targets (cells/nuclei). Barcodes are considered background if either they fall below
##' the top \code{top_n_cand} barcodes, or if the number of reads is less than bg_max.
##' 
##' total read/UMI count is greater than or equal to \code{bg_min} and less than or equal to \code{bg_max}. 
##' Candidate cells/nuclei are the top \code{top_n_cand} droplets ranked by read counts. This should be 
##' at least the expected number of targets recovered from the experiment. If the 
##' top \code{top_n_cand} droplets includes those with read/UMI counts < \code{bg_max}, the \code{bg_max} 
##' parameter is lowered so that \code{top_n_cand} candidates are always output.
##'
##' @param x SCE. An SCE object
##' @param bg_max Numeric. Maximum number of counts for a droplet to be considered background. Ignored if top_n_cand is set.
##' @param top_n_cand Numeric. Specifies the top number of droplets by expression to considered candidates
##'
##' @return SCE object with bg_info filled, a list with components
##' \item{Background}{Character vector of column names of x that are considered background droplets}
##' \item{Candidate}{Character vector of column names of x that are considered candidate droplets}
##' \item{Labels}{Numeric vector for droplets specifying 0 for candidate and 1 for background}
##' @export
#get_bg_cand <- function(x, bg_max=NULL, top_n_cand=NULL){
#	drplt_ids <- rownames(x@dropl_info)
#	labels <- rep(1, length(drplt_ids))
#	names(labels) <- drplt_ids
#
#	if (!is.null(top_n_cand)){
#		top_n_ix <- order(x@dropl_info[,"total_counts"], na.last = TRUE, decreasing = TRUE)[1:top_n_cand]
#		labels[top_n_ix] <- 0
#	} else if (!is.null(bg_max)){
#		labels[x@dropl_info[,"total_counts"] > bg_max] <- 0
#	}
#	else {
#		stop("Either top_n_cand or bg_max must be set")
#	}
#	x@dropl_info[,"background"] <- labels
#	return(x)
#}

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
					   n_pcs=1, 
					   min_iter=5, 
					   max_iter=1000, 
					   eps=1e-10, 
					   n_runs=10, 
					   seedn=NULL, 
					   verbose=TRUE){
	runs <- list()
	if (is.null(n_pcs)) {n_pcs <- ncol(x@pcs$x)}
	df <- x@pcs$x[,1:n_pcs,drop=FALSE]
	for (i in 1:n_runs){
		runs[[i]] <- run_mv_em_diag(df, k=2, min_iter=min_iter, max_iter=max_iter,
									labels = x@dropl_info[,"background"], eps=eps, seedn=seedn, verbose=verbose)
		if (!is.null(seedn)) seedn <- seedn + 1
	}
	for (i in 1:n_runs){
		if (!is.null(n_pcs)){
			runs[[i]] <- run_mv_em_diag(x@pcs$x[,1:n_pcs], k=2, min_iter=min_iter, max_iter=max_iter,
										labels = x@dropl_info[,"background"], eps=eps, seedn=seedn, verbose=verbose)
		} else {
			runs[[i]] <- run_mv_em_diag(x@pcs$x, k=2, min_iter=min_iter, max_iter=max_iter, 
										labels = x@dropl_info[,"background"], eps=eps, seedn=seedn, verbose=verbose)
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

#' Run EM on PCs + BgScore
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
run_em_bg_score_pcs <- function(x, 
					   n_pcs=3, 
					   min_iter=5, 
					   max_iter=1000, 
					   eps=1e-10, 
					   n_runs=10, 
					   seedn=NULL, 
					   verbose=TRUE){
	runs <- list()
	if (is.null(n_pcs)) {n_pcs <- ncol(x@pcs$x)}
	df <- cbind(x@dropl_info$bg_score, x@pcs$x[,1:n_pcs])
	for (i in 1:n_runs){
		runs[[i]] <- run_mv_em_diag(df, k=2, min_iter=min_iter, max_iter=max_iter,
									labels = x@dropl_info[,"background"], eps=eps, seedn=seedn, verbose=verbose)
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
#' Call targets (typically cells or nuclei) from droplets if the 
#' log-likelihood membership probability is higher than \code{p}.
#'
#' @param x SCE. SCE object
#' @param p Numeric. Proportion of likelihood to call a target
#'
#' @return SCE object
#' @export
call_targets <- function(x, p=0.95){
	keep <- (x@emo$Z[,2] > p) & (x@dropl_info[,"total_counts"] > x@limits$min_tg_count) & (x@dropl_info[,"n_genes"] > x@limits$min_tg_gene)
	x@dropl_info[,"Target"] <- keep
	return(x)
}

#' Return target IDs
#'
#' @param x SCE. SCE object
#'
#' @return Character vector with droplet IDs of called targets
#' @export
targets_ids <- function(x){
	return(rownames(x@dropl_info)[x@dropl_info[,"Target"]])
}

#' Get percent of reads aligning to MT genome and MALAT1 gene
#'
#' Places percent of reads aligning to MT genome and MALAT1 gene in 
#' \code{x@dropl_info}, under columns MALAT1 and MT_PCT. The gene 
#' names in the rows of \code{x@raw} should be MALAT1 and only mitochondrial 
#' genes should start with MT-.
#'
#' @param x SCE.
#'
#' @return SCE object
#' @export
get_mt_malat1 <- function(x){
	n_drop <- nrow(x@dropl_info)
	
	malat_gene <- grep("malat1", rownames(x@raw), ignore.case=T, value=TRUE)
	if ( length(malat_gene) == 0 ) malat1 <- NA
	else malat1 <- x@raw[malat_gene,] / x@dropl_info[,"total_counts"]

	mt_genes <- grep("MT-", rownames(x@raw), ignore.case=T, value=TRUE)
	mt_pct <- Matrix::colSums( x@raw[mt_genes,] ) / x@dropl_info[,"total_counts"]

	x@dropl_info[,"MALAT1"] <- malat1
	x@dropl_info[,"MT_PCT"] <- mt_pct
	return(x)
}
