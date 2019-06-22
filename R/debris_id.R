
#' Get PCs from expression data
#'
#' Runs irlba implementation of prcomp. Takes normalized expression data in 
#' \code{x@sim@norm}, subsets using genes in the \code{x@de@deg}.
#' The first \code{n_pcs} are calculated.
#' The return object for \code{prcomp_irlba} is placed in \code{x@pcs}. The 
#' principal component scores for the droplets are thus stored in \code{x@pcs$x}.
#'
#' @param x SCE. An SCE object.
#' @param n_pcs Integer. Number of PCs to calculate.
#' @param genes Vector. Use these genes specified.
#' @param verbose Boolean.
#'
#' @return An SCE object with \code{\link[irlba]{prcomp_irlba}} output in 
#' the \code{pcs} slot
#' @importFrom irlba prcomp_irlba
#' @importFrom Matrix colSums
#' @export
get_pcs <- function(x, 
					n_pcs=30, 
					genes=NULL, 
					verbose=FALSE){
	if (verbose) cat("Running PCA\n")
	if (is.null(genes)){
		gn <- x@de@deg
	} else {
		gn <- genes
	}
	if (n_pcs < 30) n_pcs <- 30
	if (length(gn) <= n_pcs) n_pcs <- length(gn)-1
	if (length(gn) == 0){
		stop(paste0("No genes listed in slot ", genes))
	}
	expr <- Matrix::t(x@sim@norm[gn,])
	expr <- expr[,Matrix::colSums(expr) > 0] # remove genes with no counts
	if (length(gn) > 100){
		pcs <- irlba::prcomp_irlba(expr, n=n_pcs, retx=TRUE, center=TRUE, scale.=TRUE)
	} else {
		pcs <- prcomp(expr, retx=TRUE, center=TRUE, scale.=TRUE)
	}
	rownames(pcs$x) <- rownames(expr)
	rownames(pcs$rotation) <- colnames(expr)
	colnames(pcs$x) <- paste0("PC", as.character(1:ncol(pcs$x)))
	x@sim@pcs <- pcs
	if (verbose) cat("Found PCs\n")
	return(x)
}

#' @export
sweep_cols <- function(x, vec){
	d <- Matrix::Diagonal(x = 1/vec)
	x <- x %*% d
}

#' Get proportion of counts in \code{rows}in x.
#'
#' @param x Matrix. Matrix object.
#' @param rows Character. rownames in x to calculate proportions of.
#'
#' @return Numeric vector.
#' @importFrom Matrix Diagonal colSums
#' @export
get_col_props <- function(x, rows){
	x <- x  %*% Matrix::Diagonal(x = 1/Matrix::colSums(x))
	col_props <- Matrix::colSums(x[rows,])
	names(col_props) <- colnames(x)
	return(col_props)
}

#' Get background score for each cell
#'
#' Background score is calculated as the percentage of reads mapping to a 
#' DE gene enriched in the background droplets. Differentially expressed genes 
#' must be calculated beforehand.
#'
#' @param x SCE. SCE object.
#'
#' @return SCE object.
#' @export
get_bgscore <- function(x){
	log2fc <- x@de@table[x@de@deg,"log2fc"]
	log2fc <- Matrix::Matrix(log2fc)
	bg_score <- Matrix::t(log2fc) %*% x@sim@norm[x@de@deg,]
	bg_score <- as.numeric(bg_score)
	names(bg_score) <- colnames(x@sim@norm)
    # bg_genes <- x@de@deg_low
    # bg_score <- get_col_props(x@sim@counts, bg_genes)
    x@sim@bg_score <- bg_score
    # Add bg score to dropl_info slot for candidate droplets
    bg_score <- bg_score[ names(bg_score) %in% rownames(x@dropl_info) ]
    x@dropl_info[names(bg_score),"bg_score"] <- bg_score
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
#' Selects only PCs that are correlated with the background score to avoid clustering based on biological cell type 
#' variation. PCs with a pearson correlation coefficient greater than or equal to \code{min_r} are used.
#' EM is run \code{n_runs} times with random initializations, and parameters with 
#' the best log likelihood are output. Within EM, iterate a minimum of \code{min_iter} times and 
#' a maximum of \code{max_iter} times. The \code{eps} value gives the threshold of percentage change in 
#' log likelihood until convergence is reached. \code{seedn} gives the seed for random number generation, 
#' which can be used to reproduce results.
#'
#' @param x SCE. SCE object with PCs and labels.
#' @param min_r Numeric. Only run PCs correlated to the background score. PCs are correlated with the background 
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
run_em_pcs <- function(x, 
					   min_r=0.5, 
					   min_iter=5, 
					   max_iter=1000, 
					   eps=1e-10, 
					   n_runs=10, 
					   seedn=NULL, 
					   verbose=TRUE){
	runs <- list()

	# Select PCs that are correlated with bg score
	cors <- cor(x@sim@pcs$x, x@sim@bg_score)
	cors <- abs(cors)
	cord_pcs <- (1:ncol(x@sim@pcs$x))[which(cors >= min_r)]
	if (length(cord_pcs) == 0){
		stop("No PCs are correlated with background score. Cannot run EM on uncorrelated PCs.")
	}
	df <- as.data.frame(cbind(x@sim@bg_score, x@sim@pcs$x[,cord_pcs]))
	if (verbose){
		cat(paste0("Using PCs ", paste0(as.character(cord_pcs), collapse=","), " for EM\n"))
	}
	df <- as.data.frame(x@sim@pcs$x[,cord_pcs])

	# Run EM
	if (verbose){
		cat("Running semi-supervised EM\n")
	}
	for (i in 1:n_runs){
		runs[[i]] <- run_mv_em_diag(df, k=2, min_iter=min_iter, max_iter=max_iter,
									labels = x@sim@labels, eps=eps, seedn=seedn, verbose=FALSE)
		if (!is.null(seedn)) seedn <- seedn + 1
	}
	if (verbose){
		cat("Finished EM\n")
	}

	# Get run with max likelihood of parameters
	llks <- sapply(runs, function(x) x@llks[length(x@llks)])
	if (sum(!is.na(llks)) == 0){
		stop("No EM runs converged successfully.\n")
	}
	max_i <- which.max(llks)
	x@sim@emo <- runs[[max_i]]
	colnames(x@sim@emo@Z) <- c("Background", "Target")
	x@sim@emo@pcs <- cord_pcs
	return(x)
}

##' Run EM on PCs + BgScore
##' 
##' Run expectation maximization (EM) using diagonal covariance on expression PCs in an SCE object with k=2 groups. 
##' EM is run \code{n_runs} times with random initializations, and the EM output parameters with 
##' the best log likelihood is output. Within EM, iterate a minimum of \code{min_iter} times and 
##' a maximum of \code{max_iter} times. The \code{eps} value gives the threshold of percentage change in 
##' log likelihood until covnergence is called. \code{seedn} gives the seed for random number generation, 
##' which can be used to reproduce results.
##'
##' @param x SCE. SCE object with PCs and labels
##' @param n_pcs Integer. Number of expression PCs to use for running EM
##' @param min_iter Integer. Minimum number of iterations
##' @param max_iter Integer. Maximum number of iterations
##' @param eps Numeric. Epsilon of log-likilhood change to call convergence
##' @param n_runs Integer. Number of EM runs, with best result returned
##' @param seedn Numeric. Seed for random number generation
##' @param verbose Boolean. Print out logging information
##'
##' @return SCE object with EM output in the \code{emo} slot. See \code{\link{run_mv_em_diag}} for details
##' @export
#run_em_bg_score_pcs <- function(x, 
#								select_pcs=TRUE, 
#					   n_pcs=3, 
#					   min_iter=5, 
#					   max_iter=1000, 
#					   eps=1e-10, 
#					   n_runs=10, 
#					   seedn=NULL, 
#					   verbose=TRUE){
#	runs <- list()
#	if (is.null(n_pcs)) {n_pcs <- ncol(x@sim@pcs$x)}
#	if (n_pcs > 0){
#		df <- cbind(x@sim@bg_score, x@sim@pcs$x[,1:n_pcs])
#	} else {
#		df <- data.frame(bg_score=x@sim@bg_score)
#	}
#	if (select_pcs){
#		cord_pcs <- c()
#		pv_thresh <- 0.05 / ncol(x@sim@pcs$x)
#		for (i in 1:ncol(x@sim@pcs$x)){
#			ct <- cor.test(x@sim@bg_score, x@sim@pcs$x[,i], use="pairwise", method="pearson")
#			if (ct$p.value < pv_thresh){
#				cord_pcs <- c(cord_pcs, i)
#			}
#		}
#		if (length(cord_pcs) == 0){
#			stop("No PCs are correlated with background score. Cannot run EM on uncorrelated PCs.")
#		}
#		df <- cbind(x@sim@bg_score, x@sim@pcs$x[,cord_pcs])
#	}
#	if (verbose){
#		cat("Running semi-supervised EM\n")
#	}
#	for (i in 1:n_runs){
#		runs[[i]] <- run_mv_em_diag(df, k=2, min_iter=min_iter, max_iter=max_iter,
#									labels = x@sim@labels, eps=eps, seedn=seedn, verbose=FALSE)
#		if (!is.null(seedn)) seedn <- seedn + 1
#	}
#	if (verbose){
#		cat("Finished EM\n")
#	}
#	llks <- sapply(runs, function(x) x@llks[length(x@llks)])
#	if (sum(!is.na(llks)) == 0){
#		stop("No EM runs converged successfully.\n")
#	}
#	max_i <- which.max(llks)
#	x@sim@emo <- runs[[max_i]]
#	colnames(x@sim@emo@Z) <- c("Background", "Target")
#	x@sim@emo@pcs <- cord_pcs
#	return(x)
#}

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
	keep_names <- rownames(x@sim@emo@Z)[ x@sim@emo@Z[,2] > lk_fraction ]
	x@dropl_info[,"Target"] <- rownames(x@dropl_info) %in% keep_names
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
