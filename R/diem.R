
#' Evaluate results from iteration for convergence
#'
#' @param x SCE object.
#' @param thresh Numeric. The percent change (a numeric from 0 to 1).
#' @param verbose Boolean. Verbosity.
#'
#' @return A boolean indicating whether the change in signal droplets is under the threshold.
eval_iter <- function(x, thresh=0.01, verbose=FALSE){
	iteration <- length(x@diem)
	if (iteration < 2){
		return(FALSE)
	}
	Zc <- x@diem[[iteration]]@emo[["Z"]]
	Zp <- x@diem[[iteration-1]]@emo[["Z"]]

	diffs <- abs(Zc - Zp)
	pct_change <- mean(diffs)
	if (pct_change < thresh) {
		return(TRUE)
	} else {
		return(FALSE)
	}
}

diem_iter <- function(sce, 
					  cpm_thresh=3, 
					  log2fc_thresh=0, 
					  init_em_tails=0.05, 
					  verbose=TRUE){
	iteration <- length(sce@diem)+1
	sce@diem[[iteration]] <- DIEM()

	sce <- get_log2fc(sce, cpm_thresh=cpm_thresh, log2fc_thresh=log2fc_thresh)
	sce <- get_pi(sce)
	sce <- run_em(sce, init_em_tails=init_em_tails, verbose=verbose)
	sce <- call_targets(sce)

	return(sce)
}

#' Initialize DIEM
#'
#' Initialize DIEM pipeline on an SCE object, e.g. one returned from \code{\link{read_10x}} function.
#'
#' @param sce SCE. SCE object with raw counts.
#' @param init_de_log_base Numeric. Take droplets with counts within 2 logs below the max droplet size using this log base
#' for initializing the DE cutpoint.
#' @param init_de_cutpoint Numeric. Manually set the initialization of the DE cutpoint. If not NULL, it is calculated from the data.
#' @param top_n Numeric. Only the top \code{top_n} droplets ranked by total counts are classified.
#' @param min_counts Integer. Only droplets with at least this number of total counts are classified. 
#'  Can override \code{top_n} if \code{top_n} includes droplets with less than \code{min_counts}.
#' @param cpm_thresh Numeric. Calculate pi from genes with counts per million mapped reads (CPM) above this threshold.
#' @param log2fc_thresh Numeric. Calculate pi from genes with log2 fold changes above this threshold.
#' @param lk_fraction Numeric. Select droplets that have a fraction of target likelihood greater than this number.
#' @param verbose Boolean. Verbosity.
#' 
#' @return SCE object.
#' @export
diem_init <- function(sce, 
				 top_n=1e4, 
				 min_counts=100, 
				 init_de_log_base=50, 
				 init_de_cutpoint=NULL, 
				 cpm_thresh=3, 
				 log2fc_thresh=0, 
				 init_em_tails=0.05, 
				 lk_fraction=0.95, 
				 ss_pct=0, 
				 verbose=TRUE){
	if (verbose) cat("Initializating with iteration 1\n")

	sce@diem <- list()
	sce@p_thresh <- lk_fraction
	sce@converged <- FALSE
	sce@iter_use <- numeric()

	# Initializations
	sce <- set_test_set(sce, top_n=top_n, min_counts=min_counts) # Set @test_IDs
	sce <- init_de_cutpoint(sce, log_base=init_de_log_base, de_cutpoint=init_de_cutpoint, verbose=verbose) # Set @de_cut_init
	sce <- normalize(sce, logt=TRUE, verbose=verbose) # Set @norm
	sce <- fix_tails(sce, pct=ss_pct) # Set @labels

	# First iteration
	sce <- diem_iter(sce, 
					 cpm_thresh=cpm_thresh, 
					 log2fc_thresh=log2fc_thresh, 
					 init_em_tails=init_em_tails, 
					 verbose=verbose)

	# Check separation of tails
	check_init_sep(sce)

	nsig <- length(targets_ids(sce))
	if (verbose){
		cat(paste0("Number of targets identified: ", as.character(nsig), "\n"))
	}

	return(sce)
}


#' Run diem pipeline
#'
#' Run DIEM pipeline on an SCE object, typically returned from \code{\link{read_10x}} function.
#'
#' @param sce SCE. SCE object with raw counts.
#' @param init_de_log_base Numeric. Take droplets with counts within 2 logs below the max droplet size using this log base
#' for initializing the DE cutpoint.
#' @param de_cutpoint Numeric. Manually set the initialization of the DE cutpoint. If not NULL, it is calculated from the data.
#' @param top_n Numeric. Only the top \code{top_n} droplets ranked by total counts are classified.
#' @param min_counts Integer. Only droplets with at least this number of total counts are classified. 
#'  Can override \code{top_n} if \code{top_n} includes droplets with less than \code{min_counts}.
#' @param cpm_thresh Numeric. Calculate pi from genes with counts per million mapped reads (CPM) above this threshold.
#' @param log2fc_thresh Numeric. Calculate pi from genes with log2 fold changes above this threshold.
#' @param lk_fraction Numeric. Select droplets that have a fraction of target likelihood greater than this number.
#' @param max_iter_diem Numeric. Maximum number of iterations during DIEM until convergence of labeling.
#' @param eps_iter_diem Numeric. Percent threshold for change in labeling until convergence.
#' @param verbose Boolean. Verbosity.
#' 
#' @return SCE object.
#' @export
diem <- function(sce,
				 top_n=1e4, 
				 min_counts=100, 
				 init_de_log_base=50, 
				 init_de_cutpoint=NULL, 
				 cpm_thresh=3,  
				 log2fc_thresh=0, 
				 init_em_tails=0.05, 
				 lk_fraction=0.95, 
				 max_iter_diem=25, 
				 eps_iter_diem=0.001, 
				 ss_pct=0, 
				 verbose=TRUE){

	if (length(sce@diem) == 0){
		sce <- diem_init(sce, top_n=top_n, ss_pct=ss_pct, 
						 min_counts=min_counts, init_de_log_base=init_de_log_base,
						 init_de_cutpoint=init_de_cutpoint, 
						 cpm_thresh=cpm_thresh, log2fc_thresh=log2fc_thresh, 
						 lk_fraction=lk_fraction, verbose=verbose)
	} else {
		sce@diem <- sce@diem[1]
	}

	if (length(targets_ids(sce)) == 0){
		cat("Warning: Converged to 0 targets. Check iterations\nReturning unfinished SCE object.\n")
		return(sce)
	}

	#===================================================================
	# Run through DIEM iterations
	#===================================================================
	ix <- length(sce@diem)
	while ( (ix < max_iter_diem) & (!sce@converged) ){
		ix <- ix + 1
		if (verbose) cat(paste0("Iteration ", as.character(ix), "\n"))

		# Run DIEM algorithm
		sce <- diem_iter(sce, 
						 cpm_thresh=cpm_thresh, 
						 log2fc_thresh=log2fc_thresh, 
						 init_em_tails=init_em_tails, 
						 verbose=verbose)

		nsig <- length(targets_ids(sce))

		# If nothing found
		if (nsig == 0){
			cat("Warning: Converged to 0 targets. Check iterations.\nReturning unfinished SCE object.\n")
			return(sce)
		}

		# Check if converged
		if (eval_iter(sce, thresh=eps_iter_diem, verbose=verbose)){
			nsig <- as.character(length(targets_ids(sce)))
			if (verbose) cat(paste0("Converged!\nNumber of targets identified: ", as.character(nsig), ".\n"))
			sce@converged <- TRUE
		} else {
			if (verbose) cat(paste0("Number of targets identified: ", as.character(length(targets_ids(sce))), "\n"))
		}
	}
	#===================================================================
	#===================================================================
	if (!sce@converged){
		ix <- 1
		cat(paste0("Warning: DIEM did not converge within ", as.character(max_iter_diem), " iterations.\n"))
		cat(paste0("Warning: Using calls from initialization for conservative filtering.\n"))
	}

	sce@iter_use <- ix
	sce <- fill_dropl_info(sce, iteration=ix)

	return(sce)
}
