
#' Evaluate results from iteration
#'
#' @param x SCE object.
#' @param thresh Numeric. The percent change (a numeric from 0 to 1).
#' @param verbose Boolean. Verbosity.
#'
#' @return A boolean indicating whether the change in nuclei droplets is under the threshold.
eval_iter <- function(x, thresh=0.05, verbose=FALSE){
	if (length(x@diem@prev_nuclei) == 0){
		return(FALSE)
	}
	nuc1 <- sort(x@diem@prev_nuclei)
	nuc2 <- sort(targets_ids(x))
	# What is in the new set that isn't in the old set
	s_diff1 <- setdiff(nuc2, nuc1)
	# What is in the old set that isn't in the new set
	s_diff2 <- setdiff(nuc1, nuc2)
	# Combine
	s_diff <- c(s_diff1, s_diff2)
	if (verbose){
		cat(paste0(as.character(length(s_diff1)), " added, ", as.character(length(s_diff2)), " removed.\n"))
	}
	pct_change <- length(s_diff)/length(nuc1)
	# If number in new is less than pct thresh, return TRUE
	if (pct_change < thresh) {
		return(TRUE)
	} else {
		return(FALSE)
	}
}

diem_iter <- function(sce, 
					  cpm_thresh=3, 
					  log2fc_thresh=0.25, 
					  pct_tails=0.05, 
					  min_iter_em=5, 
					  max_iter_em=1000, 
					  eps_em=1e-10, 
					  lk_fraction=0.95, 
					  verbose=TRUE){
	sce <- get_log2fc(sce, cpm_thresh=cpm_thresh, log2fc_thresh=log2fc_thresh)
	sce <- get_pi(sce)
	sce <- run_em(sce, 
				  pct_tails=pct_tails, 
				  min_iter=min_iter_em, 
				  max_iter=max_iter_em, 
				  eps=eps_em, 
				  verbose=verbose)
	sce <- call_targets(sce, lk_fraction=lk_fraction)
	sce <- set_de_groups(sce)
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
#' @param top_n Numeric. Only run EM classification on the top top_n droplets ranked by total counts.
#' @param logt Boolean. Log1p transform counts during normalization.
#' @param cpm_thresh Numeric. Calculate pi from genes with counts per million mapped reads (CPM) above this threshold.
#' @param log2fc_thresh Numeric. Calculate pi from genes with log2 fold changes above this threshold.
#' @param pct_tails Numeric. Value from 0 to 0.5. Initialize means of low and high groups using bottom and top percentiles, 
#' respectively
#' @param min_iter_em Numeric. Minimum number of iterations during EM.
#' @param max_iter_em Numeric. Maximum number of iterations during EM.
#' @param eps_em Numeric. Epsilon of log-likilhood change to call convergence in EM.
#' @param lk_fraction Numeric. Select droplets that have a fraction of target likelihood greater than this number.
#' @param verbose Boolean. Verbosity.
#' 
#' @return SCE object.
#' @export
diem_init <- function(sce, 
				 top_n=1e4, 
				 init_de_log_base=50, 
				 init_de_cutpoint=NULL, 
				 logt=TRUE, 
				 cpm_thresh=3, 
				 log2fc_thresh=0.25, 
				 pct_tails=0.05, 
				 min_iter_em=5, 
				 max_iter_em=1000, 
				 eps_em=1e-10, 
				 lk_fraction=0.99, 
				 verbose=TRUE){
	if (verbose) cat("Initializating with iteration 1\n")
	
	# Initializations
	sce <- set_test_set(sce, top_n=top_n)
	sce <- init_de_cutpoint(sce, log_base=init_de_log_base, de_cutpoint=init_de_cutpoint, verbose=verbose)
	sce <- normalize(sce, logt=logt, verbose=verbose)
	sce <- fix_tails(sce, pct=pct_tails)

	# Iteration run
	sce <- diem_iter(sce, 
					 cpm_thresh=cpm_thresh, 
					 log2fc_thresh=log2fc_thresh, 
					 pct_tails=pct_tails, 
					 min_iter_em=min_iter_em, 
					 max_iter_em=max_iter_em, 
					 eps_em=eps_em, 
					 lk_fraction=lk_fraction, 
					 verbose=verbose)


	# Save initialization as first run
	sce@prev_iter[[1]] <- sce@dropl_info[!is.na(sce@dropl_info$pi_l), c("pi_l", "pi_h", "Call")]
	n_nuc <- length(sce@high_droplets)

	if (verbose){
		cat(paste0("Number of nuclei identified: ", as.character(n_nuc), "\n"))
	}

	if (length(sce@high_droplets) == 0){
		cat("Warning: No nuclei found in initialization\n")
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
#' @param top_n Numeric. Only run EM classification on the top top_n droplets ranked by total counts.
#' @param logt Boolean. Log1p transform counts during normalization.
#' @param cpm_thresh Numeric. Calculate pi from genes with counts per million mapped reads (CPM) above this threshold.
#' @param log2fc_thresh Numeric. Calculate pi from genes with log2 fold changes above this threshold.
#' @param pct_tails Numeric. Value from 0 to 0.5. Initialize means of low and high groups using bottom and top percentiles, 
#' respectively
#' @param min_iter_em Numeric. Minimum number of iterations during EM.
#' @param max_iter_em Numeric. Maximum number of iterations during EM.
#' @param eps_em Numeric. Epsilon of log-likilhood change to call convergence in EM.
#' @param lk_fraction Numeric. Select droplets that have a fraction of target likelihood greater than this number.
#' @param max_iter_diem Numeric. Maximum number of iterations during DIEM until convergence of labeling.
#' @param eps_iter_diem Numeric. Percent threshold for change in labeling until convergence.
#' @param verbose Boolean. Verbosity.
#' 
#' @return SCE object.
#' @export
diem <- function(sce,
				 top_n=1e4,  
				 init_de_log_base=50, 
				 init_de_cutpoint=NULL, 
				 logt=TRUE, 
				 cpm_thresh=3,  
				 log2fc_thresh=0.25, 
				 pct_tails=0.05, 
				 min_iter_em=5, 
				 max_iter_em=1000, 
				 eps_em=1e-10, 
				 lk_fraction=0.99, 
				 max_iter_diem=25, 
				 eps_iter_diem=0.05, 
				 verbose=TRUE){
	# Initialize if not already
	if (length(sce@prev_iter) == 0){
		sce <- diem_init(sce, top_n=top_n, init_de_log_base=init_de_log_base,
			init_de_cutpoint=init_de_cutpoint, logt=logt, 
			cpm_thresh=cpm_thresh, log2fc_thresh=log2fc_thresh, pct_tails=pct_tails, 
			min_iter_em=min_iter_em, max_iter_em=max_iter_em, eps_em=eps_em,
			lk_fraction=lk_fraction, verbose=verbose)
	}

	if (length(sce@high_droplets) == 0){
		cat("Warning: Converged to 0 nuclei. Check iterations\n")
		return(sce)
	}

	# Run through DIEM iterations
	ix <- 1
	converged <- FALSE
	while (ix < max_iter_diem){
		ix <- ix + 1
		if (verbose) cat(paste0("Iteration ", as.character(ix), "\n"))

		# Run DIEM algorithm
		sce <- diem_iter(sce, 
						 cpm_thresh=cpm_thresh, 
						 log2fc_thresh=log2fc_thresh, 
						 pct_tails=pct_tails, 
						 min_iter_em=min_iter_em, 
						 max_iter_em=max_iter_em, 
						 eps_em=eps_em, 
						 lk_fraction=lk_fraction, 
						 verbose=verbose)

		# If no nuclei found
		if (length(sce@high_droplets) == 0){
			cat("Warning: Converged to 0 nuclei. Check iterations\n")
			return(sce)
		}

		# Check if converged
		if (eval_iter(sce, thresh=eps_iter_diem, verbose=verbose)){
			sce@prev_iter[[ix]] <- sce@dropl_info[!is.na(sce@dropl_info$pi_l),c("pi_l", "pi_h", "Call")]
			sce@diem@prev_nuclei <- sce@high_droplets
			n_nuc <- as.character(length(sce@high_droplets))
			if (verbose) cat(paste0("Converged!\nNumber of nuclei identified: ", n_nuc, ".\n"))
			converged <- TRUE
			break
		} else {
			sce@prev_iter[[ix]] <- sce@dropl_info[!is.na(sce@dropl_info$pi_l),c("pi_l", "pi_h", "Call")]
			sce@diem@prev_nuclei <- sce@high_droplets
			if (verbose) cat(paste0("Number of nuclei identified: ", as.character(length(sce@high_droplets)), "\n"))
		}
	}
	if (!converged){
		cat(paste0("Warning: DIEM did not converge within ", as.character(max_iter_diem), " iterations.\n"))
	}
	return(sce)
}
