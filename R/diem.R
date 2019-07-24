
#' Run diem pipeline
#'
#' Run DIEM pipeline on an SCE object, typically returned from \code{\link{read_10x}} function.
#'
#' @param sce SCE. SCE object with raw counts.
#' @param init_de_log_base Numeric. Take droplets with counts within 2 logs below the max droplet size using this log base
#' for initializing the DE cutpoint.
#' @param diff_thresh Numeric. Manually set the initialization of the DE cutpoint. If not NULL, it is calculated from the data.
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
				 min_counts=30, 
				 top_n_de=1e4, 
				 diff_thresh=NULL, 
				 cpm_thresh=3,  
				 log2fc_thresh=0, 
				 n_genes=100, 
				 init_em_tails=0.05, 
				 pp_thresh=0.95, 
				 max_iter=100, 
				 eps=1e-6, 
				 ss_pct_low=0.00, 
				 ss_pct_high=0.00, 
				 min_gene_call=200, 
				 verbose=TRUE){

	# sce@diem <- DIEM()
	if ((pp_thresh < 0) | (pp_thresh > 1)){
		stop("pp_thresh must be between 0 and 1")
	}
	sce@pp_thresh <- pp_thresh

	sce <- set_test_set(sce, top_n=top_n, min_counts=min_counts) # Set @test_IDs
	sce <- get_diff_thresh(sce, top_n=top_n_de, diff_thresh=diff_thresh, verbose=verbose) # Set @de_cut_init
	sce <- set_diff_genes(sce, n_genes=n_genes, cpm_thresh=cpm_thresh)
	sce <- fix_tails(sce, pct_low=ss_pct_low, pct_high=ss_pct_high) # Set @labels
	
	#===================================================================
	# Run EM
	#===================================================================
	sce <- run_em(x=sce, eps=eps, max_iter=max_iter, labels=sce@labels, verbose=verbose)

	#===================================================================
	#===================================================================
	sce <- call_targets(sce, pp_thresh=sce@pp_thresh, min_genes=min_gene_call)

	if (!sce@diem@converged){
		cat(paste0("Warning: DIEM did not converge within ", as.character(max_iter), " iterations.\n"))
	}

	sce <- fill_dropl_info(sce)

	return(sce)
}



