
#' Run diem pipeline
#'
#' Run DIEM pipeline on an SCE object, typically returned from \code{\link{read_10x}} function.
#'
#' @param sce SCE. SCE object with raw counts.
#' @param log_base_em Numeric. Take droplets with counts above 2 logs below the max droplet count using this base
#' for EM classification
#' @param top_n Numeric. Only run EM to classify droplets if they are in the 
#' top number of droplets ranked by total counts.
#' @param min_counts Integer. Only run EM to classify droplets with at least this number of total counts.
#' @param log_base_label Numeric. Fix the labels as nuclei for droplets with counts above 1 log below 
#' the max droplet count using this base.
#' @param de_cutpoint Numeric. Manually set the DE cutpoint. If not NULL, it is calculated from the data.
#' @param log2fc_thresh Numeric. Test only genes that have an absolute log2 fold-change above this threshold.
#' @param de_p_thresh Numeric. Adjusted p-value threshold.
#' @param de_correct Character. p-value correction method.
#' @param scale_factor Numeric. A scaling factor to multiply values after division by total count.
#' @param min_iter Integer. Minimum number of iterations in EM.
#' @param max_iter Integer. Maximum number of iterations in EM.
#' @param eps Numeric. Epsilon of log-likilhood change to call convergence in EM.
#' @param n_runs Integer. Number of EM initializations to run, with best result returned.
#' @param lk_fraction Numeric. Select droplets that have a fraction of target likelihood greater than this number.
#' @param seedn Numeric. Seed for random number generation.
#' @param verbose Boolean. Print out logging information.
#' 
#' @return SCE object.
#' @useDynLib diem
#' @export
diem <- function(sce, 
				 log_base_em=10, 
				 top_n=NULL, 
				 min_counts=NULL, 
				 log_base_label=5, 
				 de_cutpoint=NULL, 
				 log2fc_thresh=0, 
				 de_p_thresh=0.05, 
				 de_correct="fdr", 
				 scale_factor=1e4, 
				 min_iter=5, 
				 max_iter=1000, 
				 eps=1e-10, 
				 n_runs=10, 
				 lk_fraction=0.95, 
				 seedn=NULL, 
				 verbose=TRUE){
	sce <- subset_for_em(sce, log_base=log_base_em, top_n=top_n, min_counts=min_counts)
	sce <- set_labels(sce, log_base=log_base_label)
	sce <- get_de_genes(sce, de_cutpoint=de_cutpoint, log2fc_thresh=log2fc_thresh, de_p_thresh=de_p_thresh, de_correct=de_correct, verbose=verbose)
	sce <- normalize(sce, scale_factor=scale_factor)
	sce <- get_pi(sce)
	sce <- run_em(sce, min_iter=min_iter, max_iter=max_iter, eps=eps, n_runs=n_runs, seedn=seedn, verbose=verbose)
	sce <- call_targets(sce, lk_fraction=lk_fraction)
	return(sce)
}

