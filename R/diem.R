
#' Run diem pipeline
#'
#' Run DIEM pipeline on an SCE object, typically returned from \code{\link{read_10x}} function.
#'
#' @param sce SCE. SCE object with raw counts.
#' @param min_counts Integer. Only droplets with at least this number of total counts are used for testing.
#' @param top_n Numeric. Only the top \code{top_n} droplets ranked by total counts are used for testing.
#' @param de_cutpoint Numeric. Manually set the DE cutpoint.
#' @param log2fc_thresh Numeric. Test only genes that have this log2 fold-change threshold between low and high
#' droplets will be considered as differentially expressed. Fold change is for both comparison low vs high and high vs low
#' @param de_p_thresh Numeric. Adjusted p-value threshold for differential expression.
#' @param de_correct Character. Correction method.
#' @param simf Numeric. Generate number of candidates * \code{simf} random samples.
#' @param gma Numeric. Multiply p by this number to generate alpha for Dirichlet-multinomial random sampling.
#' @param scale_factor Numeric. A scaling factor to multiply values after division by total count.
#' @param n_pcs Integer. Number of PCs to calculate.
#' @param genes Vector. Use these genes specified for PCA.
#' @param min_r Numeric. Only run PCs correlated to the background score. PCs are correlated with the background 
#' score and the absolute value of the R coefficient must be greater than this threshold.
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
				 min_counts=NULL, 
				 top_n=NULL, 
				 de_cutpoint=NULL, 
				 log2fc_thresh=0.5, 
				 de_p_thresh=0.05, 
				 de_correct="bonferroni", 
				 simf=0.5, 
				 gma=10, 
				 scale_factor=1, 
				 n_pcs=30, 
				 genes=NULL, 
				 min_r=0.5, 
				 min_iter=5, 
				 max_iter=1000, 
				 eps=1e-10, 
				 n_runs=10, 
				 lk_fraction=0.95, 
				 seedn=NULL, 
				 verbose=TRUE){
	sce <- set_min_counts(sce, min_counts=min_counts, top_n=top_n)
	sce <- set_de_cutpoint(sce, de_cutpoint=de_cutpoint, verbose=verbose)
	sce <- get_de(sce, log2fc_thresh=log2fc_thresh, de_p_thresh=de_p_thresh, de_correct=de_correct, verbose=verbose)
	sce <- sim_bg(sce, simf=simf, gma=gma, seedn=seedn, verbose=verbose)
	sce <- normalize(sce, scale_factor=scale_factor)
	sce <- get_bgscore(sce)
	sce <- get_pcs(sce, n_pcs=n_pcs, genes=genes, verbose=verbose)
	sce <- run_em_pcs(sce, min_r=min_r, min_iter=min_iter, max_iter=max_iter, eps=eps, n_runs=n_runs, seedn=seedn, verbose=verbose)
	sce <- call_targets(sce, lk_fraction=lk_fraction)
	return(sce)
}

