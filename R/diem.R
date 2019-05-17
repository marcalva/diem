
#' Run diem pipeline
#'
#' Run DIEM on an SCE object, typically returned from \code{\link{read_10x}} function.
#' DIEM first simulates background droplets by randomly sampling from a multinomial 
#' distribution. The parameters (gene probabilities) from the distribution are derived 
#' from droplets with low read counts, as the majority of these are likely to represent 
#' background RNA contamination. The number of simulated background droplets is given by 
#' the parameter \code{simf}, which multiplies the number of candidate droplets tested. 
#' PCA is run on the expression matrix containing the simulated and candidate droplets. 
#' DIEM then runs semi-supervised EM on the first \code{n_pcs} PCs. In order to call 
#' droplets as either background or target, the likelihood fraction must be greater than 
#' \code{p} (0.95 by default). This also calculates the percentage of reads that originate 
#' from mitochondrial genes, or from MALAT1.
#'
#' @param sce SCE. ACE object with raw counts
#' @param n_pcs Integer. Number of expression PCs to use for running EM
#' @param n_runs Integer. Number of random starts for EM run, with best result returned
#' @param min_bg_count Numeric. Minimum number of counts for a droplet to be considered background
#' @param max_bg_count Numeric. Maximum number of counts for a droplet to be considered background
#' @param min_tg_count Numeric. Minimum number of counts for a droplet to be considered target
#' @param max_tg_count Numeric. Maximum number of counts for a droplet to be considered target
#' @param min_tg_gene Numeric. Minimum number of genes detected for a droplet to be considered target
#' @param max_tg_gene Numeric. Maximum number of genes detected for a droplet to be considered target
#' @param top_n_tg Numeric. Set the top n ranked by counts as the threshold for considering targets
#' @param top_n_bg Numeric. Set the top n ranked by counts as the threshold for considering background
#' @param n_deg Integer. Number of differentially expressed genes between low/high count droplets to use for calculating PCA
#' @param simf Numeric. Factor to multiply the number of candidates to get the number of background droplets simulated
#' @param n_bin Numeric. Number of bins to sample from between min_bg_count and max_bg_count
#' @param shape_mult Numeric. Number to multiply gene probabilities for sampling from dirichlet
#' @param p Numeric. Minimum membership probability to call a target
#' @param seedn Numeric. Seed for random number generation
#' @param verbose Boolean. Print out logging information
#' 
#' @return SCE object with droplet IDs labeled as background or nucleus.
#' @useDynLib diem
#' @export
diem <- function(sce, 
				 n_pcs=1, 
				 n_runs=10, 
				 min_bg_count=0, 
				 max_bg_count=150, 
				 min_tg_count=199, 
				 max_tg_count=Inf,
				 min_tg_gene=199,
				 max_tg_gene=Inf,
				 top_n_tg=NULL,
				 top_n_bg=NULL, 
				 n_deg=2000, 
				 simf=1, 
				 n_bin=1,
				 shape_mult=100,
				 p=0.5, 
				 seedn=NULL, 
				 verbose=TRUE){
	sce <- set_limits(sce, min_bg_count=min_bg_count, max_bg_count=max_bg_count, 
					  min_tg_count=min_tg_count, max_tg_count=max_tg_count, 
					  min_tg_gene=min_tg_gene, max_tg_gene=max_tg_gene,
					  top_n_tg=top_n_tg, top_n_bg=top_n_bg)
	sce <- get_de_genes(sce, n_genes=n_deg, verbose=verbose)
	sce <- get_var_genes(sce, nf=n_deg)
	sce <- set_expression(sce, simf=simf, n_bin=1, seedn=seedn, shape_mult=shape_mult, verbose=verbose)
	sce <- normalize(sce, verbose=verbose)
	sce <- get_bgscore(sce)
	sce <- get_pcs(sce, genes="deg", verbose=verbose)
	sce <- run_em_bg_score_pcs(sce, n_pcs=n_pcs, n_runs=n_runs, seedn=seedn, verbose=verbose)
	sce <- call_targets(sce, p=p)
	sce <- get_mt_malat1(sce)
	return(sce)
}

#' Run diem pipeline for 10X samples
#'
#' Wrapper to run diem pipeline on 10X samples, returning \code{SCE} objects.
#'
#' @return List of SCE objects.
#' @export
run_diem_10x_sim <- function(path_10X, verbose=FALSE, ...){
	sce_list <- list()
	for (path in path_10X){
		if (verbose) cat(paste0("Running ", path, "\n"))
		sce <- read_10x(path)
		sce_list[[path]] <- diem(sce, verbose=verbose, ...)
	}
	return(sce_list)
}

