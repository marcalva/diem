
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
#' @param min_count Numeric. Minimum number of counts for a droplet to be considered a cell/nucleus
#' @param min_gene Numeric. Minimum number of genes for a droplet to be considered a cell/nucleus
#' @param simf Numeric. Factor to multiply the number of candidates to get the number of background droplets simulated
#' @param p Numeric. Minimum membership probability to call a target
#' @param seedn Numeric. Seed for random number generation
#' @param verbose Boolean. Print out logging information
#' 
#' @return SCE object with droplet IDs labeled as background or nucleus.
#' @useDynLib diem
#' @export
diem <- function(sce, 
				 n_pcs=5, 
				 n_runs=10, 
				 min_bg_count=0, 
				 max_bg_count=150, 
				 min_count=199, 
				 min_gene=199, 
				 simf=2, 
				 p=0.95, 
				 seedn=NULL, 
				 verbose=TRUE){
	sce <- get_de_genes(sce, low_count=c(min_bg_count,max_bg_count), high_count=c(min_count, Inf), verbose=verbose)
	sce <- set_expression(sce, count_range=c(min_count, Inf), gene_range=c(min_gene, Inf), simf=simf, seedn=seedn, verbose=verbose)
	sce <- normalize(sce, verbose=verbose)
	sce <- get_bgscore(sce)
	sce <- get_pcs(sce, n_pcs, genes="deg", verbose=verbose)
	sce <- run_em_pcs(sce, n_pcs=n_pcs, n_runs=n_runs, seedn=seedn, verbose=verbose)
	sce <- call_targets(sce, min_gene=min_gene, p=p)
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

