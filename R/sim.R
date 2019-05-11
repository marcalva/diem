
#' Bin by expression counts
#' 
#' Sum read counts for genes across droplets in SCE object. 
#' Given start, end, and by parameters, sums reads across 
#' genes in each bin. There are (end - start)/by bins. If 
#' this is a fraction, the number of bins are rounded down, 
#' with the same width for each bin.
#'
#' @param x SCE object
#' @param start_bin Integer. Read count for starting the bins
#' @param end_bin Integer. Read count for ending the bins
#' @param by_bin Numeric. Step size, the width of each bin
#'
#' @param Matrix with read counts summed across the rows
#' @export
bin_by_counts <- function(x, start_bin=0, end_bin=150, by_bin=150){
	n_bin <- floor((end_bin - start_bin)/by_bin)
	counts <- matrix(nrow=nrow(x@raw), ncol=n_bin)
	rownames(counts) <- rownames(x@raw)
	bins <- lapply(1:n_bin, function(x) c(start_bin + (x*by_bin) - by_bin, start_bin + (x*by_bin)))

	keep <- sapply(bins, function(bin) any(x@dropl_info[,"total_counts"] > bin[1] & x@dropl_info[,"total_counts"] <= bin[2]) )
	bins <- bins[keep]

	colnames(counts) <- sapply(bins, function(x) paste0("Bin", as.character(x[1]), "_", as.character(x[2])))

	for (i in 1:n_bin){
		counts[,i] <- Matrix::rowSums( x@raw[, x@dropl_info[,"total_counts"] > bins[[i]][1] & x@dropl_info[,"total_counts"] <= bins[[i]][2], drop=FALSE] )
	
		# Total sum to greatest counts in bin
		counts[,i] <- counts[,i]/sum(counts[,i])
		# counts[,i] <- bins[[i]][2] * counts[,i]
	}
	return(counts)
}

#' Create matrix with multinomial draws
#'
#' Draws random samples from a multinomial distribution. Takes \code{n_draw} draws 
#' from a multinomial with parameter p (probability) and k (size). The parameter p 
#' is stored in the columns of matrix \code{bins}. The column to select for each 
#' draw is given in the vector \code{bin_sampled}, which should be the same size as 
#' \code{n_draw}. The sizes of the draws k are given in the vector \code{size_sampled}, 
#' which also needs to be of length \code{n_draws}.
#'
#' @param bins Matrix. Contains probabilities for multinomial draws
#' @param bin_sampled Vector. Column indices of bins to use for each draw
#' @param size_sampled Vector. Sizes to use for each draw
#' @param n_draw Integer. Number of random samples to take
#' @param seedn Boolean. Seed number
#'
#' @return Matrix with each column containing a random sample from the multinomial.
#' @export
draw_multinom <- function(bins, bin_sampled, size_sampled, n_draw, seedn=NULL){
	set.seed(seedn)
	sim <- matrix(nrow=nrow(bins), ncol=n_draw)
	for (i in 1:n_draw){
		sim[,i] <- rmultinom(n=1, size=size_sampled[i], prob=bins[,bin_sampled[i]])
	}
	rownames(sim) <- rownames(bins)
	colnames(sim) <- paste0("SIM", as.character(1:n_draw))
	return(sim)
}

#' Simulate background droplets with higher read counts
#'
#' Given an SCE object with raw expression data containing droplets from background 
#' and targets (cells/nuclei), generate simulated background droplets with the same 
#' count distribution as candidate targets. Then replace the expression data with 
#' counts only from the simulated droplets and the candidate target droplets. Labels 
#' the simulated background droplets in the column "background" of \code{x@dropl_info}. 
#' Candidate targets are the droplets with total counts within the \code{count_range} 
#' range and total genes detected in the \code{genes_detected} range. The simulated 
#' background droplets are generated from a multinomial distribution, which is estimated 
#' from the background droplets. The background droplets are identified as those with
#' total reads counts in the range between \code{start_bin} and \code{end_bin}. By default, 
#' only 1 multinomial distribution is estimated. However, several distributions can 
#' be estimated from bins of \code{by_bin} width between \code{start_bin} and \code{end_bin}.
#' 
#' @param x SCE.
#' @param count_range Vector. Range of read counts to consider as candidate targets
#' @param gene_range Vector. Range of genes detected to consider as candidate targets
#' @param simf Numeric. Generate number of candidates * simf random samples
#' @param start_bin Integer. Start range for read counts to generate multinomial distribution(s)
#' @param end_bin Integer. End range for read counts to generate multinomial distribution(s)
#' @param by_bin Numeric. Bin size to use for generating multinomial distribution(s)
#' @param seedn Numeric. Number to set seed
#' @param verbose Boolean
#'
#' @return SCE object with raw expression data replaced by simulated background droplets and 
#' candidate targets. Simulated background droplets are labeled in the background column 
#' of \code{x@dropl_info} as a 1, while candidate targets are labeled as 0.
#' @export
set_expression <- function(x, 
						   count_range=c(199,Inf), 
						   gene_range=c(0, Inf), 
						   simf=2, 
						   start_bin=0, 
						   end_bin=150, 
						   by_bin=NULL, 
						   seedn=NULL,
						   verbose=FALSE){
	if (verbose) cat("Simulating background droplets\n")

	# Use 1 bin by default.

	if (is.null(by_bin)) by_bin <- (end_bin - start_bin)
	# Select candidate targets from sce object x
	candidate_bool <- ( x@dropl_info[,"total_counts"] > count_range[1] ) &
	( x@dropl_info[,"total_counts"] <= count_range[2] ) &
	( x@dropl_info[,"n_genes"] > gene_range[1] ) &
	( x@dropl_info[,"n_genes"] <= gene_range[2] )

	if ( sum(candidate_bool) == 0 ){
		stop("No candidate droplets meet count_range and gene_range requirements")
	}

	candd <- x@raw[, candidate_bool]
	count_sizes <- x@dropl_info[candidate_bool,"total_counts"] # distribution of read counts for multinomial simulation
	
	# Multinomial distributions binned from count ranges. Used to simulate background droplets
	if (verbose) cat("Binning background droplets...\n")
	binned <- bin_by_counts(x, start_bin=start_bin, end_bin=end_bin, by_bin=by_bin)
	probs <- rep( 1/ncol(binned), times=ncol(binned))

	n_draw <- simf*ncol(candd)

	if (verbose) cat("Drawing from multinomial...\n")
	sizes_sampled <- sample(x=count_sizes, size=n_draw, replace=TRUE)
	bins_sampled <- sample(1:ncol(binned), size=n_draw, prob=probs, replace=TRUE)
	simM <- draw_multinom(bins=binned, bin_sampled=bins_sampled, size_sampled=sizes_sampled, n_draw, seedn=seedn)

	simM <- Matrix::Matrix(simM, sparse=TRUE)
	all_expr <- cbind(simM, candd)
	x@raw <- all_expr
	x <- fill_counts(x)

	# Label background
	labels <- rep(0, ncol(x@raw)); names(labels) = colnames(x@raw)
	labels[1:n_draw] <- 1
	x@dropl_info[,"background"] <- labels

	if (verbose) cat("Simulated background droplets\n")
	return(x)
}
