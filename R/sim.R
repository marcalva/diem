
#' Bin by expression counts
#' 
#' Bins expression data across columns and sums values.
#' Given start, end, and by parameters, sums reads across 
#' genes in each bin. There are (end - start)/by bins. If 
#' this is a fraction, the number of bins are rounded down, 
#' with the same width for each bin. Start and end are calculated as 
#' greater than and less than or equal to, respectively.
#'
#' @param x Sparse matrix. Expression data in sparse matrix format
#' @param start_bin Integer. Read count for starting the bins
#' @param end_bin Integer. Read count for ending the bins
#' @param n_bin Numeric. Number of bins
#' @param bins List. If not null, use these bins. Each element is a bin with start and end.
#'
#' @param Matrix with read counts summed across the rows
#' @export
bin_by_counts <- function(x, start_bin=0, end_bin=NULL, n_bin=1, min_in_bin=1, bins=NULL){
	total_counts <- Matrix::colSums(x)
	if (is.null(bins)){
		if (is.null(end_bin)) end_bin <- max(total_counts)
		by_bin <- (end_bin - start_bin)/n_bin
		by_bin <- ceiling(by_bin)
		bins <- lapply(1:n_bin, function(x) c(start_bin + (x*by_bin) - by_bin, start_bin + (x*by_bin)))
		bins[[n_bin]][2] <- end_bin # Lengthen last bin if uneven

	}

	keep <- sapply(bins, function(bin) sum(total_counts > bin[1] & total_counts <= bin[2]) >= min_in_bin )
	bins <- bins[keep]
	n_bin <- length(bins)

	counts <- matrix(nrow=nrow(x), ncol=n_bin)
	rownames(counts) <- rownames(x)
	colnames(counts) <- sapply(bins, function(x) paste0("Bin", as.character(x[1]), "_", as.character(x[2])))

	for (i in 1:n_bin){
		counts[,i] <- Matrix::rowSums( x[, total_counts > bins[[i]][1] & total_counts <= bins[[i]][2], drop=FALSE] )
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

#' Estimate gamma
estimate_gamma <- function(x, ranges=NULL){
	p <- Matrix::colSums(x); p <- p / sum(p)
	if (is.null(ranges)){
		ranges <- seq(50,1000, 50)
	}
	loglik <- function(m, alpha){
		alpha_0 <- sum(alpha)
		logllks <- apply(m, 1, function(xi){
						 n <- sum(xi)
						 lfactorial(n) + 
						 lgamma(alpha_0) - 
						 lgamma(n + alpha_0) + 
						 sum(lgamma(xi + alpha)) - 
						 sum(lfactorial(xi)) - 
						 sum(lgamma(alpha))})
		return(sum(logllks))
	}
	llks <- sapply(ranges, function(a) {print(a); loglik(x, a*p)})
	return(ranges[which.max(llks)])
}

#' Create matrix with dirichlet-multinomial draws
#'
#' Draws random samples from a dirichlet-multinomial distribution. The \code{bins} 
#' matrix contains the probability values for sampling in the columns. The 
#' columns to use for each draw are given in the vector \code{bin_sampled}. This should 
#' be the same length as the value of \code{n_draw}. The dirichlet is sampled by 
#' taking a sample from the gamma distribution, using the probability in \code{bins} 
#' as the shape. The parameter \code{alpha} multiplies the probability vector in 
#' \code{bins}. The larger the value, the tighter the dirichlet sample, leading to 
#' the dirichlet-multinomial approaching the multinomial. After sampling the dirichlet, 
#' a multinomial of size taken from the vector \code{size_sampled} is taken.
#'
#' @param bins Matrix. Contains probabilities for random samples for dirichlet
#' @param bin_sampled Vector. Column indices of bins to use for each draw
#' @param size_sampled Vector. Sizes to use for each multinomial draw
#' @param n_draw Integer. Number of random samples to take
#' @param seedn Boolean. Seed number
#'
#' @return Matrix with each column containing a random sample from the dirichlet-multinomial.
#' @export
draw_diri_multinom <- function(dmp, seedn=NULL){
	set.seed(seedn)
	nv <- nrow(dmp@alphas)
	n_draw <- ncol(dmp@alphas)
	sim <- matrix(nrow=nv, ncol=n_draw)
	rownames(sim) <- rownames(dmp@alphas)
	colnames(sim) <- paste0("SIM", as.character(1:n_draw))
	for (i in 1:n_draw){
		sim[,i] <- rmultinom(n=1, size=dmp@sizes[i], prob=dmp@alphas[,i])
	}
	return(sim)
}

#' Simulate background droplets with higher read counts
#'
#' Given an SCE object with raw expression data containing droplets from background 
#' and targets (cells/nuclei), generate simulated background droplets with the same 
#' count distribution as candidate targets. Then replace the expression data with 
#' counts only from the simulated droplets and the candidate target droplets. Labels 
#' the simulated background droplets in the column "background" of \code{x@dropl_info}. 
#' Candidate targets are the droplets with total counts and genes detected as 
#' specified in the SCE object limits. The simulated 
#' background droplets are generated from a multinomial distribution, which is estimated 
#' from the background droplets. The background droplets are identified as those with
#' total reads counts in the range between \code{start_bin} and \code{end_bin}. By default, 
#' only 1 multinomial distribution is estimated. However, several distributions can 
#' be estimated from bins of \code{by_bin} width between \code{start_bin} and \code{end_bin}.
#' 
#' @param x SCE.
#' @param simf Numeric. Generate number of candidates * \code{simf} random samples.
#' @param gma Numeric. Multiply p by this number to generate alpha for Dirichlet-multinomial random sampling.
#' @param seedn Numeric. Number to set seed.
#' @param verbose Boolean
#'
#' @return SCE object.
#' @importFrom Matrix colSums Matrix
#' @export
sim_bg <- function(x, 
				   simf=0.5, 
				   gma=10, 
				   seedn=NULL, 
				   verbose=FALSE){
	if (verbose) cat("Simulating background droplets\n")

	tc <- Matrix::colSums(x@counts)

	set.seed(seedn)

	candd <- get_candd_counts(x)
	candd <- candd[x@de@deg,]
	nv <- nrow(candd)

	n_draw <- simf*ncol(candd)

	# Get p from DE test
	p <- x@de@low_means[x@de@deg]
	p <- p/sum(p)
	p <- p*gma
	# Generate alphas from gamma
	alphas <- matrix(nrow=nv, ncol=n_draw)
	rownames(alphas) <- rownames(x@de@deg)
	for (i in 1:n_draw){
		alphas[,i] <- rgamma(n=nv, shape=p)
	}

	# Get sizes
	candd_counts <- Matrix::colSums(candd)
	sizes_sampled <- sample(x=candd_counts, size=n_draw, replace=TRUE)

	# Set Dirichlet-multinomial parameters
	dmp <- DM_params(alphas=alphas, sizes=sizes_sampled, gamma=gma)
	x@sim@dm_params <- dmp

	# Draw from Dirichlet-multinomial
	sim <- draw_diri_multinom(dmp, seedn=seedn)
	sim <- Matrix::Matrix(sim, sparse=TRUE)
	all_expr <- cbind(sim, candd)
	# Fix labels, 0 is unknown, 1 is background
	labels <- c(rep(1, ncol(sim)), rep(0, ncol(candd)))
	x@sim@counts <- all_expr
	x@sim@labels <- labels

	if (verbose) cat("Simulated background droplets\n")
	return(x)
}
