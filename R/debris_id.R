
#' Initialize EM parameters for k=2 low and high groups
#'
#' Initialize the means of the two multinomial groups.
#' 
#' @param x SCE. SCE object with pi calculated
#' @param pct Numeric. Percent of droplets at each tail to fix. A number from 0 to 1.
#'
#' @return List with mu0 a p x 2 matrix and tau0 a length 2 numeric
#'
#' @importFrom Matrix colSums
#' @export
init_tails <- function(x, counts, pct=0.1){
	
	dc <- Matrix::colSums(x@counts)
	dc <- dc[colnames(counts)]
	n <- ncol(counts)
	n_tail <- round(n*pct)
	topn <- (names(dc)[order(dc, decreasing=TRUE)])[1:n_tail]
	botn <- (names(dc)[order(dc, decreasing=FALSE)])[1:n_tail]

	mu_top <- Matrix::rowSums(counts[,topn,drop=FALSE])
	mu_top <- mu_top/sum(mu_top)
	mu_bot <- Matrix::rowSums(counts[,botn,drop=FALSE])
	mu_bot <- mu_bot/sum(mu_bot)

	mu0 <- matrix(cbind(mu_top, mu_bot), ncol=2)
	colnames(mu0) <- c("Signal", "Background")
	rownames(mu0) <- rownames(counts)
	tau0 <- c(0.5, 0.5)
	ret <- list(Mu=mu0, Pi=tau0)
	return(ret)
}

#' x is a vector
#' @export
fraction_log <- function(x){
	x_c = x - max(x)
	x_c = exp(x_c)
	frac = x_c / sum(x_c);
	return(frac);
}

#' x is a vector
#' @export
sum_log <- function(x){
	max_x <- max(x)
	x_c = x - max_x
	x_sum <- log(sum(exp(x_c))) + max_x
	return(x_sum);
}

#' Get multinomial log prob of columns in X 
#'
#' @param X sparseMatrix. A sample by feature matrix of counts.
#' @param prob Numeric
#'
#' @import Matrix
#' @export
dmultinom_sparse <- function(X, prob){
	if (length(prob) != ncol(X)){
		stop("Length of prob and rows in X must be the same.")
	}
	X <- as(X, "dgCMatrix")
	Xlg <- X
	Xlg@x <- lgamma(Xlg@x + 1) # Get log gamma

	m <- lgamma(Matrix::rowSums(X) + 1)
	xls <- Matrix::rowSums(Xlg)

	px <- as.numeric(X %*% log(prob))

	prob <- m - xls + px
	return(prob)
}

#' Computed density of multinomial mixture for a matrix
#'
#' Given an n x p matrix of count data \code{x}, as well as a 
#' p x k matrix of probabilities \code{p}, return the probabilities 
#' of each sample for the k groups in a n x k matrix. If \code{x} 
#' has any non-integer values, they are rounded to the nearest integer.
#'
#' @param x sparseMatrix. A sample by feature matrix of counts.
#' @param p numeric matrix. A feature by group matrix of probabilities.
#' @param log Logical. If TRUE, return matrix with log probabilities.
#'
#' @return a numeric matrix with n samples by x variables.
#' @export
dmmn <- function(x, p, tau, labels=NULL){
	if (ncol(x) != nrow(p)) stop("Number of columns in x must be the same as the number of rows in p.")
	if (length(tau) != ncol(p)) stop("Length of tau must be the same as the number of columns in p.")
	if ( any(round(colSums(p), 1) != 1) ) stop("Probability columns in p must sum to 1.")
	if (round(sum(tau), 1) != 1) stop("Mixture probabilities in tau must sum to 1.")
	if ( any(Matrix::rowSums(x) == 0) ) stop("Sample rows of x must have at least 1 count.")

	n <- nrow(x)
	k <- ncol(p)

	if (is.null(labels)) labels <- numeric(n)
	if (length(labels) != n){
		stop("Size of labels must be the same as the number of rows in x.")
	}
	if ( any(is.na(labels)) ){
		stop("Labels must not contain any NA values.")
	}
	ks <- sort(unique(labels))
	if (max(ks) > k){
		stop("Labels must have a max integer value no more than the number of columns in p.")
	}
	ks <- ks[ks != 0]

	r <- matrix(nrow=n, ncol=k)
	colnames(r) <- colnames(p)
	rownames(r) <- rownames(x)
	for (ki in 1:k){
		r[,ki] <- dmultinom_sparse(X=x, prob=p[,ki])
		tau_v <- rep(x=log(tau[[ki]]), times=n) # Add log of mixture coefficient
		tau_v[(labels != 0) & (labels != ki)] <- -Inf
		r[,ki] <- r[,ki] + tau_v
	}
	return(r)
}

#' Computed log likelihood of multinomial mixture for a matrix
#'
#' Given an n x p matrix of count data \code{x}, as well as a 
#' p x k matrix of probabilities \code{p}, return the probabilities 
#' of each sample for the k groups in a n x k matrix. If \code{x} 
#' has any non-integer values, they are rounded to the nearest integer.
#'
#' @param x sparseMatrix. A sample by feature matrix of counts.
#' @param p numeric matrix. A feature by group matrix of probabilities.
#' @param log Logical. If TRUE, return matrix with log probabilities.
#'
#' @return a numeric matrix with n samples by x variables.
#' @export
dmmn_llk <- function(x, p, tau, probs=NULL, labels=NULL){
	if (is.null(probs)) probs <- dmmn(x=x, p=p, tau=tau, labels=labels)
	probs[apply(probs, 1, function(i) all(is.infinite(i))), ] <- c(0,0)
	probs_sum <- apply(probs, 1, sum_log)
	llk <- sum(probs_sum)
	return(llk)
}

#' Computed expected log likelihood of multinomial mixture
#'
#' Given an n x p matrix of count data \code{x}, as well as a 
#' p x k matrix of probabilities \code{p}, return the responsibilities 
#' of each sample for the k groups in a n x k matrix. If \code{x} 
#' has any non-integer values, they are rounded to the nearest integer.
#'
#' @param x sparseMatrix. A sample by feature matrix of counts.
#' @param p numeric matrix. A feature by group matrix of probabilities.
#' @param tau numeric. Mixture coefficients
#' @param log Logical. If TRUE, return matrix with log probabilities.
#'
#' @return a numeric matrix with n samples by k groups.
#' @export
e_step_mn <- function(x, p, tau, probs=NULL, labels=NULL){
	if (is.null(probs)) probs <- dmmn(x=x, p=p, tau=tau, labels=labels)
	probs[apply(probs, 1, function(i) all(is.infinite(i))), ] <- c(0,0)
	r <- t(apply(probs, 1, fraction_log))
	if (any(is.infinite(r))) stop("One or more responsibilities are infinite, probably from dividing by 0.")
	return(r)
}

#' Maximize multinomial parameters in EM
#'
#' Given an n x p matrix of count data, as well as a 
#' n x k matrix of responsibilities, return the maximum 
#' likelihood estimate (MLE) of p for each k group in 
#' a p x k matrix.
#'
#'
#' @export
m_step_mn <- function(x, r, labels=NULL){
	n <- nrow(x)
	k <- ncol(r)
	# Semi-supervised
	if (!is.null(labels)){
		if (length(labels) != nrow(r)){
			stop("Size of labels must be the same as the number of rows in r.")
		}
		ks <- sort(unique(labels))
		if (max(ks) > k){
			stop("Labels must have a max integer value no more than the number of columns in p.")
		}
		ks <- ks[ks != 0]
		for (k in ks){
			p <- numeric(ncol(r))
			p[k] <- 1
			r[labels == k, ] <- p
		}
	}
	wm <- as.matrix(Matrix::t(x) %*% r)
	Mu <- sweep(wm, 2, colSums(wm), "/")
	tau <- colSums(r)
	tau <- tau/sum(tau)
	return(list(Mu=Mu, Pi=tau))
}

#' Run EM on pi features
#' 
#' Run expectation maximization (EM) using diagonal covariance on pi matrix in an SCE object with k=2 groups. 
#' EM is run \code{n_runs} times with random initializations, and parameters with 
#' the best log likelihood are output. Within EM, iterate a minimum of \code{min_iter} times and 
#' a maximum of \code{max_iter} times. The \code{eps} value gives the threshold of percentage change in 
#' log likelihood until convergence is reached.
#'
#' @param x SCE object
#'
#' @return SCE object with EM output in the \code{emo} slot. See \code{\link{run_mv_em_diag}} for details
#' @importFrom EMCluster emcluster
#' @importFrom mvtnorm dmvnorm
#' @export
run_em <- function(x, eps=1e-8, max_iter=1e3, labels=NULL, verbose=TRUE){
	
	# Store EM output
	emo <- list()
	x@diem@converged  <- FALSE

	genes <- x@diem@de@diff_genes
	drops <- x@test_IDs
	counts <- Matrix::Matrix(x@counts[genes, drops], sparse=TRUE)

	if (verbose){
		cat("Running EM on ", as.character(nrow(counts)), " genes and ", as.character(ncol(counts)), " droplets.\n")
	}

	no_counts <- Matrix::colSums(counts) == 0
	if (sum(no_counts) > 0){
		cat(paste0("Warning: removing ", as.character(sum(no_counts)), " droplets with no diff genes detected.\n"))
		x@test_IDs <- x@test_IDs[!no_counts]
		counts <- counts[,!no_counts]
		labels <- labels[!no_counts]
	}

	# Run EM
	mn_params <- init_tails(x, counts)
	loglk <- -Inf
	loglks <- c(loglk)
	iter <- 1
	counts <- Matrix::t(counts)
	probs <- dmmn(x=counts, p=mn_params$Mu, tau=mn_params$Pi, labels=labels)
	while (iter < max_iter){
		# E step
		resp <- e_step_mn(x=counts, p=mn_params$Mu, tau=mn_params$Pi, probs=probs, labels=labels)

		# Store in EMO
		# M step
#		bg_IDs <- base::setdiff(colnames(x@counts), x@test_IDs)
#		bg_counts <- x@counts[,bg_IDs]
#		bg_sums <- Matrix::rowSums(bg_counts)
#		respm <- as.matrix(data.frame(Signal=rep(0, ncol(x@counts)), Background=rep(1, ncol(x@counts))), ncol=2)
#		rownames(respm) <- colnames(x@counts)
#		respm[rownames(resp),] <- resp
#		labelsm <- rep(2, ncol(x@counts)); names(labelsm) <- colnames(x@counts)
#		labelsm[names(labels)] <- labels
#		mn_params <- m_step_mn(x=Matrix::t(x@counts[genes,]), r=respm)
		mn_params <- m_step_mn(x=counts, r=resp)
		probs <- dmmn(x=counts, p=mn_params$Mu, tau=mn_params$Pi, labels=labels)

		# Evaluate log likelihood
		loglk <- loglks[[iter]] <- dmmn_llk(x=counts, p=mn_params$Mu, tau=mn_params$Pi, probs=probs, labels=labels)

		emo[[iter]] <- list(Z=resp, Mu=mn_params$Mu, Pi=mn_params$Pi, loglk=loglk)

		if (verbose){
			cat(paste0("Iteration ", as.character(iter), "; llk ", as.character(loglk), "\n"))
		}

		if (iter > 1){
			dloglk <- (loglks[[iter]] - loglks[[iter-1]])/abs(loglks[[iter-1]])
			if (dloglk < 0){
				cat("Warning: Likelihood decreased.\n")
			}
			if (dloglk < eps){
				x@diem@converged  <- TRUE
				if (verbose) cat("Converged!\n")
				break
			}
		}
		iter <- iter + 1
	}
	probs <- dmmn(x=counts, p=mn_params$Mu, tau=mn_params$Pi)
	x@diem@PP <- e_step_mn(x=counts, p=mn_params$Mu, tau=mn_params$Pi, probs=probs) # No labels

	x@diem@assign <- c("Signal"=1, "Background"=2)

	x@diem@emo <- emo
	if (verbose){
		cat("Finished EM\n")
	}
	return(x)
}

#' Call targets after EM
#'
#' Call targets from droplets if the log-likelihood membership 
#' probability is higher than \code{lk_fraction} given during 
#' initialization.
#'
#' @param x SCE. SCE object.
#' @param interation Integer. DIEM iteration number to use for calliing targets
#'
#' @return SCE object
#' @export
call_targets <- function(x, pp_thresh=0.95, min_genes=200){

	Z <- x@diem@PP
	asgn <- x@diem@assign

	calls <- rep("Debris", nrow(Z))
	names(calls) <- rownames(Z)

	tb <- (Z[,asgn["Signal"]] > pp_thresh) & (x@dropl_info[rownames(Z),"n_genes"] >= min_genes)
	target_names <- rownames(Z)[tb]
	calls[tb] <- "Signal"

	x@diem@calls <- calls

	return(x)
}

#' Return target IDs
#'
#' @param x SCE. SCE object
#'
#' @return Character vector with droplet IDs of called targets
#' @export
targets_ids <- function(x){
	calls <- x@diem@calls

	return(names(calls)[calls == "Signal"])
}

#' Insert pi and calls to dropl_info
#'
#' @param x SCE. SCE object
#'
#' @return Character vector with droplet IDs of called targets
#' @export
fill_dropl_info <- function(x){
	# Fill calls
	x@dropl_info[names(x@diem@calls),"Call"] <- x@diem@calls

	# Fill PP
	x@dropl_info[rownames(x@diem@PP),"PP"] <- x@diem@PP[,1,drop=FALSE]

	return(x)

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

#' Get percent of reads aligning to given genes
#'
#' @param x SCE.
#' @param genes Character. Genes to calculate percentage of in counts
#' @param name Character. Column name to place in dropl_info
#'
#' @return SCE object
#' @export
get_gene_pct <- function(x, genes, name){
	n_drop <- nrow(x@dropl_info)
	expr <- x@counts[genes,]
	if (length(expr) == 0){
		stop("None of genes found in counts.")
	}
	gene_pct <- Matrix::colSums(x@counts[genes,]) / x@dropl_info[,"total_counts"]
	x@dropl_info[names(gene_pct),] <- gene_pct
	return(x)
}
