
#' Initialize parameters for EM
#'
#' @param x Matrix or data frame. Data observations with n rows and p columns
#' @param k Integer. Number of clusters
#'
#' @return Parameters in a list with components
#' \item{mu}{List of length k, with each element a numeric vector of mus of length p}
#' \item{sgma}{List of length k, with each element a numeric vector of covariance}
#' \item{tau}{Numeric vector of length k}
init_param <- function(x, k, labels=NULL){
	n <- nrow(x)
	p <- ncol(x)
	if (is.null(labels)){
		mins <- apply(x, 2, min)
		maxs <- apply(x, 2, max)
	} else {
		mins <- apply(x[labels == 0,,drop=FALSE], 2, min)
		maxs <- apply(x[labels == 0,,drop=FALSE], 2, max)
	}
	mu <- list()
	for (i in 1:k){
		mu[[i]] <- rep(NA, p)
		for (j in 1:p){
			mu[[i]][j] <- runif(n=1, min=mins[[j]], max=maxs[[j]])
		}
	}
	xvar <- apply(x, 2, var)
	sgma <- lapply(1:k, function(x) xvar)
	# sgma <- lapply(1:k, rep, x=1, length.out=p)
	tau <- rep(1/k, k)
	# If labeled
	if (!is.null(labels)){
		label_k <- sort(unique(labels))
		for (i in label_k){
			if (i == 0) next
			mu[[i]] <- apply(x[labels == i,,drop=FALSE], 2, mean)
			sgma[[i]] <- apply(x[labels == i,,drop=FALSE], 2, var)
		}
	}

	return(list(mu=mu, sgma=sgma, tau=tau))
}

#' M Step
#'
#' @param x Matrix or data frame. Data observations with n rows and p columns
#' @param Z
#' 
#' @return list with updated mu, sigma, and tau parameters
m_step <- function(x, Z){
	n <- nrow(x)
	p <- ncol(x)
	k <- ncol(Z)
	tau <- colSums(Z)/n
	mu <- lapply(1:k, function(i) as.vector((t(Z[,i]) %*% x)/sum(Z[,i])) )
	sgma <- lapply(1:k, function(i) {
				   xc <- apply(x, 1, function(j) j - mu[[i]])
				   if (p > 1) xc <- t(xc)
				   covar <- t(Z[,i]) %*% (xc * xc) / sum(Z[,i])
				   return(as.vector(covar))
	})
	return(list(mu=mu, sgma=sgma, tau=tau))
}

#' Round x up to double precision value if R rounds x below the minimum number for double
set_double <- function(x){
	prec <- .Machine[["double.xmin"]]
	for (i in 1:length(x)){
		if (x[i] > 0 & x[i] < prec){
			x[i] <- prec
		}
	}
	return(x)
}

#' Run multivariate EM restricted to diagonal covariance
#'
#' Run multivariate un-, semi-, or fully supervised EM.
#'
#' @param x	matrix. Matrix with samples on the rows and features on the columns
#' @param k Numeric. Number of clusters. Must be pre-specified
#' @param min_iter Numeric. Minimum number of iterations before calling convergence
#' @param max_iter Numeric. Maximum number of iterations
#' @param labels Numeric vector. Fixed label IDs for semi-supervised clustering. 0 is unknown (to be estimated)
#' @param mu0 Numeric. Starting value for mu
#' @param sgma0 Numeric. Starting value for sigma
#' @param eps Numeric. Threshold of change in log likelihood to stop iterations
#' @param seedn Numeric. Numeric value specifying seed for random initialization
#' @param verbose Boolean. Print logging information
#'
#' @return list with output
#' @export
run_mv_em_diag <- function(x, 
						   k=2, 
						   min_iter=5, 
						   max_iter=1000, 
						   labels=NULL, 
						   mu0=NULL, 
						   sgma0=NULL, 
						   eps=1e-10, 
						   seedn=NULL, 
						   verbose=TRUE){
	x <- as.matrix(x)
	p <- ncol(x)
	n <- nrow(x)
	
	if (is.null(labels)) {
		labels <- integer(n)
		semisup <- FALSE
		if (verbose) cat("Running unsupervised EM\n")
	} else {
		labels <- as.integer(labels)
		semisup <- TRUE
		if (verbose) cat("Running semisupervised EM\n")
	}
	if (any(labels > k) | any(labels < 0)) stop("labels must only contain integers in 0 through k")
	if (!is.integer(labels)) stop("labels must be a vector integers")
	if (length(labels) < n) stop(paste0("Length of labels must be same as number of rows of x, not ", str(length(labels))))
	
	# 1) Initialize parameters
	# Initialize mu with random uniform from min to max of each variable p
	set.seed(seedn)
	params <- init_param(x, k, labels)
	mu <- params$mu
	sgma <- params$sgma
	tau <- params$tau
	if (!is.null(mu0)) mu <- mu0
	if (!is.null(sgma0)) sgma <- sgma0
	llk <- -Inf
	llks <- c()
	
	# Iterations
	converged <- FALSE
	n_iter <- NULL
	for (iter in 1:max_iter){
		# Populate Z membership probs
		Z <- e_stepCPP(x, k, mu, sgma, tau, semisup, labels)
		
		# Maximize observed-data log-likelihood
		mp <- m_step(x, Z)
		mu_n <- mp$mu
		sgma_n <- mp$sgma
		tau_n <- mp$tau

		# These 2 checks exit the iterations when 1 group contains either 0 or 1 observation.
		# Should this be changed to return NA, instead of converged?
		if (any(unlist(sgma_n) == 0) | any(tau_n == 0) ) {
			tau <- tau_n
			mu <- mu_n
			sgma <- sgma_n
			llks[iter] <- llk
			cat("Warning: converged after 1 or less observation in a group. Can't calculate likelihood\n")
			print(mp)
			if (verbose) cat(paste0("Converged after ", as.character(iter), " iterations.\n"))
			converged <- TRUE
			n_iter <- iter
			break
		}

		# Set double precision to match with C++ if x is below double minimum
		mu_n <- lapply(mu_n, set_double)
		sgma_n <- lapply(sgma_n, set_double)
		tau_n <- set_double(tau_n)

		# Evaluate data likelihood if possible
		if (any(unlist(sgma_n) == 0) | any(tau_n == 0) ) {
			llk_n <- llk
		} else {
			llk_n <- get_llkCPP(x, mu_n, sgma_n, tau_n, semisup=semisup, labels);
		}

		if (is.na(llk_n)) stop("Returned likelihood is NA")
		llk_diff <- (llk_n - llk) / abs(llk_n)
		if (verbose) cat(paste0("Iteration ", as.character(iter), ", llk: ", as.character(llk_n), '\n'))
		# if (llk_diff < 0) stop("Likelihood decreased during EM iteration. Probably a bug...")
		
		# Update paramters
		tau <- tau_n
		mu <- mu_n
		sgma <- sgma_n
		llk <- llk_n
		llks[iter] <- llk
		
		# Check if likelihood difference doesn't change
		if (iter >= min_iter & abs(llk_diff) < eps) {
			if (verbose) cat(paste0("Converged after ", as.character(iter), " iterations.\n"))
			converged <- TRUE
			n_iter <- iter
			break
		}
	}
	rownames(Z) <- rownames(x)

	# Return
	if (converged){
		ret <- list(Z = Z,
					mu = mu,
					sgma = sgma,
					tau = tau,
					llks = llks,
					n_iter = n_iter)
	}
	else{
		ret <- list(Z = NULL,
					mu = NULL,
					sgma = NULL,
					tau = NULL,
					llks = NULL, 
					n_iter = n_iter)
	}
	return(ret)
}

##' Log-likelihod of multivariate normal with diagonal covariance
##'
##' Calculate log-likelihood of sample x given mu and diagonal sigma.
##'
##' @param x numeric vector or matrix. The data observations
##' @param mu numeric vector. Mu average of each variable in the multivariate normal
##' @param sgma numeric vector or diagonal matrix. If matrix, extracts diagonal only.
##'
##' @return numeric value of the log likelihood of the given observation x under
##'  parameters mu and sgma
#mvn_logllk_diag <- function(x, mu, sgma){
#	# If sgma is a matrix, change to vector
#	if (inherits(sgma, "matrix") | inherits(sgma, "Matrix")) {
#		sgma = diag(sgma)
#	}
#	x <- as.vector(x)
#	mu <- as.vector(mu)
#	sgma <- as.vector(sgma)
#	same_len = (length(x) == length(mu)) & (length(mu) == length(sgma))
#	if (!same_len) stop("x, mu, and sgma must have the same length")
#
#	logdet <- sum(log(sgma))
#	
#	sgma_inv <- 1/sgma
#	x <- x - mu
#	xe <- x %*% (sgma_inv * x) # Xt %*% sgma_INV %*% X
#	
#	con <- length(x) * log(2*pi)
#	
#	ret <- -0.5 * (logdet + xe + con)
#	return(ret)
#}
#
##' Fraction of logs
##'
##' Given log values x1, .., xn, calculate for each value, the relative fraction
##' e.g. x1 / (x1 + ... + xn)
##'
##' @param x numeric vector. Values in log space
##'
##' @return numeric vector of relative fractions of logs
#fraction_log <- function(x){
#	x <- x - max(x)
#	x <- exp(x)
#	frac <- x / sum(x)
#	return(frac)
#}
#
##' Sum of logs
##'
##' Given log values x1, .., xn, calculate the sum: 
##' log( sum( exp(x1) + ... + exp(xn) ) )
##'
##' @param x numeric vector. Values in log space
##'
##' @return numeric value of the sum of the logs
#sum_log <- function(x){
#	max_x <- max(x)
#	x <- x - max(x)
#	x <- log(sum(exp(x))) + max_x
#	return(x)
#}
#
## Function to get expected likelihood under mvn model
## x is a sample x feature matrix (n x f)
## mu is a list of length k containing the vector of averages for each k cluster
## sgma is a list of length k containing the diagonals of the covariance matrix for each cluster
## tau is a vector of length k containing the mixing parameter for each cluster k
#
##' log likelihood of multivariate normal under 1 or more latent clusters.
##'
##' Calculates log likelihood of multivariate distribution given data and parameters. 
##' For each of k clusters, calculate log likelihood, multiply by its mixing parameter tau, and 
##' sum across clusters.
##'
##' @param x matrix or data frame. Observations with which to calculate likelihood under parameters.
##' @param mu list of numeric vectors. Averages for each of the k clusters.
##' @param sgma list of numeric vectors. Covariance for each of the k clusters.
##' @param tau numeric vector. Mixing parameters for each of the k clusters.
##' 
##' @return numeric value of the observed log likelihood under k clusters.
#get_llk <- function(x, mu, sgma, tau){
#	eps <- 0.01
#	if (sum(tau) < (1 - eps) | sum(tau) > (1+eps)) stop(paste0("Parameter tau must sum to 1, not ", as.character(tau)))
#	same_len <- (length(mu) == length(sgma)) & (length(sgma) == length(tau))
#	if (!same_len) stop("mu, sgma, and tau must have same length")
#
#	n <- nrow(x)
#	k <- length(mu)
#	llk <- 0
#	for (i in seq(n)){
#		llks <- c()
#		for (j in seq(k)){
#			if (i < 11){
#				print(mvn_logllk_diag(x[i,], mu[[j]], sgma[[j]]))
#			}
#			llks[j] <- log(tau[j]) + mvn_logllk_diag(x[i,], mu[[j]], sgma[[j]])
#		}
#		llk <- llk + sum_log(llks)
#		if (i < 11){
#			# print(llk)
#		}
#	}
#	return(llk)
#}
#
##' Semisupervised log likelihood of multivariate normal with labeled and unlabeled data
##'
##' Calculates log likelihood of multivariate distribution given data and parameters. 
##' For each of k clusters, calculate log likelihood, multiply by its mixing parameter tau, and 
##' sum across clusters.
##'
##' @param x matrix or data frame. Observations with which to calculate likelihood under parameters.
##' @param mu list of numeric vectors. Averages for each of the k clusters.
##' @param sgma list of numeric vectors. Covariance for each of the k clusters.
##' @param tau numeric vector. Mixing parameters for each of the k clusters.
##' @param labels numeric vector. Labels of observations. 0 for unknown.
##' 
##' @return numeric value of the observed log likelihood under k clusters.
#get_llk_ss <- function(x, mu, sgma, tau, labels){
#	eps <- 0.01
#	if (sum(tau) < (1 - eps) | sum(tau) > (1+eps)) stop(paste0("Parameter tau must sum to 1, not ", as.character(tau)))
#	same_len <- (length(mu) == length(sgma)) & (length(sgma) == length(tau))
#	if (!same_len) stop("mu, sgma, and tau must have same length")
#
#	llk <- 0
#
#	obs <- x[labels != 0, , drop=FALSE]
#	obs_labels <- labels[labels != 0]
#	for (i in 1:nrow(obs)){
#		j <- obs_labels[i]
#		llk <- llk + log(tau[j]) + mvn_logllk_diag(obs[i,], mu[[j]], sgma[[j]])
#	}
#
#	unobs <- x[labels == 0, , drop=FALSE]
#	n <- nrow(unobs)
#	k <- length(mu)
#	for (i in seq(n)){
#		llks <- c()
#		for (j in seq(k)){
#			llks[j] <- log(tau[j]) + mvn_logllk_diag(unobs[i,], mu[[j]], sgma[[j]])
#		}
#		llk <- llk + sum_log(llks)
#	}
#	return(llk)
#}
#
##' Expected log likelihood of multivariate normal under 1 or more latent clusters.
##'
##' Calculates expected log likelihood of multivariate distribution given data and parameters. 
##' For each of k clusters, calculate log likelihood, multiply by its mixing parameter tau, and 
##' sum across clusters.
##'
##' @param x matrix or data frame. Observations with which to calculate likelihood under parameters.
##' @param mu list of numeric vectors. Averages for each of the k clusters.
##' @param sgma list of numeric vectors. Covariance for each of the k clusters.
##' @param tau numeric vector. Mixing parameters for each of the k clusters.
##' @param Z matrix or data frame. Latent membership probabilities. n x k matrix.
##' 
##' @return numeric value of the expected log likelihood under k clusters.
#get_exp_llk <- function(x, mu, sgma, tau, Z){
#	eps <- 0.01
#	if (sum(tau) < (1 - eps) | sum(tau) > (1+eps)) stop(paste0("Parameter tau must sum to 1, not ", as.character(tau)))
#	same_len <- (length(mu) == length(sgma)) & (length(sgma) == length(tau))
#	if (!same_len) stop("mu, sgma, and tau must have same length")
#
#	n <- nrow(x)
#	k <- length(mu)
#	exp_llk <- 0
#	for (i in seq(n)){
#		for (j in seq(k)){
#			this_llk <- Z[i,j] * (log(tau[j]) + mvn_logllk_diag(x[i,], mu[[j]], sgma[[j]]))
#			exp_llk <- exp_llk + this_llk
#		}
#	}
#	return(exp_llk)
#}
#
##' E Step
##'
##' @param x Matrix or data frame. Data observations with n rows and p columns
##' @param k Integer. Number of clusters
##' @param mu List of numeric vectors. Contains averages
##' @param sgma List. List of numeric vectors containing diagonal elements of sigma
##' @param tau Numeric. A vector of length \code{2}; the mixing parameter
##' @param semisup Boolean. Semi-supervised EM
##' @param labels Integer. If \code{semisup}, then give the labels
##' 
##' @return Z numeric matrix with membership probabilities
#e_step <- function(x, k, mu, sgma, tau, semisup=FALSE, labels=NULL){
#	n <- nrow(x)
#	Z <- matrix(0, nrow=n, ncol=k)
#	for (i in 1:n){
#		# Get log likelihood for each of k clusters
#		f_i <- rep(-Inf, k)
#		for (j in 1:k){
#			f_i[j] <- log(tau[[j]]) + mvn_logllk_diag(x[i,], mu[[j]], sgma[[j]])
#		}
#		if (semisup & labels[i] != 0){
#			Z[i,] <- rep(0,k)
#			Z[i,labels[i]] <- 1
#		} else {
#			Z[i,] <- fraction_log(f_i)
#		}
#	}
#	return(Z)
#}
