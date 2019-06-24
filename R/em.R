
#' Initialize parameters for EM using K-means
#'
#' @param x Matrix or data frame. Data observations with n rows and p columns
#' @param k Integer. Number of clusters
#' @param labels Vector. Integer vector of labels. If semi-supervised, samples with 0 are 
#' unknown, wheres positive integers are fixed labels
#'
#' @return Parameters in a list with components
#' \item{mu}{List of length k, with each element a numeric vector of mus of length p}
#' \item{sgma}{List of length k, with each element a numeric vector of covariance}
#' \item{tau}{Numeric vector of length k}
#' @export
init_param <- function(x, k, labels=NULL){
	n <- nrow(x)
	p <- ncol(x)

	if (is.null(labels)){
		km <- kmeans(x, centers=k,nstart=10)
		mu <- list()
		sgma <- list()
		for (i in 1:k){
			mu[[i]] <- km$centers[i,]
			var_i <- var(x[km$cluster == i,])
			if (p > 1) sgma[[i]] <- diag(var_i)
			else sgma[[i]] <- var_i
		}
		tau <- sapply(1:k, function(l) sum(km$cluster == l)/n)
		return(list(mu=mu, sgma=sgma, tau=tau))
	} else {
		all_labels <- labels # For storing both known and k-means clusters
		labels_known <- setdiff(sort(unique(labels)), 0)
		labels_unknown <- setdiff(1:k, labels_known)

		mu <- list()
		sgma <- list()

		km <- kmeans(x, centers=k, nstart=10)
		centers <- km$centers
		kclusters <- km$cluster
		kclusters_av <- 1:k
		for (i in labels_known){
			cmean <- colMeans(x[labels == i,,drop=FALSE])
			# Find closest of kclusters_av
			centers_m <- as.matrix(centers[kclusters_av,], nrow=length(kclusters_av))
			distances <- dist( rbind(cmean, centers_m) )
			nearest_kclust <- kclusters_av[which.min(distances[1:length(kclusters_av)])]
			mu[[i]] <- centers[nearest_kclust,]
			var_i <- var(x[kclusters == nearest_kclust,])
			if (p > 1) {
				sgma[[i]] <- diag(var_i)
			} else {
				sgma[[i]] <- var_i
			}
			all_labels[kclusters == nearest_kclust] <- i
			kclusters_av <- setdiff(kclusters_av, nearest_kclust)
		}
		for (i in seq_along(labels_unknown)){
			mu[[labels_unknown[i]]] <- centers[kclusters_av[i],]
			var_i <- var(x[kclusters == kclusters_av[i],])
			if (p > 1) sgma[[labels_unknown[i]]] <- diag(var_i)
			else sgma[[labels_unknown[i]]] <- var_i
			all_labels[kclusters == kclusters_av[i]] <- labels_unknown[i]
		}
		tau <- sapply(1:k, function(l) sum(all_labels == l)/n)
		return(list(mu=mu, sgma=sgma, tau=tau))
	}
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
#' 
#' For compatibility with Rcpp
#'
#' @param x Numeric.
#'
#' @return A numeric (same as x).
set_double <- function(x){
	prec <- .Machine[["double.xmin"]]
	x <- sapply(x, function(i) {
				if (i > 0 & i < prec) i <- prec
				return(i)})
	return(x)
}

#' Run multivariate EM restricted to diagonal covariance
#'
#' Run multivariate un-, semi-, or fully supervised EM.
#'
#' @param x	matrix. Matrix with samples on the rows and features on the columns.
#' @param k Numeric. Number of clusters. Must be pre-specified.
#' @param min_iter Numeric. Minimum number of iterations before calling convergence.
#' @param max_iter Numeric. Maximum number of iterations.
#' @param labels Numeric vector. Fixed label IDs for semi-supervised clustering. 0 is unknown (to be estimated).
#' @param mu0 Numeric. Starting value for mu. If set, all mu0, sgma0, and tau0 must be set.
#' @param sgma0 Numeric. Starting value for sigma. If set, all mu0, sgma0, and tau0 must be set.
#' @param tau0 Numeric. Starting value for tau. If set, all mu0, sgma0, and tau0 must be set.
#' @param eps Numeric. Threshold of change in log likelihood to stop iterations.
#' @param seedn Numeric. Numeric value specifying seed for random initialization.
#' @param verbose Boolean. Print logging information.
#' @param trace Boolean. Print out parameter and likelihod values during iterations.
#'
#' @return EMO object
#' @export
run_mv_em_diag <- function(x, 
						   k=2, 
						   min_iter=5, 
						   max_iter=1000, 
						   labels=NULL, 
						   mu0=NULL, 
						   sgma0=NULL, 
						   tau0=NULL, 
						   eps=1e-10, 
						   seedn=NULL, 
						   verbose=TRUE, 
						   trace=FALSE){
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
	set.seed(seedn)
	params <- init_param(x, k, labels)
	mu <- params$mu
	sgma <- params$sgma
	tau <- params$tau
	if (!is.null(mu0)){
		mu <- mu0
		sgma <- sgma0
		tau <- tau0
	}
	llk <- -Inf
	llks <- c()
	
	# Debugging
	if (trace){
		cat("Initialized parameters:\n")
		cat("sigma:\n"); print(sgma);cat("mu:\n");  print(mu); cat("tau:\n"); print(tau)
	}
	# Iterations
	converged <- FALSE
	n_iter <- NA
	for (iter in 1:max_iter){
		# Populate responsibilities (E-step)
		Z <- e_stepCPP(x, k, mu, sgma, tau, semisup, labels)
		
		# Maximize observed-data log-likelihood (M-step)
		mp <- m_step(x, Z)
		mu_n <- mp$mu
		sgma_n <- mp$sgma
		tau_n <- mp$tau

		# Debugging
		if (trace){
			cat("sigma:\n"); print(mu_n); cat("mu:\n"); print(sgma_n); cat("tau:\n"); print(tau_n)
		}

		# Set double precision to match with C++ if x is below double minimum
		mu_n <- lapply(mu_n, set_double)
		sgma_n <- lapply(sgma_n, set_double)
		tau_n <- set_double(tau_n)

		# These 2 checks exit the iterations when 1 group contains either 0 or 1 observation.
		# Should this be changed to return NA, instead of converged?
		if (any(unlist(sgma_n) == 0) | any(tau_n == 0) ) {
			tau <- tau_n
			mu <- mu_n
			sgma <- sgma_n
			llks[iter] <- llk
			if (verbose) cat("Warning: converged after 1 or less observation in a group. Can't calculate likelihood\n")
			print(mp)
			converged <- TRUE
			n_iter <- iter
			break
		}

		llk_n <- get_llkCPP(x, mu_n, sgma_n, tau_n, semisup=semisup, labels);

		if (is.na(llk_n)) stop("Returned likelihood is NA")
		llk_diff <- (llk_n - llk) / abs(llk_n)
		if (trace){
			cat(paste0("Iteration ", as.character(iter), ", llk: ", as.character(llk_n), '\n'))
			if (llk_diff < 0) cat("Likelihood decreased during EM iteration. Probably a bug...")
		}
		
		# Update paramters
		mu <- mu_n
		sgma <- sgma_n
		tau <- tau_n
		llk <- llk_n
		llks[iter] <- llk
		
		# Check if likelihood difference doesn't change
		if (iter >= min_iter & abs(llk_diff) < eps) {
			if (trace){
				cat(paste0("Converged after ", as.character(iter), " iterations.\n"))
			}
			converged <- TRUE
			n_iter <- iter
			break
		}
	}

	# Return
	if (converged){
		rownames(Z) <- rownames(x)
		for (i in 1:k){
			names(mu[[i]]) <- colnames(x)
			names(sgma[[i]]) <- colnames(x)
		}
		ret <- EMO(Z = Z,
				   mu = mu,
				   sgma = sgma,
				   tau = tau,
				   llks = llks,
				   n_iter = n_iter)
	}
	else{
		ret <- EMO(Z = NA,
				   mu = NA,
				   sgma = NA,
				   tau = NA,
				   llks = NA, 
				   n_iter = n_iter)
	}
	if (verbose) cat("Finished EM\n")
	return(ret)
}

