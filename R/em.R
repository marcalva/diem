
#' Get log likelihood under Dirichlet-multinomial
#'
#' Get log likelihood from \code{counts} under a Dirichlet-multinomial 
#' with parameter \code{Alpha}.
#'
#' @param counts A gene by droplet data matrix of read counts.
#' @param Alpha A gene by cluster matrix of parameter values for the 
#'  Dirichlet-multinomial. Each entry must be greater than 0.
#' @param droplets Optional vector of droplet IDs the likelihood is 
#'  calculated for (instead of all droplets in \code{counts}).
#' @param sizes Optional vector of droplet (column) sizes.
#' @param threads Number of threads for parallel execution. Default is 1.
#' @param verbose Verbosity.
#' 
#' @return A droplet by cluster matrix of likelihood values
#' 
#' @importFrom Matrix colSums
#'
get_llk <- function(counts, Alpha, droplets = NULL, sizes = NULL, threads = 1, verbose = FALSE){

    if (any(is.na(Alpha))){
        stop("Alpha has NA value(s)")
    }
    if (any(Alpha <= 0)){
        stop("Alpha has value(s) <= 0")
    }

    if (is.null(sizes)){
        sizes <- colSums(counts)
    }
    if (verbose) message("estimating likelihoods")
    if (is.null(droplets)){
        llk <- LlkDirMultSparsePar(counts, sizes, alpha = Alpha, threads = threads, display_progress = verbose)
    } else {
        llk <- LlkDirMultSparsePar(counts[,droplets], sizes[droplets], alpha = Alpha, threads = threads, display_progress = verbose)
    }
    rownames(llk) <- colnames(counts)
    return(llk)
}

#' Estimate alpha using method of moments
#' 
#' Estimate Alpha parameters of DM using method of moments. The input 
#' data \code{counts} is in sparse matrix format. The parameters for 
#' DM can be estimated for multiple centers, where both \code{clust_mem} and 
#' \code{weights} specify the cluster's members and their weights, 
#' respectively.
#'
#' @param counts A gene by droplet sparse count matrix.
#' @param clust_mem A droplet by cluster boolean matrix, indicating whether 
#'  a droplet will count towards estimation of the cluster's parameters.
#' @param weights A droplet by cluster matrix of weights, or the posterior 
#'  probability the droplet belongs to the cluster
#' @param psc pseudocount to add to the final alpha. This avoids alpha values 
#'  of 0 that are outside the domain of the parameter.
#' @param threads Number of threads for parallel execution. Default is 1.
#'
#' @return A gene by cluster matrix of alpha values for the cluster's 
#'  DM parameters
#'
#' @importFrom Matrix t
alpha_mom <- function(counts, 
                      clust_mem, 
                      weights, 
                      psc = 1e-10, 
                      threads = 1){
    clusts <- 1:ncol(clust_mem)
    p <- divide_by_colsum(counts)
    p_bar <- sapply(clusts, function(k){
                    pb <- fast_wmeanCPP(x = t(p[,clust_mem[,k],drop=FALSE]), 
                                        weights = weights[clust_mem[,k], k], 
                                        threads = threads)
                    return(pb) })

    p_var <- sapply(clusts, function(k){
                    pv <- fast_wvarCPP(x = t(p[,clust_mem[,k],drop=FALSE]), 
                                       mu = p_bar[,k], 
                                       weights = weights[clust_mem[,k], k], 
                                       threads = threads)
                    return(pv) })

    a0 <- sapply(clusts, function(k){
                 nm <- p_bar[,k] * (1 - p_bar[,k])
                 a <- (nm / p_var[,k]) - 1
                 keep <- !( is.infinite(a) | is.na(a) )
                 a[!keep] <- 0
                 a[ a < 0 ] <- 0
                 return(max(a))
                    })
    alpha0 <- p_bar %*% diag(a0)
    alpha0 <- alpha0 + psc
    return(alpha0)
}

#' Maximize leave-one-out with respect to alpha
#' 
#' Estimate Alpha parameters of DM using method of moments. The input 
#' data \code{counts} is in sparse matrix format. The parameters for 
#' DM can be estimated for multiple centers, where both \code{clust_mem} and 
#' \code{weights} specify the cluster's members and their weights, 
#' respectively.
#'
#' @param counts A gene by droplet sparse count matrix.
#' @param alpha0 A gene by cluster matrix of initial alpha parameter values.
#' @param clust_mem A droplet by cluster boolean matrix, indicating whether 
#'  a droplet will count towards estimation of the cluster's parameters.
#' @param weights A droplet by cluster matrix of weights, or the posterior 
#'  probability the droplet belongs to the cluster
#' @param psc pseudocount to add to the final alpha. This avoids alpha values 
#'  of 0 that are outside the domain of the parameter.
#' @param eps The threshold of change in the parameter to stop the 
#'  LOO iterations.
#' @param max_iter The maximum number of LOO iterations.
#' @param threads Number of threads for parallel execution. Default is 1.
#'
#' @return A gene by cluster matrix of alpha values for the cluster's 
#'  DM parameters
#'
#' @importFrom Matrix colSums t
alpha_max_loo <- function(counts, 
                          alpha0, 
                          clust_mem, 
                          weights, 
                          psc = 1e-10, 
                          eps = 1e-4, 
                          max_iter = 1e4, 
                          threads = 1){
    sizes <- colSums(counts)
    clusts <- 1:ncol(clust_mem)

    # LOO optimize
    Alpha <- sapply(clusts, function(k){
                    Ak <- max_loo(x = t(counts[,clust_mem[,k],drop=FALSE]), 
                                  sizes = sizes[clust_mem[,k]], 
                                  weights = weights[clust_mem[,k], k], 
                                  alpha = alpha0[,k], 
                                  eps = eps, 
                                  psc = psc, 
                                  threads = threads, 
                                  max_iter = max_iter)
                    return(Ak)
                          })
    rownames(Alpha) <- rownames(counts)
    return(Alpha)
}

#' Get MLE of alphas
#' 
#' Given the count data and posterior probabilities (membership), 
#' update alpha by finding the MLE. This initializes alpha to 
#' the mean proportion of the cluster's gene expression and 
#' maximizes the likelihood with respect to the precision (the sum 
#' of the alphas). Then, the leave-one-out likelihood is maximized, 
#' providing the final MLE.
#'
#' @param counts A gene by droplet matrix
#' @param Z A droplet by cluster matrix of the posterior probability a 
#'  droplet belongs to cluster k. Rows must sum to 1 and entries must 
#'  be between 0 and 1.
#' @param test_set Character vector of droplet IDs that belong to the 
#'  test set
#' @param bg_set Character vector of droplet IDs that belong to the 
#'  background set
#' @param eps The threshold of change in the parameter to stop the 
#'  LOO iterations.
#' @param max_iter_loo The maximum number of LOO iterations
#' @param psc The pseudocount to add to the final MLE estimate. This avoids 
#'  inappropriate alpha values of 0.
#' @param ignore The posterior probability must be larger than this value 
#'  for it to count towards alpha estimation in a cluster.
#' @param threads Number of threads for parallel execution. Default is 1.
#' @param verbose Verbosity.
#' 
#' @importFrom Matrix rowSums t
#' @importFrom stats optimize
#
get_alpha <- function(counts, 
                      Z, 
                      test_set, 
                      bg_set, 
                      eps = 1e-4, 
                      max_iter_loo = 1e4, 
                      psc = 1e-10, 
                      ignore = 1e-10, 
                      threads = 1, 
                      verbose = FALSE){
    if (verbose) message("estimating alpha")

    clusts <- 1:ncol(Z)

    if (nrow(Z) != length(test_set)){
        stop("Length of test set is not equal to number of rows in Z")
    }

    # Get droplet membership
    Z_all <- matrix(0, nrow = ncol(counts), ncol = length(clusts))
    Z_all[,1] <- 1
    rownames(Z_all) <- colnames(counts)
    colnames(Z_all) <- clusts
    Z_all[test_set,] <- Z[test_set,]

    clust_mem <- Z_all > ignore

    alpha0 <- alpha_mom(counts = counts, 
                        clust_mem = clust_mem, 
                        weights = Z_all,  
                        psc = 0, 
                        threads = threads)

    Alpha <- alpha_max_loo(counts = counts, 
                           alpha0 = alpha0, 
                           clust_mem = clust_mem, 
                           weights = Z_all, 
                           eps = eps, 
                           psc = psc, 
                           max_iter = max_iter_loo, 
                           threads = threads)

    rownames(Alpha) <- rownames(counts)
    rm(Z_all)

    return(Alpha)
}

#' Get mixing coefficients of Mixture model from posterior probabilities Z
#'
#' @param Z A droplet by cluster matrix of the posterior probability a 
#'  droplet belongs to cluster k. Rows must sum to 1 and entries must 
#'  be between 0 and 1.
#' @param add A number, or weight, to add to the column sums of the 
#'  posterior probabilities. This number is added to the \code{add_to} index.
#' @param add_to The index to add the number \code{add} to.
#' 
#' @return A cluster-length numeric vector that sums to 1 and entries 
#'  between 0 and 1.
get_pi <- function(Z, add = 0, add_to = 1){
    Pi <- colSums(Z)
    Pi[add_to] <- Pi[add_to] + add
    Pi <- Pi / sum(Pi)
    return(Pi)
}

#' Get posterior probabilities from log likelihoods and mixing coefficients
#'
#' @param llk A droplet by cluster matrix that gives the log liklihood of 
#'  each droplet against each cluster.
#' @param Pi A cluster-length numeric vector of mixing coefficients.
#'
#' @return A droplet by cluster matrix of posterior probabilities
#'
get_z <- function(llk, Pi){
    if (length(Pi) != ncol(llk)){
        stop("Length of Pi must match the number of columns in llk")
    }

    K <- length(Pi)
    if (K == 1){
        Z <- matrix(1, nrow = K, ncol = 1)
        rownames(Z) <- rownames(llk)
        return(Z)
    }

    llk_pi <- t(apply(llk, 1, function(j) j + log(Pi)))
    Z <- as.matrix(t(apply(llk_pi, 1, fraction_log)))
    rownames(Z) <- rownames(llk)
    return(Z)
}

#' Run EM
#' 
#' @param counts A gene by droplet matrix
#' @param params A list containing
#' \describe{
#'      \item{Alpha}{ A gene by cluster matrix with parameters of 
#'          Dirichlet-multinomial.}
#'      \item{Pi}{ A vector of mixing coefficients. Must sum to 1 
#'          and be between 0 and 1.}
#'      }
#'  while 1 means background.
#' @param llk A droplet (row) by cluster (column) matrix of log-likelihoods 
#'  of a droplet belonging to each cluster. The log-likelihoods are only 
#'  calculated for the test set droplets and should correspond to the 
#'  likelihoods under the current values of the parameters given in 
#'  \code{params}.
#' @param test_set Character vector of droplet IDs that belong to the 
#'  test set
#' @param bg_set Character vector of droplet IDs that belong to the 
#'  background set
#' @param eps The threshold of change in the parameter to stop the 
#'  EM. The parameter checked is the average change in posterior 
#'  probabilities Z.
#' @param max_iter Maximum number of iterations.
#' @param threads Number of threads for parallel execution. Default is 1.
#' @param verbose Verbosity.
#'
#' @importFrom Matrix t colSums
em <- function(counts, 
               params, 
               llk, 
               test_set, 
               bg_set, 
               eps = 1e-4, 
               max_iter = 200, 
               threads = 1, 
               verbose = TRUE){
    sizes <- colSums(counts)

    Alpha <- params$Alpha
    Pi <- params$Pi
    Z <- get_z(llk, Pi)

    K <- ncol(Alpha)

    delta <- Inf
    iter <- 1
    while (iter <= max_iter & delta > eps){
        # Estimate pi
        Pi <- get_pi(Z, add = length(bg_set), add_to = 1)

        # Estimate alpha
        Alpha <- get_alpha(counts, 
                           Z, 
                           test_set = test_set, 
                           bg_set = bg_set)

        # Evaluate llk and estimate Z
        llk <- get_llk(counts = counts[,test_set], 
                       Alpha = Alpha, 
                       sizes = sizes[test_set], 
                       threads = threads)
        rownames(llk) <- test_set
        Z <- get_z(llk, Pi)

        # Evaluate delta
        if (iter > 1) {
            delta <- sum(abs(Z - Z_old)) / sum(Z_old)
        }
        if (verbose) message("iteration ", iter, "; delta = ", round(delta, 10))

        # Update
        iter <- iter + 1
        Alpha_old <- Alpha
        Z_old <- Z
        Pi_old <- Pi
    }

    converged <- FALSE
    if (iter > max_iter) {
        warning("warning: failed to converge after ", max_iter, " iterations")
    } else {
        converged <- TRUE
        if (verbose)
            message("converged after ", iter-1, " iterations")
    }

    Z <- get_z(llk, Pi)
    clust_max <- apply(Z, 1, which.max)
    names(clust_max) <- rownames(Z)

    params <- list("Alpha" = Alpha, 
                   "Pi" = Pi)
    ret <- list("params" = params, 
                "llk" = llk, 
                "converged" = converged, 
                "cluster" = clust_max)
    return(ret)
}

#' Run EM
#' 
#' Estimate the parameters of the Dirichlet-multinomial mixture model, 
#' and estimate the posterior probability a droplet belongs to each 
#' of the clusters.
#' 
#' @param x An SCE object.
#' @param eps The delta threshold for when to call convergence for 
#'  the EM estimation of the Dirichlet-multinomial mixture model. The EM 
#'  stops when delta falls below this value. We define delta as the 
#'  average change in posterior probability. By default this is set to 
#'  1e4, so that the EM converges when less than 1 in 10,000 labels 
#'  change on average.
#' @param max_iter_dm Maximum number of iterations for the EM estimation 
#'  of the Dirichlet-multinomial mixture model.
#' @param k_init Run EM on the \code{k_init} initialization(s). 
#'  If NULL (default), run on all \code{k_init} initializations.
#' @param threads Number of threads for parallel execution. Default is 1.
#' @param verbose Verbosity.
#'
#' @return An SCE object.
#'
#' @importFrom Matrix t colSums
#' @export
run_em <- function(x, 
                   eps = 1e-4, 
                   max_iter_dm = 200, 
                   k_init = NULL, 
                   threads = 1, 
                   verbose = TRUE){

    genes.use <- rownames(x@gene_data)[x@gene_data$exprsd]
    droplets.use <- colnames(x@counts)

    if (length(genes.use) == 0 | length(droplets.use) == 0) stop("specify test set and filter genes before running EM.")
    if (length(x@kruns) == 0) stop("initialize parameters before running run_em")

    # counts is a droplet by gene matrix
    sizes <- colSums(x@counts[genes.use,droplets.use])

    N <- length(droplets.use)
    G <- length(genes.use)

    # Run EM for each K start value
    k_init <- check_k_init(x, k_init)

    emo <- list()
    for (k in k_init){
        k <- as.character(k)
        if (verbose) message("running EM for k_init = ", k)
        x@kruns[[k]] <- em(counts = x@counts[genes.use,droplets.use],
                           params = x@kruns[[k]]$params, 
                           llk = x@kruns[[k]]$llk, 
                           test_set = x@test_set, 
                           bg_set = x@bg_set, 
                           eps = eps, 
                           max_iter = max_iter_dm, 
                           threads = threads, 
                           verbose = verbose)
        nclusters <- ncol(x@kruns[[k]]$llk)
        message("using a final k of 1 background and ", 
                nclusters, " cell type clusters for k_init ", k_init)
    }
    if (verbose){
        message("finished EM")
    }
    return(x)
}

