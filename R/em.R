
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
get_llk <- function(counts, Alpha, droplets = NULL, sizes = NULL, threads = 1, verbose = TRUE){

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
#' @param max_loo The maximum number of LOO iterations
#' @param psc The pseudocount to add to the final MLE estimate. This avoids 
#'  inappropriate alpha values of 0.
#' @param threads Number of threads for parallel execution. Default is 1.
#' @param tol The \code{tol} parameter for \code{\link[stats]{optimize}} 
#' used in the  initialization.
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
                      max_loo = 500, 
                      psc = 1e-10, 
                      threads = 1, 
                      tol = 1e2, 
                      verbose = TRUE){
    K <- ncol(Z)
    ks <- 1:K
    sizes <- colSums(counts)

    if (verbose) message("estimating alpha")

    clusts <- 1:ncol(Z)

    if (nrow(Z) != length(test_set)){
        stop("Length of test set is not equal to number of rows in Z")
    }

    # Get droplet membership
    clust_mem <- matrix(FALSE, nrow = ncol(counts), ncol = length(clusts))
    rownames(clust_mem) <- colnames(counts)
    colnames(clust_mem) <- clusts
    clust_mem[bg_set, 1] <- TRUE
    clust_mem[test_set,] <- apply(Z, 2, function(i) i > 0.01)

    zl <- lapply(1:K, function(k){
                 v1 <- rep(1, ncol(counts))
                 names(v1) <- colnames(counts)
                 v1[test_set] <- Z[,k]
                 return(v1[clust_mem[,k]])
                      })
    names(zl) <- clusts

    # Loop over droplets with high enough p(z)
    p_bar <- sapply(clusts, function(k){
                      p <- counts[,clust_mem[,k],drop=FALSE] %*% zl[[k]]
                      p <- as.matrix(p) + psc
                      p / sum(p) })

    # Initialize
    a_range <- c(10, 1e6)
    a0 <- c()
    a0 <- sapply(clusts, function(k){
                   f <- function(x){
                       r <- get_llk(counts = counts[,clust_mem[,k],drop=FALSE], 
                                    Alpha = x * p_bar[,k,drop=FALSE], 
                                    sizes = sizes[clust_mem[,k]], 
                                    threads = threads, 
                                    verbose = FALSE)
                       r <- as.numeric(t(r) %*% zl[[k]])
                       return(r)
                   }
                   a_max <- optimize(f, interval = a_range, tol = tol, maximum = TRUE)
                   return(a_max$maximum)
                      })

    # LOO optimize
    Alpha <- sapply(ks, function(k){
                      Ak <- compute_LOO_step_all(x = t(counts[,clust_mem[,k],drop=FALSE]), 
                                                 sizes = sizes[clust_mem[,k]], 
                                                 weights = zl[[k]], 
                                                 alpha = p_bar[,k] * a0[k], 
                                                 eps = eps, 
                                                 psc = psc, 
                                                 threads = threads, 
                                                 max_loo = max_loo)
                      return(Ak)})

    rownames(Alpha) <- rownames(counts)

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
               max_iter = 1e2, 
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

    params <- list("Alpha" = Alpha, 
                   "Pi" = Pi)
    ret <- list("params" = params, 
                "llk" = llk, 
                "converged" = converged)
    return(ret)
}

#' Run EM
#' 
#' Estimate the parameters of the Dirichlet-multinomial mixture model, filter 
#' out clusters close to the background distribution, and estimate the 
#' posterior probability a droplet belongs to each of the clusters. The 
#' \code{fltr} parameter controls the distance threshold to remove 
#' clusters.
#' 
#' @param x An SCE object.
#' @param eps The delta threshold for when to call convergence for 
#'  the EM estimation of the Dirichlet-multinomial mixture model. The EM 
#'  stops when delta falls below this value. We define delta as the 
#'  average change in posterior probability. By default this is set to 
#'  1e4, so that the EM converges when less than 1 in 10,000 labels 
#'  change on average.
#' @param fltr The filter threshold between 0 and 1 
#'  that controls the minimum distance to 
#'  the background distribution that a cluster can have. Remove  
#'  centers with a distance less than this value.
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
                   fltr = 0.1, 
                   max_iter_dm = 1e2, 
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
        while (TRUE){
            x@kruns[[k]] <- em(counts = x@counts[genes.use,droplets.use],
                               params = x@kruns[[k]]$params, 
                               llk = x@kruns[[k]]$llk, 
                               test_set = x@test_set, 
                               bg_set = x@bg_set, 
                               eps = eps, 
                               max_iter = max_iter_dm, 
                               threads = threads, 
                               verbose = verbose)
            prev_k <- length(x@kruns[[k]]$params$Pi)
            x <- get_dist(x, verbose = verbose)
            x <- rm_close(x, 
                          k_init = k, 
                          fltr = fltr, 
                          verbose = verbose)
            merged_k <- length(x@kruns[[k]]$params$Pi)
            if (merged_k == prev_k){
                break
            } else {
                if (verbose){
                    message("restarting EM")
                    nclusters <- merged_k - 1
                    message("using a final k of 1 background and ", 
                            nclusters, " cell type clusters")
                }
            }
        }
        if (verbose){
            message("finished EM")
        }
    }

    return(x)
}

