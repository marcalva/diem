
#' Get log likelihood under Dirichlet-multinomial
#'
#' Get log likelihood from \code{counts} under a Dirichlet-multinomial 
#' with parameter \code{Alpha}.
#'
#' @param counts A gene by droplet data matrix of read counts.
#' @param Alpha A gene by cluster matrix of parameter values for the 
#'  Dirichlet-multinomial. Each entry must be greater than 0.
#' @param labs Optional vector of fixed labels. Must be the same 
#'  length as number of columns in counts. Set llk of labeled to 
#'  infinity.
#' @param sizes Optional vector of droplet (column) sizes.
#' 
#' @importFrom future.apply future_sapply
#' @importFrom future nbrOfWorkers
#' @importFrom Matrix colSums
#'
get_llk <- function(counts, Alpha, labs = NULL, sizes = NULL){
    sapply_f <- ifelse(nbrOfWorkers() == 1, sapply, future_sapply)

    if (any(is.na(Alpha))){
        stop("Alpha has NA value(s)")
    }
    if (any(Alpha <= 0)){
        stop("Alpha has value(s) <= 0")
    }

    if (is.null(sizes)){
        sizes <- colSums(counts)
    }
    if (is.null(labs)){
        llks <- sapply_f(1:ncol(Alpha), function(j){
                         LlkDirMultSparse(counts, sizes, alpha = Alpha[,j,drop=FALSE]) })
        llks <- as.matrix(llks)
        rownames(llks) <- colnames(counts)
    } else {
        if (length(labs) != ncol(counts)){
            stop("labs must be same length as columns of counts")
        }
        labeld <- labs != 0
        unlabeld <- ! labeld
        groups <- setdiff(unique(labs), 0)
        
        llks <- matrix(0, nrow = ncol(counts), ncol = ncol(Alpha))
        rownames(llks) <- colnames(counts)
        for (g in groups){
            ix <- which(labs == g)
            llks[ix + ( nrow(llks) * (g-1) )] <- Inf
        }
        llks[unlabeld,] <- sapply_f(1:ncol(Alpha), function(j) {
                                    LlkDirMultSparse(counts[,unlabeld], sizes[unlabeld], alpha = Alpha[,j,drop=FALSE]) })

    }
    return(llks)
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
#' @param eps The threshold of change in the parameter to stop the 
#'  LOO iterations.
#' @param max_loo The maximum number of LOO iterations
#' @param psc The pseudocount to add to the final MLE estimate. This avoids 
#'  inappropriate alpha values of 0.
#' @param tol The \code{tol} parameter for \code{\link[stats]{optimize}} 
#' used in the  initialization.
#' 
#' @importFrom Matrix rowSums t
#' @importFrom stats optimize
#' @importFrom future.apply future_sapply
#' @importFrom future nbrOfWorkers
#
get_alpha <- function(counts, 
                      Z, 
                      eps = 1e-4, 
                      max_loo = 500, 
                      psc = 1e-10, 
                      tol = 1e2){
    sapply_f <- ifelse(nbrOfWorkers() == 1, sapply, future_sapply)
    countst <- t(counts)
    K <- ncol(Z)
    ks <- 1:K
    sizes <- colSums(counts)

    clusts <- 1:ncol(Z)

    # Loop over droplets with high enough p(z)
    clust_mem <- apply(Z, 2, function(i) i > 0.01)
    p_bar <- sapply_f(clusts, function(k){
                      p <- counts[,clust_mem[,k],drop=FALSE] %*% Z[clust_mem[,k],k]
                      p <- as.matrix(p) + psc
                      p / sum(p) })

    # Initialize
    a_range <- c(10, 1e6)
    a0 <- c()
    a0 <- sapply_f(clusts, function(k){
                   ct <- counts[,clust_mem[,k],drop=FALSE]
                   si <- sizes[clust_mem[,k]]
                   Zk <- Z[clust_mem[,k],k]
                   f <- function(x){
                       r <- get_llk(ct, 
                                    x * p_bar[,k,drop=FALSE], 
                                    sizes = si)
                       r <- as.numeric(t(r) %*% Zk)
                       return(r)
                   }
                   a_max <- optimize(f, interval = a_range, tol = tol, maximum = TRUE)
                   rm(Zk)
                   return(a_max$maximum)
                      })

    # LOO optimize
    Alpha <- sapply_f(ks, function(k){
                      ct <- countst[clust_mem[,k],,drop=FALSE]
                      si <- sizes[clust_mem[,k]]
                      Zk <- Z[clust_mem[,k],k]
                      Ak <- p_bar[,k] * a0[k]
                      Ak <- compute_LOO_step_all(ct, si, Zk, Ak, eps = eps, psc = psc, max_loo = max_loo)
                      rm(Zk)
                      return(Ak)})

    rownames(Alpha) <- rownames(counts)

    return(Alpha)
}

#' Get mixing coefficients of Mixture model from posterior probabilities Z
#'
#' @param Z A droplet by cluster matrix of the posterior probability a 
#'  droplet belongs to cluster k. Rows must sum to 1 and entries must 
#'  be between 0 and 1.
#' 
#' @return A cluster-length numeric vector that sums to 1 and entries 
#'  between 0 and 1.
get_pi <- function(Z){
    Pi <- colSums(Z)
    Pi <- Pi / sum(Pi)
    return(Pi)
}

#' Get posterior probabilities from log likelihoods and mixing coefficients
#'
#' @param llks A droplet by cluster matrix that gives the log liklihood of 
#'  each droplet against each cluster.
#' @param Pi A cluster-length numeric vector of mixing coefficients.
#'
#' @return A droplet by cluster matrix of posterior probabilities
#'
#' @importFrom future.apply future_apply
#' @importFrom future nbrOfWorkers
get_z <- function(llks, Pi){
    if (length(Pi) != ncol(llks)){
        stop("Length of Pi must match the number of columns in llks")
    }

    K <- length(Pi)
    if (K == 1){
        Z <- matrix(1, nrow = K, ncol = 1)
        rownames(Z) <- rownames(llks)
        return(Z)
    }

    apply_f <- ifelse(nbrOfWorkers() == 1, apply, future_apply)

    llks_pi <- t(apply_f(llks, 1, function(j) j + log(Pi)))
    Z <- as.matrix(t(apply_f(llks_pi, 1, fraction_log)))
    rownames(Z) <- rownames(llks)
    return(Z)
}

#' Fix posterior probabilities of labeled to 1
#'
#' @param Z A droplet by cluster matrix of the posterior probability a 
#'  droplet belongs to cluster k. Rows must sum to 1 and entries must 
#'  be between 0 and 1.
#' @param labs A droplet-length vector containing labels of the droplets. 
#'  Non-zero values indicate unlabeled droplets.
#' 
#' @return The droplet by cluster matrix with labeled droplets fixed 
#'  to ther corresponding cluster.
fix_z <- function(Z, labs = NULL){
    if (is.null(labs)) return(Z)
    if (sum(labs) == 0) return(Z)
    if (length(labs) != nrow(Z)){
        stop("length of labs must equal to number of rows in Z")
    }

    if (max(labs) > ncol(Z)){
        stop("Max label in labs must be less than or equal to column in Z")
    }

    f <- setdiff(sort(unique(labs)), 0)
    for (i in f){
        Z[labs == i,] <- 0
        Z[labs == i, i] <- 1
    }
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
#' @param labs Vector of labels for droplets. A label of 0 means unknwon 
#'  while 1 means background.
#' @param eps The threshold of change in the parameter to stop the 
#'  EM. The parameter checked is the average change in posterior 
#'  probabilities Z.
#' @param max_iter Maximum number of iterations.
#' @param verbose Verbosity.
#'
#' @importFrom Matrix t colSums
em <- function(counts, 
               params, 
               labs, 
               eps = 1e-4, 
               max_iter = 1e2, 
               verbose = TRUE){
    sizes <- colSums(counts)

    Alpha <- params$Alpha
    Pi <- params$Pi
    unl <- labs == 0

    K <- ncol(Alpha)

    delta <- Inf
    iter <- 1
    while (iter <= max_iter & delta > eps){

        # Estimate prob of Z
        llks <- get_llk(counts, Alpha, labs = labs, sizes)
        rownames(llks) <- colnames(counts)
        Z <- get_z(llks, Pi)

        # Estimate pi
        Pi <- get_pi(Z)

        # Estimate alpha
        Alpha <- get_alpha(counts, Z)

        # Evaluate delta
        if (iter > 1) {
            delta <- sum(abs(Z[unl,] - Z_old[unl,])) / sum(Z_old[unl,])
        }
        if (verbose) message("iteration ", iter, "; delta = ", round(delta, 10))
        iter <- iter + 1
        Alpha_old <- Alpha
        Z_old <- Z
        Pi_old <- Pi
    }

    converged <- FALSE
    if (iter > max_iter) {
        warning("Warning: failed to converge after ", max_iter, " iterations")
    } else {
        converged <- TRUE
        if (verbose)
            message("converged after ", iter-1, " iterations")
    }

    params <- list("Alpha" = Alpha, 
                   "Pi" = Pi)
    ret <- list("params" = params, 
                "llk" = llks, 
                "Z" = Z,  
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
                   verbose = TRUE){

    genes.use <- rownames(x@gene_data)[x@gene_data$exprsd]
    droplets.use <- colnames(x@counts)

    if (length(genes.use) == 0 | length(droplets.use) == 0) stop("Specify test set and filter genes before running EM.")

    # counts is a droplet by gene matrix
    sizes <- colSums(x@counts[genes.use,droplets.use])

    N <- length(droplets.use)
    G <- length(genes.use)

    # Fix labels of debris
    labs <- rep(0, N)
    names(labs) <- droplets.use
    labs[x@bg_set] <- 1

    # Run EM for each K start value
    k_init <- check_k_init(x, k_init)

    emo <- list()
    for (k in k_init){
        k <- as.character(k)
        if (verbose) message("running EM for k_init = ", k)
        while (TRUE){
            x@kruns[[k]] <- em(x@counts[genes.use,droplets.use],
                               x@kruns[[k]]$params, 
                               labs = labs, 
                               eps = eps, 
                               max_iter = max_iter_dm, 
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

