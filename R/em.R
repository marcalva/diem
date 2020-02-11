
#' fraction of logs
#'
#' @param x numeric vector
#' @export
fraction_log <- function(x){
    x_c = x - max(x)
    x_c = exp(x_c)
    frac = x_c / sum(x_c);
    return(frac);
}

#' sum of logs
#'
#' @param x numeric vector
#' @export
sum_log <- function(x){
    max_x <- max(x)
    x_c = x - max_x
    x_sum <- log(sum(exp(x_c))) + max_x
    return(x_sum)
}

#' Get total log likelihood
#' @importFrom future.apply future_sapply
#' @importFrom future nbrOfWorkers
#
# @export
get_llk <- function(counts, Alpha, sizes = NULL){
    sapply_f <- ifelse(nbrOfWorkers() == 1, sapply, future_sapply)

    if (is.null(sizes)){
        sizes <- colSums(counts)
    }

    llks <- sapply_f(1:ncol(Alpha), function(j) LlkDirMultSparse(counts, sizes, alpha = Alpha[,j,drop=FALSE]) )
    llks <- matrix(llks, ncol = ncol(Alpha))
    rownames(llks) <- colnames(counts)
    return(llks)
}

#' LOO optimization
#' 
#' counts is a droplet by gene matrix
#' sizes is total droplet size
#' a is current value of alpha parameter, gene by k matrix
#' z is current membership
#' @importFrom future.apply future_sapply
#' @importFrom future nbrOfWorkers
#
# @export
get_alpha <- function(counts, 
                      Z, 
                      eps = 1e-4, 
                      max_loo = 500, 
                      psc = 1e-10, 
                      debug = FALSE){
    sapply_f <- ifelse(nbrOfWorkers() == 1, sapply, future_sapply)
    countst <- t(counts)
    K <- ncol(Z)
    ks <- 1:K
    sizes <- rowSums(counts)

    clusts <- 1:ncol(Z)

    # Loop over droplets with high enough p(z)
    clust_mem <- apply(Z, 2, function(i) i > 0.01)
    p_bar <- sapply_f(clusts, function(k){
                      p <- countst[,clust_mem[,k],drop=FALSE] %*% Z[clust_mem[,k],k]
                      p <- as.matrix(p) + psc
                      p / sum(p) })

    a_range <- c(10, 1e6)
    a0 <- c()
    a0 <- sapply_f(clusts, function(k){
                   ct <- countst[,clust_mem[,k],drop=FALSE]
                   si <- sizes[clust_mem[,k]]
                   Zk <- Z[clust_mem[,k],k]
                   f <- function(x){
                       r <- get_llk(ct, 
                                    x * p_bar[,k,drop=FALSE], 
                                    sizes = si)
                       r <- as.numeric(t(r) %*% Zk)
                       return(r)
                   }
                   a_max <- optimize(f, interval = a_range, tol = 100, maximum = TRUE)
                   a_max$maximum
                      })

    Alpha <- sapply_f(ks, function(k){
                      ct <- counts[clust_mem[,k],,drop=FALSE]
                      si <- sizes[clust_mem[,k]]
                      Zk <- Z[clust_mem[,k],k]
                      Ak <- p_bar[,k] * a0[k]
                      Ak_old <- Ak
                      delt <- Inf
                      iter <- 1
                      while (delt > eps && iter <= max_loo){
                          Ak <- compute_LOO_step(ct, si, Zk, Ak)
                          Ak[Ak < 0] <- 0
                          Ak <- Ak + psc
                          delt <- sum(abs(Ak - Ak_old)) / sum(Ak_old)
                          if (debug){
                              message("iteration ", iter, "; sum AK = ", sum(Ak), "; delta = ", delt)
                          }
                          Ak_old <- Ak
                          iter <- iter + 1
                      }
                      return(Ak)
                      })

    rownames(Alpha) <- colnames(counts)

    return(Alpha)
}

get_pi <- function(Z){
    Pi <- colSums(Z)
    Pi <- Pi / sum(Pi)
    return(Pi)
}

#' Get Z
#'
#' @importFrom future.apply future_apply
#' @importFrom future nbrOfWorkers
#
# @export
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
    Z <- t(apply_f(llks_pi, 1, fraction_log))
    rownames(Z) <- rownames(llks)
    return(Z)
}

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

#' Simulate from dirichlet-multinoimal
rdirm <- function(n, size, a){
    if (length(size) == 1) size <- rep(size, n)
    nt <- length(a) * n
    xi <- rgamma(n = nt, shape = a, scale = 1)
    datf <- matrix(xi, nrow = length(a))
    datf <- apply(datf, 2, function(i) i / sum(i))
    ret <- sapply(1:n, function(i) rmultinom(1, size = size[i], prob = datf[,i]))
    return(ret)
}

#' Get cluster distances to background disrtibution
#' 
#' Get the likelihood-based distance of the clusters to the background.
#'
#' @params x An SCE object.
#' @params llks An optional droplet by cluster matrix containing the 
#'  log likelihoods of the droplet given the cluster's parameters. If 
#'  give, avoids re-calcluating the log likelihoods.
#' @params verbose Verbosity.
#'
#' @return An SCE object
get_dist <- function(x, llks = NULL, verbose = TRUE){
    if (verbose) message("checking distances to background distribution...")

    genes.use <- rownames(x@gene_data)[x@gene_data$exprsd]
    counts <- x@counts[genes.use,]

    k_init <- names(x@init)

    for (k in k_init){
        kc <- as.character(k)

        i <- length(x@init[[kc]])
        ic <- x@init[[kc]][[i]]
        params <- ic$params
        Alpha <- params$Alpha
        Pi <- params$Pi
        Z <- ic$Z

        if (ncol(Alpha) == 1){
            if (verbose){
                message("only the debris distribution is present, no clusters to remove")
            }
            return(params)
        }

        if (is.null(llks)){
            llks <- get_llk(counts, Alpha)
        }

        d <- c(0)
        for (k in 2:ncol(Alpha)){
            diffs <- llks[,k] - llks[,1]
            diffs <- diffs * Z[,k]
            wsum <- sum(Z[,k])
            if (wsum == 0){
                d[k] <- 0
            } else {
                diffs <- sum(diffs) / wsum
                d[k] <- diffs
            }
        }
        d[d < 0] <- 0
        d[is.na(d)] <- 0
        if (max(d) == 0){
            warnings("Warning: all distances are 0")
            return(params)
        }
        ds <- d / max(d)
        x@init[[kc]][[i]]$Dist <- ds
    }
    return(x)

}

#' @importFrom future.apply future_sapply
#' @importFrom future nbrOfWorkers
#' @export
em <- function(counts, 
               params, 
               labs, 
               eps = 1e-4, 
               max_iter = 1e2, 
               verbose = TRUE){

    countst <- t(counts)
    sizes <- colSums(countst)

    Alpha <- params$Alpha
    Pi <- params$Pi
    unl <- labs == 0

    K <- ncol(Alpha)

    delta <- Inf
    iter <- 1
    while (iter <= max_iter & delta > eps){

        # Estimate prob of Z
        llks <- get_llk(countst, Alpha, sizes)
        rownames(llks) <- colnames(countst)
        Z <- get_z(llks, Pi)
        Z <- fix_z(Z, labs)

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

    if (iter > max_iter) {
        warning("Warning: failed to converge after ", max_iter, " iterations")
    } else {
        if (verbose)
            message("converged after ", iter-1, " iterations")
    }

    clust_max <- apply(Z, 1, which.max)
    clust_prob <- apply(Z, 1, function(i) i[which.max(i)])
    ellk <- sum(sapply(1:nrow(Z), function(i) llks[i, clust_max[i]]))
    params <- list("Alpha" = Alpha, 
                   "Pi" = Pi)
    ret <- list("params" = params, 
                "llk" = llks, 
                "Z" = Z,  
                "Cluster" = clust_max, 
                "ClusterProb" = clust_prob, 
                "ellk" = ellk)
    return(ret)
}

#' Run EM
#' Run EM for each initialized K, storing output in SCE object.
#' @export
run_em <- function(x, 
                   k_init = NULL, 
                   eps = 1e-4, 
                   fltr = 0.1, 
                   max_iter = 1e2, 
                   verbose = TRUE){

    genes.use <- rownames(x@gene_data)[x@gene_data$exprsd]
    droplets.use <- colnames(x@counts)

    if (length(genes.use) == 0 | length(droplets.use) == 0) stop("Specify test set and filter genes before running EM.")

    # counts is a droplet by gene matrix
    countst <- x@counts[genes.use,droplets.use]
    counts <- t(x@counts[genes.use,droplets.use])
    sizes <- colSums(countst)

    N <- nrow(counts)
    G <- ncol(counts)

    # Fix labels of debris
    labs <- rep(0, N)
    names(labs) <- droplets.use
    labs[x@bg_set] <- 1

    # Run EM for each K start value
    if (is.null(k_init))
        k_init <- names(x@init)
    else {
        k_init <- as.character(k_init)
        if ( any(!k_init %in% names(x@init)) ){
            stop("Value for k_init ", k_init, " is not found in ",
                 "initialized k_init values: ", 
                 paste(k_init, collapse = " "))
        }
    }

    emo <- list()
    for (k in k_init){
        k <- as.character(k)
        if (verbose) message("running EM for k_init = ", k)
        initr <- x@init[[k]][[length(x@init[[k]])]]
        params <- initr$params
        while (TRUE){
            emo[[k]] <- em(counts, 
                           params, 
                           labs = labs, 
                           eps = eps, 
                           max_iter = max_iter, 
                           verbose = verbose)
            # Remove distributions close to debris
            prev_k <- ncol(emo[[k]]$params$Alpha)
            x <- get_dist(x, llks = emo[[k]]$llk, verbose = verbose)
            x <- rm_close(x, fltr = fltr, verbose = verbose)
            initr <- x@init[[k]][[length(x@init[[k]])]]
            params <- initr$params
            emo[[k]][["Dist"]] <- params$Z
            emo[[k]]$params <- params
            merged_k <- ncol(emo[[k]]$params$Alpha)
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
    }

    x@emo <- emo

    if (verbose){
        message("finished EM")
    }

    return(x)
}

#' Select k
#' @export
test_k <- function(x, k_init = 1:30, fltr = 0.1, pct_debris = 5){
    #xk <- x
    #top_n <- round(length(xk@test_set) * (pct_debris/100))
    #tc <- order(xk@droplet_data[xk@bg_set,"total_counts"], decreasing=T)[1:top_n]
    #l2 <- xk@bg_set[tc]
    #xk@bg_set <- setdiff(xk@bg_set, l2)

    x <- init(x, k_init = k_init)
    x <- run_em(x)

    ret <- c()
    #ret <- matrix(nrow = length(K), ncol = 3)
    #rownames(ret) <- as.character(K)
    #colnames(ret) <- c("Pass", "Fail", "PercentFail")
    for (k in k_init){
        kc <- as.character(k)
        x <- call_targets(x, k_init = kc)
        emo <- x@emo[[kc]]
        n_pass <- sum(x@droplet_data[,"Call"] == "Clean")
        ret[kc] <- n_pass
        #nf <- sum(emo$cluster[l2] != 1)
        #ret[kc,"Pass"] <- n_pass
        #ret[kc,"Fail"] <- nf
        #ret[kc,"PercentFail"] <- nf / top_n
    }
    x@test_k <- ret
    return(x)
}

