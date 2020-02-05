
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

get_z <- function(llks, Pi, labs){
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
    Z[labs == 1,1] <- 1
    Z[labs == 1,2:K] <- 0
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

rm_debris_clust <- function(counts, Alpha, Z, labs, Pi=NULL, thresh = 0.10, verbose = TRUE){
    if (verbose) message("checking to remove clusters similar to debris...")
    llks <- get_llk(counts, Alpha)
    if (is.null(Pi)) {
        llksf <- t(apply(llks, 1, fraction_log))
        pis <- colSums(llksf)
        Pi <- pis / sum(pis)
    }
    Z <- get_z(llks, Pi, labs)
    s = c()
    for (k in 2:ncol(Alpha)){
        diffs = llks[,k] - llks[,1]
        diffs = diffs * Z[,k]
        diffs = sum(diffs) / sum(Z[,k])
        s[k-1] = diffs
    }
    s_p <- s / max(s)
    ckeep <- s_p >= thresh
    ckeep <- c(TRUE, ckeep)

    if (verbose){
        torm <- length(ckeep) - sum(ckeep)
        if (torm > 0){
            message("removed ", torm, 
                    ngettext(torm, " cluster", " clusters"),
                    " similar to debris")
        }
    }
    Alpha <- Alpha[,ckeep,drop=FALSE]
    Pi <- Pi[ckeep] / sum(Pi[ckeep])
    Z <- t(apply(Z[,ckeep], 1, function(i) i / sum(i)))
    return(list("Alpha" = Alpha, "Pi" = Pi, "s" = s_p, "Z" = Z))
}

#' @importFrom future.apply future_sapply
#' @importFrom future nbrOfWorkers
#' @export
em <- function(counts, 
               Alpha, 
               Pi,
               labs, 
               eps = 1, 
               max_iter = 1e2, 
               verbose = TRUE){

    countst <- t(counts)
    sizes <- colSums(countst)
    K <- ncol(Alpha)

    delta <- Inf
    Z_old <- -Inf
    iter <- 1
    while (iter <= max_iter & delta > eps){

        # Estimate prob of Z
        llks <- get_llk(countst, Alpha, sizes)
        rownames(llks) <- colnames(countst)
        Z <- get_z(llks, Pi, labs)

        # Estimate pi
        Pi <- colSums(Z)
        Pi <- Pi / sum(Pi)

        # Estimate alpha
        Alpha <- dm_loo(counts, Alpha, Z, eps = 1e-10)

        # Evaluate change in Z
        delta <- sum(abs(Z - Z_old))
        if (verbose) message("iteration ", iter, "; delta = ", round(delta, 3))
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
                   "Pi" = Pi, 
                   "llk" = llks, 
                   "Z" = Z,  
                   "Cluster" = clust_max, 
                   "ClusterProb" = clust_prob, 
                   "ellk" = ellk)
    return(params)
}

#' Run EM
#' Run EM for each initialized K, storing output in SCE object.
#' @export
run_em <- function(x, 
                   K = NULL, 
                   eps = 1, 
                   thresh = 0.1, 
                   max_iter = 1e2, 
                   seedn = 1, 
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
    if (is.null(K))
        K <- names(x@init)
    else {
        K <- as.character(K)
        if ( any(!K %in% names(x@init)) ){
            stop("Value for K ", K, " is not found in ",
                 "initialized K values: ", 
                 paste(K, collapse = " "))
        }
    }

    set.seed(seedn)
    emo <- list()
    for (k in K){
        k <- as.character(k)
        if (verbose) message("running EM for k = ", k)
        Alpha <- x@init[[k]]$Alpha[genes.use,]
        Pi <- x@init[[k]]$Pi
        while (TRUE){
            emo[[k]] <- em(counts, 
                           Alpha, 
                           Pi, 
                           labs = labs, 
                           eps = eps, 
                           max_iter = max_iter, 
                           verbose = verbose)
            # Remove distributions close to debris
            Alpha <- emo[[k]]$Alpha
            Pi <- emo[[k]]$Pi
            Z <- emo[[k]]$Z
            prev_k <- ncol(Alpha)
            ret <- rm_debris_clust(counts, 
                                   Alpha, 
                                   Z, 
                                   Pi, 
                                   thresh = thresh, 
                                   verbose = verbose)
            Alpha <- ret$Alpha
            Pi <- ret$Pi
            merged_k <- ncol(Alpha)
            if (merged_k == prev_k){
                break
            } else {
                if (verbose){
                    message("restarting EM")
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

    # Fill in droplet data
    #x@droplet_data[,"CleanProb"] <- 1 - Z[,1]
    #x@droplet_data$Cluster <- NULL
    #x@droplet_data$Cluster <- clust_max
    #x@droplet_data$ClusterProb <- clust_prob

#' Select k
#' @export
test_k <- function(x, K = 1:15, thresh = 0.1, pct_debris = 5, eps = 5){
    #xk <- x
    #top_n <- round(length(xk@test_set) * (pct_debris/100))
    #tc <- order(xk@droplet_data[xk@bg_set,"total_counts"], decreasing=T)[1:top_n]
    #l2 <- xk@bg_set[tc]
    #xk@bg_set <- setdiff(xk@bg_set, l2)

    x <- init(x, K = K, thresh = thresh)
    x <- run_em(x, eps = eps, thresh = thresh)

    ret <- c()
    #ret <- matrix(nrow = length(K), ncol = 3)
    #rownames(ret) <- as.character(K)
    #colnames(ret) <- c("Pass", "Fail", "PercentFail")
    for (k in K){
        kc <- as.character(k)
        x <- call_targets(x, K = kc)
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

#' LOO optimization
#' 
#' counts is a droplet by gene matrix
#' sizes is total droplet size
#' a is current value of alpha parameter, gene by k matrix
#' z is current membership
dm_loo <- function(counts, A, Z, eps = 1e-4, max_loo = 3, add_to = 1e-16, verbose = FALSE){
    N <- nrow(counts)
    G <- ncol(counts)
    A_new <- A
    A_old <- A_new
    K <- ncol(Z)
    ks <- 1:K
    sizes <- rowSums(counts)
    for (k in ks){
        delt <- Inf
        iter <- 1
        while (delt > eps && iter <= max_loo){
            A_old[,k] <- A_new[,k]
            A_new[,k] <- compute_LOO_step(counts, sizes, Z[,k], A_old[,k])
            A_new[A_new[,k] < 0,k] <- 0
            A_new[, k] <- A_new[, k] + add_to
            delt <- sum((A_new[,k] - A_old[,k])^2) / sum((A_old[,k])^2)
            if (verbose){
                message("iteration ", iter, "; delta; ", delt)
            }
            iter <- iter + 1
        }
    }
    return(A_new)
}

