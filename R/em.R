
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
    return(x_sum);
}

#' @importFrom future.apply future_sapply
#' @importFrom future nbrOfWorkers
#' @export
em <- function(counts, 
               Alpha, 
               labs, 
               eps = 1, 
               max_iter = 1e2, 
               verbose = TRUE){

    sapply_f <- ifelse(nbrOfWorkers() == 1, sapply, future_sapply)
    apply_f <- ifelse(nbrOfWorkers() == 1, apply, future_apply)
    
    countst <- t(counts)
    sizes <- colSums(countst)
    K <- ncol(Alpha)
    llks <- sapply_f(1:K, function(k) LlkDirMultSparse(countst, sizes = sizes, alpha = Alpha[,k]))
    Z <- t(apply_f(llks, 1, fraction_log))
    rownames(Z) <- rownames(counts)
    Z[labs == 1,1] <- 1
    Z[labs == 1,2:K] <- 0
    Pi <- colSums(Z)
    Pi <- Pi / sum(Pi)

    delta <- Inf
    iter <- 1
    while (iter <= max_iter & delta > eps){
        Alpha_new <- dm_loo(counts, Alpha, Z, eps = 1e-10)
        llks <- sapply_f(1:K, function(k) LlkDirMultSparse(countst, sizes = sizes, alpha = Alpha_new[,k]))
        llks_pi <- t(apply_f(llks, 1, function(j) j + log(Pi)))
        Z_new <- t(apply_f(llks_pi, 1, fraction_log))
        rownames(Z_new) <- rownames(counts)
        Z_new[labs == 1,1] <- 1
        Z_new[labs == 1,2:K] <- 0
        Pi_new <- colSums(Z_new)
        Pi_new <- Pi_new / sum(Pi_new)
        delta <- sum(abs(Z_new - Z))
        if (verbose) message("iteration ", iter, "; delta = ", round(delta, 3))
        iter <- iter + 1
        Alpha <- Alpha_new
        Z <- Z_new
        Pi <- Pi_new
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

    emo <- list()
    for (k in K){
        k <- as.character(k)
        if (verbose) message("running EM for k = ", k)
        emo[[k]] <- em(counts, 
                       x@init[[k]]$Alpha[genes.use,], 
                       labs = labs, 
                       eps = eps, 
                       max_iter = max_iter, 
                       verbose = verbose)
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
test_k <- function(x, K = 1:15, pct_debris = 5, eps = 5){
    xk <- x
    top_n <- round(length(xk@test_set) * (pct_debris/100))
    tc <- order(xk@droplet_data[xk@bg_set,"total_counts"], decreasing=T)[1:top_n]
    l2 <- xk@bg_set[tc]
    xk@bg_set <- setdiff(xk@bg_set, l2)

    xk <- init(xk, K = K)
    xk <- run_em(xk, eps = eps)

    ret <- matrix(nrow = length(K), ncol = 3)
    rownames(ret) <- as.character(K)
    colnames(ret) <- c("Pass", "Fail", "PercentFail")
    for (k in K){
        kc <- as.character(k)
        xk <- call_targets(xk, K = kc)
        emo <- xk@emo[[kc]]
        n_pass <- sum(xk@droplet_data[,"Call"] == "Clean")
        nf <- sum(emo$cluster[l2] != 1)
        ret[kc,"Pass"] <- n_pass
        ret[kc,"Fail"] <- nf
        ret[kc,"PercentFail"] <- nf / top_n
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

