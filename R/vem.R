
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

#' Normalizing constant for Dirichlet
#' @export
dir_c <- function(x, log=TRUE){
    if (log){
        lgamma(sum(x)) - sum(lgamma(x))
    } else {
        gamma(sum(x)) / prod(gamma(x))
    }
}

#' Entropy of Dirichlet
#'
hd <- function(a){
    K <- length(a)
    a0 <- sum(a)
    (digamma(a0) * (a0 - K) )- (digamma(a) %*% (a-1)) - dir_c(a)
}

init_best <- function(X, 
                      labs = NULL, 
                      K, 
                      top_n = Inf,  
                      init_start = 5, 
                      init_iter = 3, 
                      verbose = TRUE){
    N <- nrow(X)
    G <- ncol(X)
    X_unl <- X[labs == 0,]
    N_unl <- nrow(X_unl)
    o <- order(rowSums(X_unl), decreasing = TRUE)
    top_n <- min(top_n, N_unl)
    keep <- o[1:top_n]
    pl <- list()
    for (i in 1:init_start){
        pinit <- init_random(X, N, G, K, labs = labs)
        params <- dm(X, Alpha = pinit$Alpha, labs = labs, max_iter = init_iter)
        pl[[i]] <- params
    }
    llks <- sapply(pl, function(i) i$ellk)
    i_max <- which.max(llks)
    params <- pl[[i]]$params
    return(params)
}

#' Initialize randomly
#'
#' @export
init_random <- function(X, 
                        N, 
                        G, 
                        K, 
                        labs = NULL, 
                        add_to = 1e-16){
    clusts <- 1:K
    if (length(labs) != N) labs <- rep(0, N)
    if (length(labs) != N){
        stop("Size of labs must be the same as the number of rows in x.")
    }

    Xt <- t(X)

    # Get labeled groups
    clusts <- 1:K
    labeled <- setdiff(unique(labs), 0)
    unlabeled <- setdiff(clusts, labeled)
    
    labs_init <- labs
    names(labs_init) <- seq(N)

    unl_ix <- labs == 0
    labs_init[unl_ix] <- sample(1:K, size = sum(unl_ix), replace = TRUE)

    Alpha <- sapply(clusts, function(k) Matrix::colSums(X[labs_init == k,]))
    Alpha <- Alpha + add_to
    Z <- model.matrix(~ 0 + as.factor(labs_init))
    colnames(Z) <- NULL
    rownames(Z) <- rownames(X)
    return(list("Alpha" = Alpha, "Z" = Z))
}

dm <- function(counts, Alpha, labs, eps = 1, max_iter = 1e2){
    countst <- t(counts)
    sizes <- colSums(countst)
    K <- ncol(Alpha)
    llks <- sapply(1:K, function(k) LlkDirMultSparse(countst, sizes = sizes, alpha = Alpha[,k]))
    Z <- t(apply(llks, 1, fraction_log))
    rownames(Z) <- rownames(counts)
    Z[labs == 1,1] <- 1
    Z[labs == 1,2:K] <- 0

    delta <- Inf
    iter <- 1
    while (iter <= max_iter & delta > eps){
        Alpha_new <- dm_loo(counts, Alpha, Z, eps = 1e-10)
        llks <- sapply(1:K, function(k) LlkDirMultSparse(countst, sizes = sizes, alpha = Alpha_new[,k]))
        Z_new <- t(apply(llks, 1, fraction_log))
        rownames(Z_new) <- rownames(counts)
        Z_new[labs == 1,1] <- 1
        Z_new[labs == 1,2:K] <- 0
        delta <- sum(abs(Z_new - Z))
        print(delta)
        print(colSums(Z_new))
        iter <- iter + 1
        Alpha <- Alpha_new
        Z <- Z_new
    }
    clust_max <- apply(Z, 1, which.max)
    clust_prob <- apply(Z, 1, function(i) i[which.max(i)])
    ellk <- sum(sapply(1:nrow(Z), function(i) llks[i, clust_max[i]]))
    params <- list("Alpha" = Alpha, 
                   "llk" = llks, 
                   "Z" = Z,  
                   "Cluster" = clust_max, 
                   "ClusterProb" = clust_prob, 
                   "ellk" = ellk)
    return(params)
}

#' Run dm
run_dm <- function(x, 
                   K = 30, 
                   eps = 1, 
                   max_iter = 1e2, 
                   Alpha = NULL, 
                   verbose = TRUE){

    genes.use <- rownames(x@gene_data)[x@gene_data$exprsd]
    droplets.use <- colnames(x@counts)

    if (length(genes.use) == 0 | length(droplets.use) == 0) stop("Specify test set and filter genes before running EM.")

    countst <- x@counts[genes.use,droplets.use]
    counts <- t(x@counts[genes.use,droplets.use])
    sizes <- colSums(countst)

    N <- nrow(counts)
    G <- ncol(counts)

    # Fix labels of debris
    labs <- rep(0, N)
    names(labs) <- droplets.use
    labs[x@bg_set] <- 1

    # Random init?
    # pinit <- init_random(counts, N, G, K)
    # Initialize alphas and Z
    # A_old <- pinit$Alpha
    if (is.null(Alpha)){
        A_old <- init_dm(counts, x@ic$labels)
    } else {
        A_old <- Alpha
    }
    params <- dm(counts, A_old, labs = labs, eps = eps, max_iter = max_iter)
    Alpha <- params$Alpha
    Z <- params$Z

    clust_max <- apply(Z, 1, which.max)
    clust_prob <- apply(Z, 1, function(i) i[which.max(i)])
    ellk <- sum(sapply(1:nrow(Z), function(i) llks[i, clust_max[i]]))

    emo <- list("Alpha" = A_new, "Z" = Z, "llk" = ellk)

    # Fill in droplet data
    x@droplet_data[,"CleanProb"] <- 1 - Z[,1]
    x@droplet_data$Cluster <- NULL
    x@droplet_data$Cluster <- clust_max
    x@droplet_data$ClusterProb <- clust_prob

    x@emo <- emo

    if (verbose){
        message("Finished EM")
    }

    return(x)
}

#' Select k
get_k <- function(x, kmin = 2, kmax = 5, top_n = 500){
    ret <- matrix(nrow = length(kmin:kmax), ncol = 3)
    rownames(ret) <- kmin:kmax
    colnames(ret) <- c("Pass", "Fail", "PercentFail")
    for (k in kmin:kmax){
        xk <- x
        tc <- order(xk@droplet_data[xk@bg_set,"total_counts"], decreasing=T)[1:top_n]
        l2 <- xk@bg_set[tc]
        xk@bg_set <- setdiff(xk@bg_set, l2)
        xk <- initialize_clusters(xk, K = k, n_start = 10, km_iter = 3)
        xk <- run_dm(xk, K = k, max_iter = 5)
        n_pass <- sum(xk@droplet_data$Cluster != 1)
        nf <- sum(xk@droplet_data[l2,"Cluster"] != 1)
        ret[k,"Pass"] <- n_pass
        ret[k,"Fail"] <- nf
        ret[k,"PercentFail"] <- nf / top_n
    }
    return(ret)
}

        

#' Estimate parameters of Dirichlet
#'
#' x are the observed data values of dirichlet
#' a are the alpha parameters of dirichlet
#' x and a must be same length
est_dir <- function(x, 
                    a, 
                    N = 1, 
                    rho = .1, 
                    eps = 1e-10, 
                    max_iter = 1e5){
    if (length(x) != length(a)) stop("x and a must be the same length")

    a_old <- a
    delt <- Inf
    iter <- 1
    while (delt > eps & iter <= max_iter){
        a_old_sum <- sum(a_old)
        g <- digamma(a_old_sum) - digamma(a_old) + x
        g <- N * g
        h <- -N * trigamma(a_old)
        z <- N * trigamma(a_old_sum)
        b <- (sum(g/h)) / ( (1/z) + sum(1/h) )

        Hinv_g <- (g - b) / h
        a_new <- a_old - rho*Hinv_g
        a_new[a_new < 0] <- 1e-16
        delt <- (a_new - a_old) %*% (a_new - a_old)
        a_old <- a_new
        iter <- iter + 1
    }
    p <- lgamma(sum(a_new)) - lgamma(sum(a_new)) + (a_new-1) %*% x
    print(iter)
    print(p)
    return(a_new)
}

#' Initialize Dirichlet-Multinomial parameters from k-means init
#'
init_dm <- function(counts, labs, add_to = 1e-16){
    clusts <- unique(labs)
    rs <- sapply(clusts, function(k) colSums(counts[labs == k,]))
    # rs <- sweep(rs, 2, colSums(rs), "/")
    rs <- rs + add_to
    return(rs)
}

#' Newton-Raphson optimization for Dir-Mult with sparse matrix
#' 
#' counts is a droplet by gene matrix
#' a is current value of alpha parameter, gene by k matrix
#' z is current membership
dm_nr <- function(counts, A, Z, eps = 1e-3, max_newton = 10, rho = .1, add_to = 1e-16, tol = 100){
    N <- nrow(counts)
    G <- ncol(counts)
    A_new <- A
    A_old <- A_new
    K <- ncol(Z)
    ks <- 1:K
    sizes <- rowSums(counts)
    for (k in ks){
        message("K = ", k)
        print(k)
        delt <- Inf
        iter <- 1
        while (delt > eps && iter <= max_newton){
            A_old[,k] <- A_new[,k]
            hg <- compute_step(counts, sizes, Z[,k], A_old[,k], tol = tol)
            names(hg) <- colnames(counts)
            A_new[,k] <- A_old[,k] - (rho*hg)
            psc <- min(add_to, min(A_new[A_new[,k] > 0,k]))
            A_new[A_new[,k] < 0, k] <- psc
            # delt <- sum((A_new[,k] - A_old[,k])^2) / sum((A_old[,k])^2)
            delt <- mean(abs((A_new[,k] - A_old[,k])/A_old[,k]))
            message("Delta: ", delt)
            iter <- iter + 1
        }
    }
    return(A_new)
}

#' LOO optimization
#' 
#' counts is a droplet by gene matrix
#' a is current value of alpha parameter, gene by k matrix
#' z is current membership
dm_loo <- function(counts, A, Z, eps = 1e-4, max_loo = 3, add_to = 1e-16){
    N <- nrow(counts)
    G <- ncol(counts)
    A_new <- A
    A_old <- A_new
    K <- ncol(Z)
    ks <- 1:K
    sizes <- rowSums(counts)
    for (k in ks){
        print(k)
        delt <- Inf
        iter <- 1
        while (delt > eps && iter <= max_loo){
            A_old[,k] <- A_new[,k]
            A_new[,k] <- compute_LOO_step(counts, sizes, Z[,k], A_old[,k])
            A_new[A_new[,k] < 0,k] <- 0
            A_new[, k] <- A_new[, k] + add_to
            delt <- sum((A_new[,k] - A_old[,k])^2) / sum((A_old[,k])^2)
            iter <- iter + 1
        }
    }
    return(A_new)
}

#' Fixed point optimization
#' 
#' counts is a droplet by gene matrix
#' a is current value of alpha parameter, gene by k matrix
#' z is current membership
dm_fp <- function(counts, A, Z, eps = 1e-4, max_fp = 3, add_to = 1e-10, tol = 100){
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
        while (delt > eps && iter <= max_fp){
            A_old[,k] <- A_new[,k]
            A_new[,k] <- compute_fp_step(counts, sizes, Z[,k], A_old[,k], tol = tol)
            # A_new[,k] <- A_new[,k] + add_to
            psc <- min(add_to, min(A_new[A_new[,k] > 0,k]))
            A_new[, k] <- A_new[, k] + psc
            delt <- mean(abs((A_new[,k] - A_old[,k])/A_old[,k]))
            # delt <- sum((A_new[,k] - A_old[,k])^2) / sum((A_old[,k])^2)
            print(delt)
            iter <- iter + 1
        }
    }
    return(A_new)
}

