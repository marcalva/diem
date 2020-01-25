
get_dist <- function(X, centers){
    d <- apply(X, 1, function(x) colSums((x - centers)^2))
    return(d)
}

get_wss <- function(Xt, centers, l){
    clusts <- 1:ncol(centers)
    wss <- sapply(clusts, function(i) sum(((Xt[,l == i]) - centers[,i])^2))
    return(wss)
}


#' Run k-means on multinomial
#' @export
kmeans_ss <- function(X, 
                      K = 30, 
                      labs = NULL, 
                      max_iter = 10, 
                      scale. = FALSE){
    N <- nrow(X)
    G <- ncol(X)

    if (scale.) X <- scale(X)

    Xt <- t(X)

    if (length(labs) != N) labs <- rep(0, N)

    # Initialize centers
    centers <- matrix(0, nrow = G, ncol = K)

    clusts <- 1:K
    labeled <- setdiff(unique(labs), 0)
    unlabeled <- setdiff(clusts, labeled)
    uix <- seq(N)[labs == 0]
    Xmu <- as.matrix(X[uix,])

    # Initialize labeled
    for (l in labeled){
        centers[,l] <- rowMeans(Xt[,labs == l])
    }
    # Initialize unlabeled
    uix_s <- sample(uix, size = length(unlabeled))
    centers[,unlabeled] <- as.matrix(Xt[,uix_s,drop=FALSE])

    p_l <- labs
    # First iter
    d <- get_dist(Xmu, centers)
    asgn <- apply(d, 2, which.min)
    p_l[uix] <- asgn

    # Wss
    wss <- get_wss(Xt, centers, p_l)

    iter <- 1
    while(iter <= max_iter){
        n_l <- p_l
        # Calculate centers
        centers <- sapply(clusts, function(i) rowMeans(Xt[, p_l == i, drop=FALSE]))

        # Calculate distances
        d <- get_dist(Xmu, centers)
        asgn <- apply(d, 2, which.min)
        if (length(intersect(unique(asgn), unlabeled)) != length(unlabeled)){
            break
        }
        n_l[uix] <- asgn
        wss <- lapply(clusts, function(i) sum((t(X[p_l == i,]) - centers[,i])^2))
        wss <- do.call(sum, wss)
        if (all(p_l == n_l)){
            break
        }
        p_l <- n_l
        iter <- iter + 1
    }

    return(list("labels" = factor(n_l, levels = 1:K), 
                "centers" = centers, 
                "withinss" = wss, 
                "tot.withinss" = sum(wss)))

}

#' Initialize clustering
#'
#' Instead of randomly initializing the EM, cell types are estimated from 
#' droplets that are expected to contain cells/nuclei. The initialization 
#' is done with droplets in the cluster set. The data is then normalized 
#' by first calculating the variable genes. A loess regression line 
#' is fit between the log counts and log variance, and the only top genes 
#' ranked by residual are used to initialize the clusters. The number of 
#' genes is specified with \code{n_var}. Optionally, one can use all genes 
#' by setting \code{use_var} to FALSE. The span of the loess regression line 
#' is given by \code{lss} (default is 0.3).
#' The data is normalized by dividing counts by the total counts per droplet. 
#' Then, the counts are multiplied by a scaling factor, given by 
#' \code{sf} (the median of total counts by default). Finally, the data is 
#' log transformed after adding a constant value of 1.
#' After normalization, the k-nearest neighbors are identified in the 
#' cluster set. The number of nearest neighbors is specified by 
#' \code{nn}. Clusters are identified from the KNN graph 
#' using the Louvain algorithm. Finally, only clusters with at least 
#' \code{min_size} (20 by default) droplets are considered cell types.
#'
#' @param x An SCE object.
#' @param K Number of clusters to return from k-means. This gives 
#'  the Maximum number of clusters
#' @param verbose verbosity.
#'
#' @return An SCE object
#' @export
#' @export
initialize_clusters <- function(x, 
                                K = 30, 
                                n_start = 10, 
                                km_iter = 3, 
                                verbose = FALSE){
    if (verbose) message("initializing clusters using k-means")
    labs <- rep(0, nrow(x@pcs))
    names(labs) <- rownames(x@pcs)
    labs[intersect(names(labs), x@bg_set)] <- 1
    pcss <- scale(x@pcs)
    ks <- replicate(n_start, 
                    expr = {
                        kmeans_ss(pcss, 
                                  labs = labs, 
                                  K = K, 
                                  max_iter = km_iter, 
                                  scale. = FALSE)
                    }, 
                    simplify = FALSE)
    wss <- sapply(ks, function(i) i$tot.withinss)
    imin <- which.min(wss)
    x@ic <- ks[[imin]]
    a <- c("Debris", rep("Clean", K - 1))
    x@assignments <- factor(a)
    
    return(x)
}

#' Get PCs
#' @importFrom irlba prcomp_irlba
#' @importFrom Matrix t
#' @export
get_pcs <- function(x, n_pcs = 50){
    countsv <- t(x@counts[x@vg,])
    countsv@x <- log10(countsv@x + 1)
    prcret <- prcomp_irlba(countsv[x@test_set,], 
                           n = n_pcs, 
                           center = TRUE, 
                           scale. = FALSE)
    x@pcs <- as.matrix(countsv %*% prcret$rotation)
    return(x)
}

