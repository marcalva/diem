
#' Run k-means on multinomial
#' @export
kmeans_ss <- function(X, K = 30, labs, max_iter = 20){
    N <- nrow(X)
    G <- ncol(X)

    if (length(labs) != N) labs <- rep(0, N)

    means <- matrix(0, nrow = G, ncol = K)

    clusts <- 1:K
    labeled <- setdiff(unique(labs), 0)
    unlabeled <- setdiff(clusts, labeled)
    # Initialize
    for (l in labeled){
        means[,l] <- rowMeans(t(X[labs == l,]))
    }
    unl <- seq(N)[labs == 0]
    ksize <- min(K, floor(length(unl) / K))
    unl_s <- sample(unl, size = ksize * length(unlabeled))
    for (u in seq_along(unlabeled)){
        i1 <- (u - 1) * ksize + 1
        i2 <- u * ksize
        s <- unl_s[i1:i2]
        means[,unlabeled[u]] <- rowMeans(t(X[s,]))
    }

    p_l <- labs
    # Assignment
    unl <- seq(N)[labs == 0]
    Xmu <- as.matrix(X[unl,])
    # First iter
    d <- apply(Xmu, 1, function(x) colSums((x - means)^2))
    asgn <- apply(d, 2, which.min)
    p_l[unl] <- asgn

    # Wss
    message("WSS")
    wss <- lapply(clusts, function(i) sum((t(X[p_l == i,]) - means[,i])^2))
    wss <- do.call(sum, wss)
    message("WSS:", sum(wss))

    iter <- 1
    while(iter <= max_iter){
        n_l <- p_l
        # Calculate means
        for (i in clusts){
            means[,i] <- rowMeans(t(X[p_l == i,]))
        }
        # Calculate distances
        d <- apply(Xmu, 1, function(x) colSums((x - means)^2))
        asgn <- apply(d, 2, which.min)
        n_l[unl] <- asgn
        wss <- lapply(clusts, function(i) sum((t(X[p_l == i,]) - means[,i])^2))
        wss <- do.call(sum, wss)
        message("WSS:", sum(wss))
        if (all(p_l == n_l)) break
        p_l <- n_l
        iter <- iter + 1
    }
    if (iter > max_iter) message("Warning: did not converge")

    return(list("labels" = n_l, 
                "centers" = means, 
                "wss" = wss))

}

#' Get k-means initialization of PCs
#' @export
get_km <- function(x, K, n_start = 10, max_iter = 10){
    labs <- rep(0, nrow(x@pcs))
    names(labs) <- rownames(x@pcs)
    labs[intersect(names(labs), x@bg_set)] <- 1
    ks <- replicate(10, 
                    expr = {
                        kmeans_ss(x@pcs, 
                                  labs = labs, 
                                  K = K, 
                                  max_iter = max_iter)
                    }, 
                    simplify = FALSE)
    wss <- sapply(ks, function(i) i$wss)
    imin <- which.min(wss)
    x@ic <- ks[[imin]]
    return(x)
}

#' Get PCs
#' @importFrom irlba prcomp_irlba
#' @export
get_pcs <- function(x, K = 50){
    countsv <- t(x@counts[x@vg,])
    countsv@x <- log10(countsv@x + 1)
    prcret <- irlba::prcomp_irlba(countsv[x@test_set,], 
                                  n = K, 
                                  center = TRUE, 
                                  scale. = FALSE)
    x@pcs <- as.matrix(countsv %*% prcret$rotation)
    return(x)
}

#' Get k-nearest neighbor graph
#'
#' Run a KNN graph and store as an igraph object. The nearest neighbors 
#' are calculated using the dbscan implementation. The edges 
#' of the graph can be weighted with the inverse of the 
#' euclidean distance.
#' 
#' @param x An SCE object.
#' @param nn Number of nearest neighbors to use.
#' @param weighted Logical indicating whether to weigh the knn graph.
#' @param verbose verbosity
#'
#' @importFrom dbscan kNN
#' @importFrom igraph graph_from_data_frame simplify
#' @return An SCE object
get_knn <- function(x, nn = 30, weighted = TRUE, verbose = FALSE){
    if (length(x@norm) == 0) stop("Normalize data before running nearest neighbors.")

    datf <- t(x@norm)

    droplets.use <- rownames(datf)

    if (verbose) message("finding ", nn, " nearest neighbors")

    knn.dbscan <- kNN(datf, k=nn)

    fnames <- as.vector(sapply(rownames(datf), rep, nn))
    ti <- as.vector(t(knn.dbscan$id))
    tnames <- sapply(ti, function(i) rownames(datf)[i])

    if (weighted){
        knn.w <- 1/as.vector(t(knn.dbscan$dist))
        knn.dat <- data.frame(from=fnames, to=tnames, weight=knn.w)
    } else {
        knn.dat <- data.frame(from=fnames, to=tnames)
    }

    g <- graph_from_data_frame(knn.dat, directed=FALSE)
    g <- simplify(g)

    x@nn_graph <- g
    return(x)
}

#' Initialize clustering for EM
#'
#' Given an SCE object, identify the cell types present.
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
#' @param use_var A logical indicating whether to subset the data to include
#'  only variable genes. This overrides \code{genes.use}.
#'  The default is TRUE as it may better identify cell types.
#' @param n_var Number of variable genes to use.
#' @param lss Numeric value of the span parameter of the loess regression. 
#' @param sf Either a numeric scaling factor to multiply counts after 
#'  division by column sums, or "median" indicating to multiply by the median 
#'  number of total read/UMI counts in droplets (default).
#' @param nn Number of nearest neighbors to calculate in constructing the 
#'  graph.
#' @param min_size Numeric value giving the minimum number of droplets in 
#'  cluster for it to be used for initialization as a cell type for EM. 
#' @param verbose verbosity.
#'
#' @return An SCE object
#' @importFrom igraph cluster_louvain
#' @export
initialize_clusters <- function(x, 
                                use_var = TRUE, 
                                n_var = 2000, 
                                lss = 0.3, 
                                sf = "median", 
                                nn = 30, 
                                min_size = 20, 
                                verbose = FALSE){

    if (length(x@cluster_set) == 0) 
        stop("0 droplets in cluster_set. Set droplets for initialization with set_cluster_set first")
    if (use_var)  x <- get_var_genes(x, n_genes = n_var, lss = lss)
    x <- normalize_data(x, use_var = use_var, sf = sf)
    x <- get_knn(x, nn = nn, verbose = verbose)

    if (verbose) message("initializing clusters")

    genes.use <- rownames(x@gene_data)[x@gene_data$exprsd]

    # Set louvain clusters
    gcl <- cluster_louvain(x@nn_graph)
    graph_clust <- gcl$membership
    names(graph_clust) <- gcl$names
    if (any(is.na(graph_clust))) stop("Returned cluster values have NA.")
    graph_clust <- graph_clust+1
    tb <- table(graph_clust)
    keep <- names(tb)[tb >= min_size]
    graph_clust <- graph_clust[graph_clust %in% keep]
    if (length(graph_clust) == 0){ 
        stop("No clusters found during initialization. Try changing ", 
             sQuote("cluster_n"), sQuote("nn"), " or ", 
             sQuote("min_size"), "parameters")
    }

    # Specify cluster 1 as Debris
    all_clusters <- rep("1", length(x@bg_set))
    names(all_clusters) <- x@bg_set
    all_clusters <- c(graph_clust, all_clusters)
    all_clusters <- factor(all_clusters)
    # Change labels
    all_clusters <- factor(all_clusters, 
                           levels=levels(all_clusters), 
                           labels=seq(1, length(levels(all_clusters))))

    # Specify which labels are debris, which are cell types
    asgn <- rep("Clean", nlevels(all_clusters))
    names(asgn) <- levels(all_clusters)
    asgn[1] <- "Debris"
    asgn <- as.factor(asgn)

    if (verbose){
        message("initialized k=", nlevels(all_clusters)-1, 
                " cell types and 1 debris cluster")
    }

    x@ic <- list(clusters=all_clusters, 
                 assignments=asgn)
    return(x)
}

