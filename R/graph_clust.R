
#' @importFrom igraph graph_from_data_frame simplify cluster_louvain membership
#' @importFrom FNN get.knn
#' @export
knn <- function(mat, nn=30){
    knn.dat <- FNN::get.knn(as.matrix(mat), k=nn)
    knn.dat <- data.frame(from=rep(1:nrow(knn.dat$nn.index), nn), 
                          to=as.vector(knn.dat$nn.index), 
                          weight=1/(1 + as.vector(knn.dat$nn.dist)))
    return(knn.dat)
}

#' Get nearest neighbor graph
#'
#' @return An SCE object
#' @export
get_knn <- function(x, nn=30){
    if (length(x@pcs) == 0) stop("Get PCs before running KNN.")
    x@knn <- knn(x@pcs, nn=nn)
    return(x)
}

get_snn <- function(x, nn=50, kt=5, weighted=TRUE){
    if (length(x@pcs) == 0) stop("Get PCs before running NN.")
    if (verbose) cat(paste0("Finding shared nearest neighbors.\n"))
    
    droplets.use <- rownames(x@pcs)
    knn.dbscan <- dbscan::sNN(x@pcs[droplets.use,], k=nn, kt=kt)
    rownames(knn.dbscan$shared) <- droplets.use

    fnames <- as.vector(sapply(droplets.use, rep, nn))
    ti <- as.vector(t(knn.dbscan$id))
    tnames <- sapply(ti, function(i) droplets.use[i])

    if (weighted){
        knn.w <- as.vector(t(knn.dbscan$shared))
        knn.dat <- data.frame(from=fnames, to=tnames, weight=knn.w)
    } else {
        knn.dat <- data.frame(from=fnames, to=tnames)
    }
    knn.dat <- knn.dat[!is.na(knn.dat[,"to"]),]

    g <- igraph::graph_from_data_frame(knn.dat, directed=FALSE)
    g <- igraph::simplify(g)

    x@nn <- knn.dat
    x@nn_graph <- g
    if (verbose) cat("Done.\n")
    return(x)
}

    


#' Get louvain clusters
#'
#' @return Numeric vector of memberships
#' @importFrom igraph graph_from_data_frame simplify cluster_louvain membership
#' @importFrom FNN get.knn
#' @export
find_communities <- function(knn_df, nrun=1, method="louvain"){
    nw <- igraph::graph_from_data_frame(knn_df, directed = FALSE)
    nw <- igraph::simplify(nw)
    lc <- switch(method, 
                 louvain = igraph::cluster_louvain(nw), 
                 label_prop = igraph::cluster_label_prop(nw), 
                 walktrap = igraph::cluster_walktrap(nw), 
                 fast_greedy = igraph::cluster_fast_greedy(nw))
    return(igraph::membership(lc))
}

#' Normalize counts by simple trimmed mean of M
#' @export
tmm_counts <- function(counts, p=0.1){
    tmms <- apply(counts, 2, function(i){
                  return(mean(i[ i>=quantile(i,probs=c(p)) & i<=quantile(i,probs=c(1-p)) ]))
                 })
    counts <- sweep(counts, 2, tmms, "/")
    return(counts)
}

#' Get log2 fold change with TMM normalization
#' @export
get_logfc <- function(x){
    genes.use <- rownames(x@gene_data)[x@gene_data$exprsd]
    Test <- Matrix::rowSums(x@counts[genes.use, x@test_set])
    Debris <- Matrix::rowSums(x@counts[genes.use, x@bg_set])
    counts <- cbind(Debris, Test)
    counts <- tmm_counts(counts)
    counts <- apply(counts, 2, function(i) i/sum(i))

    logfc <- log2(counts[,"Debris"]) - log2(counts[,"Test"])
    logfc[is.na(logfc) | is.infinite(logfc)] <- 0
    ret <- cbind(counts, logfc)
    return(ret)
}

#' Initialize clusters for unlabeled data
#'
#' @return Numeric vector of memberships
#' @importFrom igraph graph_from_data_frame simplify cluster_louvain membership
#' @importFrom FNN get.knn
#' @export
initialize_clusters <- function(x, 
                                bf_thresh=10, 
                                gammas=c(10,10e6), 
                                tol_opt=100){
    x <- get_snn(x, nn=nn)

    if (length(x@nn) == 0) stop("Run nn before initializing clusters.")
    if (!"exprsd" %in% colnames(x@gene_data)) stop("Filter genes before initializing clusters.")

    genes.use <- rownames(x@gene_data)[x@gene_data$exprsd]


    gcl <- igraph::cluster_louvain(x@nn_graph)

    graph_clust <- find_communities(x@nn, method=method)
    graph_clust <- as.integer(graph_clust)
    graph_clust <- graph_clust - min(graph_clust)
    graph_clust <- graph_clust + 1 # Ensure starts at 1
    names(graph_clust) <- rownames(x@pcs)
    if (any(is.na(graph_clust))) stop("Returned cluster values have NA.")
    graph_clust <- graph_clust+1
    graph_clust <- factor(graph_clust)

    # Specify cluster 1 as Debris
    all_clusters <- rep("1", length(x@bg_set))
    names(all_clusters) <- x@bg_set
    all_clusters <- c(graph_clust, all_clusters)
    all_clusters <- as.factor(all_clusters)

    asgn <- rep("Clean", nlevels(all_clusters))
    names(asgn) <- levels(all_clusters)
    asgn[1] <- "Debris"
    asgn <- as.factor(asgn)

    genes_median <- sapply(levels(all_clusters), function(i) median(x@droplet_data[names(all_clusters)[all_clusters == i],"n_genes"]))
    print(genes_median)

    x@ic <- IC(graph=all_clusters, 
               assignments=asgn)
    return(x)
}

