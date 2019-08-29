
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

get_snn <- function(x, nn=30, kt=1){
    if (length(x@pcs) == 0) stop("Get PCs before running NN.")
    nn <- dbscan::sNN(x@pcs, k=nn, kt=kt)



#' @export
get_kclust <- function(x, K.max=20, B=500, verbose=FALSE){
    if (length(x@pcs) == 0) stop("Get PCs before running NN.")
    if (verbose) cat(paste0("Finding shared nearest neighbors.\n"))

    pcs <- x@pcs[x@cluster_set,,drop=FALSE]
    gskmn <- cluster::clusGap(pcs, FUN=kmeans, nstart=20, B=B, K.max=K.max, iter.max=100)
    best_k <- cluster::maxSE(gskmn$Tab[,"gap"], gskmn$Tab[,"SE.sim"])

    kr <- kmeans(pcs, centers=best_k, iter.max=100, nstart=100)

    # Specify cluster 1 as Debris
    all_clusters <- rep(1, length(x@bg_set))
    names(all_clusters) <- x@bg_set
    test_clust <- kr$cluster + 1
    all_clusters <- c(test_clust, all_clusters)
    all_clusters <- as.factor(all_clusters)

    asgn <- rep("Clean", nlevels(all_clusters))
    names(asgn) <- levels(all_clusters)
    asgn[1] <- "Debris"
    asgn <- as.factor(asgn)
    print(table(all_clusters))

    genes_median <- sapply(levels(all_clusters), function(i) median(x@droplet_data[names(all_clusters)[all_clusters == i],"n_genes"]))
    print(genes_median)

    x@clust_gap$Tab
    x@ic <- IC(graph=all_clusters, 
               assignments=asgn)

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

#' Initialize clusters for unlabeled data
#'
#' @return Numeric vector of memberships
#' @importFrom igraph graph_from_data_frame simplify cluster_louvain membership
#' @importFrom FNN get.knn
#' @export
initialize_clusters <- function(x, 
                                bf_thresh=10, 
                                min_size=30, 
                                gammas=c(10,10e6), 
                                tol_opt=100){
    if (length(x@nn) == 0) stop("Run nn before initializing clusters.")
    if (!"exprsd" %in% colnames(x@gene_data)) stop("Filter genes before initializing clusters.")

    genes.use <- rownames(x@gene_data)[x@gene_data$exprsd]

    gcl <- igraph::cluster_louvain(x@nn_graph)
    graph_clust <- gcl$membership
    names(graph_clust) <- gcl$names
    if (any(is.na(graph_clust))) stop("Returned cluster values have NA.")
    graph_clust <- graph_clust+1
    graph_clust <- factor(graph_clust)
    tb <- table(graph_clust)
    print(tb)
    keep <- names(tb)[tb >= min_size]
    print(keep)
    graph_clust <- graph_clust[graph_clust %in% keep]
    graph_clust <- as.factor(graph_clust)
    print(head(graph_clust))

    # Specify cluster 1 as Debris
    all_clusters <- rep("1", length(x@bg_set))
    names(all_clusters) <- x@bg_set
    all_clusters <- c(graph_clust, all_clusters)
    all_clusters <- as.factor(all_clusters)

    asgn <- rep("Clean", nlevels(all_clusters))
    names(asgn) <- levels(all_clusters)
    asgn[1] <- "Debris"
    asgn <- as.factor(asgn)
    print(table(all_clusters))

    genes_median <- sapply(levels(all_clusters), function(i) median(x@droplet_data[names(all_clusters)[all_clusters == i],"n_genes"]))
    print(genes_median)

    x@ic <- IC(graph=all_clusters, 
               assignments=asgn)
    return(x)
}

