
#' Get k-nearest neighbor graph
#'
#' @return An SCE object
#' @export
get_knn <- function(x, nn=30, weighted=TRUE, verbose=FALSE){
    if (length(x@norm) == 0) stop("Normalize data before running nearest neighbors.")

    datf <- t(x@norm)
    droplets.use <- rownames(datf)

    if (verbose) cat(paste0("Finding ", as.character(nn), " nearest neighbors.\n"))

    knn.dbscan <- dbscan::kNN(datf, k=nn)

    fnames <- as.vector(sapply(rownames(datf), rep, nn))
    ti <- as.vector(t(knn.dbscan$id))
    tnames <- sapply(ti, function(i) rownames(datf)[i])

    if (weighted){
        knn.w <- 1/as.vector(t(knn.dbscan$dist))
        knn.dat <- data.frame(from=fnames, to=tnames, weight=knn.w)
    } else {
        knn.dat <- data.frame(from=fnames, to=tnames)
    }

    g <- igraph::graph_from_data_frame(knn.dat, directed=FALSE)
    g <- igraph::simplify(g)

    x@nn_df <- knn.dat
    x@nn_graph <- g
    if (verbose) cat("Done.\n")
    return(x)
}

#' Get shared nearest neighbors
#' @export
get_snn <- function(x,
                    nn=30,
                    kt=3,
                    min_counts=10,
                    droplets.use=NULL,
                    weighted=TRUE,
                    algorithm="kd_tree",
                    verbose=FALSE){
    if (is.null(droplets.use)) droplets.use <- colnames(x@counts)
    if (length(x@norm) == 0) stop("Normalize data before running nearest neighbors.")
    datf <- as.matrix(t(x@norm))
    droplets.use <- rownames(datf)

    if (verbose) cat(paste0("Finding shared nearest neighbors.\n"))
    # knn.dbscan <- dbscan::sNN(x@pcs, k=nn, kt=2)
    knn.dbscan <- dbscan::sNN(datf, k=nn, kt=kt)
    rownames(knn.dbscan$shared) <- rownames(datf)
    knn.dbscan$id <- knn.dbscan$id
    knn.dbscan$shared <- knn.dbscan$shared

    fnames <- as.vector(sapply(rownames(datf), rep, nn))
    ti <- as.vector(t(knn.dbscan$id))
    tnames <- sapply(ti, function(i) rownames(datf)[i])

    if (weighted){
        knn.w <- as.vector(t(knn.dbscan$shared))
        knn.dat <- data.frame(from=fnames, to=tnames, weight=knn.w)
    } else {
        knn.dat <- data.frame(from=fnames, to=tnames)
    }
    knn.dat <- knn.dat[!is.na(knn.dat[,"to"]),]

    g <- igraph::graph_from_data_frame(knn.dat, directed=FALSE)
    g <- igraph::simplify(g)

    x@nn_df <- knn.dat
    x@nn_graph <- g
    if (verbose) cat("Done.\n")
    return(x)
}

norm_tmm <- function(counts, qnt=0.1){
    counts <- apply(counts, 2, function(i){
                    i_mean <-sum(i[i >= quantile(i, qnt) & i <= quantile(i, 1-qnt)])
                    return(i_mean*i/sum(i))
                    })
    return(counts)
}

score_cluster <- function(x, clust){
    counts <- x@counts[x@gene_data$exprsd,]
    clust <- as.factor(clust)
    ge <- sapply(levels(clust), function(i) Matrix::rowSums(counts[,names(clust)[clust == i]]))
    ge <- apply(ge, 2, function(i) i/sum(i))

    ref_e <- sapply(list(x@bg_set, x@test_set), function(i) Matrix::rowSums(counts[,colnames(counts) %in% i]))
    ref_e <- apply(ref_e, 2, function(i) i/sum(i))
    genes <- rownames(ref_e)[ apply(ref_e, 1, function(i) all(i > 0)) ]
    log2fc <- log2(ref_e[genes,1]) - log2(ref_e[genes,2])

    ge_scores <- apply(ge[genes,], 2, function(i) i %*% log2fc)
    return(ge_scores)
}

#' Initialize clusters for unlabeled data
#'
#' @return Numeric vector of memberships
#' @importFrom igraph graph_from_data_frame simplify cluster_louvain membership
#' @importFrom FNN get.knn
#' @export
initialize_clusters <- function(x, 
                                cluster_n=500, 
                                order_by="gene", 
                                nn=30, 
                                min_size=20, 
                                verbose=TRUE){

    x <- set_cluster_set(x, cluster_n=cluster_n, order_by=order_by, verbose=verbose)
    x <- normalize_data(x)
    x <- get_knn(x, nn=nn, verbose=verbose)
 

    if (length(x@nn_df) == 0) stop("Run nn before initializing clusters.")
    if (!"exprsd" %in% colnames(x@gene_data)) stop("Filter genes before initializing clusters.")

    if (verbose) cat("Initializing clusters.\n")

    genes.use <- rownames(x@gene_data)[x@gene_data$exprsd]

    gcl <- igraph::cluster_louvain(x@nn_graph)
    graph_clust <- gcl$membership
    names(graph_clust) <- gcl$names
    if (any(is.na(graph_clust))) stop("Returned cluster values have NA.")
    graph_clust <- graph_clust+1
    tb <- table(graph_clust)
    keep <- names(tb)[tb >= min_size]
    graph_clust <- graph_clust[graph_clust %in% keep]
    if (length(graph_clust) == 0) stop("No clusters found during initialization. Try increasing cluster_n and/or nn parameter(s)")

    # Specify cluster 1 as Debris
    all_clusters <- rep("1", length(x@bg_set))
    names(all_clusters) <- x@bg_set
    all_clusters <- c(graph_clust, all_clusters)
    all_clusters <- as.factor(all_clusters)
    all_clusters <- factor(all_clusters, levels=levels(all_clusters), labels=seq(1, length(levels(all_clusters))))

    asgn <- rep("Clean", nlevels(all_clusters))
    names(asgn) <- levels(all_clusters)
    asgn[1] <- "Debris"
    asgn <- as.factor(asgn)

    if (verbose) cat(paste0("Initialized k=", as.character(nlevels(all_clusters)-1), " cell types and 1 debris cluster.\n"))

    x@ic <- IC(graph=all_clusters, 
               assignments=asgn)
    return(x)
}

