
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
get_knn <- function(x, nn=30, weighted=TRUE, verbose=FALSE){
    if (length(x@norm) == 0) stop("Normalize data before running nearest neighbors.")

    datf <- t(x@norm)

    droplets.use <- rownames(datf)

    if (verbose) cat(paste0("Finding ", as.character(nn), " nearest neighbors.\n"))

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
    if (verbose) cat("Done.\n")
    return(x)
}

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

    x@nn_graph <- g
    if (verbose) cat("Done.\n")
    return(x)
}

#' Initialize clustering for EM
#'
#' Given an SCE object, identify the cell types present in the top 
#' \code{cluster_n} droplets.
#'
#' Instead of randomly initializing the EM, cell types are estimated from 
#' droplets that are expected to contain cells/nuclei. The top 
#' \code{cluster_n} droplets are ranked by "gene" or "count", given by the 
#' parameter \code{order_by}. Then, droplets with at least those ranked 
#' counts/genes are included in the cluster set. The data is then normalized 
#' by taking first calculating the variable genes. A loess regression line 
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
#' \code{nn}. Clusters are identified form the KNN graph 
#' using the Louvain algorithm. Finally, only clusters with at least 
#' \code{min_size} (20 by default) droplets are considered cell types.
#'
#' @param x An SCE object.
#' @param cluster_n Numeric value specifying the number of droplets to use 
#'  in the test set.
#' @param order_by Whether to order the droplets by total number of total 
#'  counts or total number of genes detected.
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
                                cluster_n=1000, 
                                order_by="gene", 
                                use_var=TRUE, 
                                n_var=2000, 
                                lss=0.3, 
                                sf = "median", 
                                nn=30, 
                                min_size=20, 
                                verbose=TRUE){

    x <- set_cluster_set(x, 
                         cluster_n=cluster_n, 
                         order_by=order_by, 
                         verbose=verbose)
    if (use_var)  x <- get_var_genes(x, n_genes=n_var, lss=lss)
    x <- normalize_data(x, use_var=use_var, sf=sf)
    x <- get_knn(x, nn=nn, verbose=verbose)

    if (verbose) cat("Initializing clusters.\n")

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
    if (length(graph_clust) == 0) 
        stop("No clusters found during initialization. Try changing cluster_n, nn, or min_size parameters")

    # Specify cluster 1 as Debris
    all_clusters <- rep("1", length(x@bg_set))
    names(all_clusters) <- x@bg_set
    all_clusters <- c(graph_clust, all_clusters)
    all_clusters <- factor(all_clusters)
    # Change labels
    all_clusters <- factor(all_clusters, levels=levels(all_clusters), labels=seq(1, length(levels(all_clusters))))

    # Specify which labels are debris, which are cell types
    asgn <- rep("Clean", nlevels(all_clusters))
    names(asgn) <- levels(all_clusters)
    asgn[1] <- "Debris"
    asgn <- as.factor(asgn)

    if (verbose) cat(paste0("Initialized k=", as.character(nlevels(all_clusters)-1), " cell types and 1 debris cluster.\n"))

    x@ic <- list(clusters=all_clusters, 
                 assignments=asgn)
    return(x)
}

