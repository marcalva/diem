
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
initialize_clusters <- function(x, method="louvain", cor_thresh=0.95){

    if (length(x@knn) == 0) stop("Run knn before initializing clusters.")
    if (!"exprsd" %in% colnames(x@gene_data)) stop("Filter genes before initializing clusters.")

    genes.use <- rownames(x@gene_data)[x@gene_data$exprsd]

    graph_clust <- find_communities(x@knn, method=method)
    graph_clust <- as.integer(graph_clust)
    graph_clust <- graph_clust - min(graph_clust)
    graph_clust <- graph_clust + 1 # Ensure starts at 1
    names(graph_clust) <- rownames(x@pcs)
    if (any(is.na(graph_clust))) stop("Returned cluster values have NA.")
    graph_clust <- graph_clust+1
    graph_clust <- factor(graph_clust)

    all_clusters <- rep("1", ncol(x@counts))
    names(all_clusters) <- colnames(x@counts)
    all_clusters[names(graph_clust)] <- as.character(graph_clust)
    all_clusters <- as.factor(all_clusters)

    logfc <- get_logfc(x)
    # logfc <- logfc[logfc[,3] > 0,]

    # Get TMM of test clusters
    print(tapply(x@droplet_data[names(graph_clust),"n_genes"], graph_clust, median))
    cg <- tapply(x@droplet_data[names(graph_clust),"n_genes"], graph_clust, median)
    cm <- tapply(x@droplet_data[names(graph_clust),"total_counts"], graph_clust, mean)
    min_clust <- names(cm)[which.min(cm)]
    clust_counts <- sapply(levels(graph_clust), function(g) {
                           dn <- names(graph_clust)[graph_clust == g]
                           Matrix::rowSums(x@counts[genes.use, dn])
                 })
    clust_tmm <- tmm_counts(clust_counts)
    clust_p <- apply(clust_tmm, 2, function(i) i/sum(i))

    # Get TMM of ref
    Debris <- Matrix::rowSums(x@counts[genes.use, x@bg_set])
    Test <- Matrix::rowSums(x@counts[genes.use, x@test_set])
    ref_counts <- cbind(Debris, Test)
    ref_tmm <- tmm_counts(ref_counts)
    ref_p <- apply(ref_tmm, 2, function(i) i/sum(i))

    # DB score
    debris_sc <- logfc[,1] %*% logfc[,3]
    clean_sc <- logfc[,2] %*% logfc[,3]

    ref_sc <- t(ref_p[rownames(logfc),]) %*% logfc[,3]
    clust_sc <- t(clust_p[rownames(logfc),]) %*% logfc[,3]

    # Correlate log1p of TMM
    cor_val <- cor(log1p(clust_tmm[,min_clust]), log1p(ref_tmm[,"Debris"]))
    print(cor_val)

    # Assign new clusters
    map <- levels(all_clusters); names(map) <- levels(all_clusters)
    # if (cor_val >= cor_thresh) map[min_clust] <- "1"
    # merged_clusters <- factor(all_clusters, levels=levels(all_clusters), labels=map)
    # merged_clusters <- factor(merged_clusters, levels=levels(merged_clusters), labels=as.character(1:nlevels(merged_clusters)))
    # Confirm re-labeling works
    # print(table(all_clusters, merged_clusters))

    merged_clusters <- all_clusters
    asgn <- rep("Clean", nlevels(merged_clusters))
    names(asgn) <- levels(merged_clusters)
    asgn[1] <- "Debris"
    asgn <- as.factor(asgn)

    # Store in IC class
    x@ic <- IC(graph=all_clusters, 
               map=map, 
               scores=cor_val, 
               merged=merged_clusters, 
               assignments=asgn)

    return(x)
}
