
#' Run a basic DE between clusters
#'
#' @param means A gene by cluster sata frame of cluster means, 
#'  of which the columns sum to 1.
#' @param top_n Number of top DE genes to return.
#' 
#' @return A data frame of DE results
de_basic <- function(means, top_n = 20){
    clusters <- 1:ncol(means)
    de <- lapply(clusters, function(cl){
                 cln <- setdiff(clusters, cl)
                 p1 <- means[,cl]
                 p2 <- rowMeans(means[,cln,drop=FALSE])
                 di <- p1 - p2
                 logfc <- log2(p1) - log2(p2)
                 datf <- data.frame("diff" = di, 
                                    "logFC" = logfc,
                                    "p1" = p1, 
                                    "p2" = p2, 
                                    "cluster" = cl, 
                                    "gene" = rownames(means))
                 o <- order(di, decreasing = TRUE)
                 datf <- datf[o,]
                 return(datf[1:top_n,])
    })
    ret <- do.call(rbind, de)
    ret$gene <- as.character(ret$gene)
    rownames(ret) <- NULL
    return(ret)
}

#' Get top up-regulated genes per cluster
#'
#' Extract marker genes for each cluster by ranking them based on the 
#' difference of their proportion means from the Dirichlet 
#' distributions estimated by DIEM. Output the top \code{top_n} 
#' genes for each cluster.
#'
#' @param x An SCE object.
#' @param top_n Number of top DE genes to return.
#' @param k_init Run EM on the \code{k_init} initialization(s). 
#'  If NULL (default), run on all \code{k_init} initializations.
#'
#' @return An SCE object.
#'
#' @export
top_genes <- function(x, top_n = 20, k_init = NULL){
    genes.use <- rownames(x@gene_data)[x@gene_data$exprsd]
    k_init <- check_k_init(x, k_init)
    means <- Alpha(x, k_init)
    means <- means[genes.use,]
    means <- sweep(means, 2, colSums(means), "/")
    de <- de_basic(means, top_n = top_n)
    return(de)
}

#' Estimate percent of cell types in droplets
#'
#' Use results from DIEM to estimate the 
#' cell type proportions within each droplet using nnls
#'
#' @param x An SCE object.
#' @param show_progress Boolean indicating whether to show a progress
#'  bar.
#' @param k_init Run EM on the \code{k_init} initialization(s). 
#'  If NULL (default), run on all \code{k_init} initializations.
#'
#' @return An SCE object, with the prop slot containing a
#'  cluster by droplet matrix of cell type proportions
#'
#' @importFrom nnls nnls
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
estim_prop_nnls <- function(x, show_progress = TRUE, k_init = NULL){
    genes.use <- rownames(x@gene_data)[x@gene_data$exprsd]
    droplets.use <- x@test_set
    k_init <- check_k_init(x, k_init)

    means <- Alpha(x, k_init)
    means <- means[genes.use,]
    means <- sweep(means, 2, colSums(means), "/")

    counts <- x@counts[genes.use, droplets.use]
    counts <- divide_by_colsum(counts)

    if (show_progress)
        pb <- txtProgressBar(min = 0, max = ncol(counts), style = 3)

    props <- lapply(1:ncol(counts), function(i){
                    if (show_progress)
                        setTxtProgressBar(pb, i)

                    ret <- nnls(means, counts[,i])$x
                    return(ret)
    })
    if (show_progress)
        close(pb)

    prop <- do.call(cbind, props)
    colnames(prop) <- droplets.use
    x@prop <- sweep(prop, 2, colSums(prop), "/")
    return(x)
}

#' Estimate percent of debris in droplets
#'
#' Given the debris clusters in deb_clust, estimate the percent 
#' of debris within each droplet. Store the percent in the 
#' test_data slot.
#'
#' @param x An SCE object.
#' @param deb_clust A numeric vector of debris clusters. Contains 
#'  only cluster 1 by default, and should always contain cluster 1.
#' @param name Column name in test_data to store debris percent in.
#' 
#' @return An SCE object
#' @export
estimate_db <- function(x, deb_clust = c(1), name = "pct.debris"){
    if ( length(x@prop) == 1 ) 
        stop("run estim_prop_nnls before estimate_db")
    if ( any( !(deb_clust %in% 1:nrow(x@prop)) ) )
        stop("deb_clust clusters must be between 1 and ", nrow(x@prop))

    x@test_data[colnames(x@prop),name] <- 100 * colSums(x@prop[deb_clust,])
    return(x)
}

#' Correct counts by subtracting estimated counts from debris
#'
#' Remove the counts predicted to originate from debris clusters. This 
#' function takes the cluster proportions for each droplet, and estimated 
#' the debris counts by adding the expresion means for the debris clusters 
#' weighted by the their proportions mentioned above, and subtracting this
#' value from the counts. The debris clusters are specified by deb_clust 
#' and default to cluster 1, which should always be included. The values to 
#' subtract are rounded according to the function give by the parameter rfunc, 
#' which is round by default. This function can be floor, so that only debris 
#' counts greater than or equal to 1 are subtracted. The rfunc function 
#' can also be set to NULL, so that no rounding is performed.
#' 
#' @param x An SCE object.
#' @param deb_clust A numeric vector of debris clusters. Contains 
#'  only cluster 1 by default, and should always contain cluster 1.
#' @param round_count Boolean indiciating whether to round the counts 
#'  to integers or not.
#' @param k_init Run EM on the \code{k_init} initialization(s). 
#'  If NULL (default), run on all \code{k_init} initializations.
#' 
#' @importFrom Matrix colSums
#' @export
correct_counts <- function(x, 
                           deb_clust = c(1), 
                           round_count = TRUE, 
                           k_init = NULL){
    genes.use <- rownames(x@gene_data)[x@gene_data$exprsd]
    droplets.use <- x@test_set
    k_init <- check_k_init(x, k_init)

    if ( length(x@prop) == 1 ) 
        stop("run estim_prop_nnls before estimate_db")

    means <- Alpha(x, k_init)
    means <- means[genes.use,]
    means <- sweep(means, 2, colSums(means), "/")

    counts <- x@counts[genes.use, droplets.use]
    total_count <- colSums(counts)

    cnorm <- fast_correct(counts, 
                          means[, deb_clust, drop=FALSE], 
                          x@prop[deb_clust, , drop=FALSE], 
                          round_count = round_count)

    rownames(cnorm) <- genes.use
    colnames(cnorm) <- droplets.use
    if (round_count)
        x@corrected <- as(cnorm, "CsparseMatrix")
    else
        x@corrected <- cnorm
    return(x)
}

