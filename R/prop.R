
#' Estimate proportions of each cluster given counts
#'
#' Estimate the proportion of each cell type within a droplet using 
#' the multinomial parameters of the fitted model.
#'
#' @param x A sparse count matrix.
#' @param means Means for each cluster under the multinomial model.
#' @param w_init Initializing values for the proportions.
#' @param lrate Learning rate.
#' @param eps Threshold for convergence.
#' @param threads Number of threads.
#' @param verbose Verbosity.
#'
#' @return A matrix containing fractions of each cluster.
#'
#' @export
get_clust_dat <- function(counts, 
                          means, 
                          w_init, 
                          lrate = 0.1,
                          eps = 0.001,
                          threads = 1,
                          verbose = TRUE){
    w <- fit_prop(counts, 
                  means, 
                  w_init, 
                  lrate = lrate, 
                  eps = eps, 
                  threads = threads, 
                  display_progress = verbose)
    colnames(w) <- colnames(means)
    rownames(w) <- colnames(counts)
    return(w)
}

#' Estimate proportions of each cluster for each droplet
#'
#' Estimate the proportion of each cell type within a droplet using 
#' the multinomial parameters of the fitted model.
#'
#' @param x An SCE object.
#' @param droplets Specify droplets to calculate proportions for. If
#'  not specified, calculates for all test droplets, which can be slow.
#' @param lrate Learning rate.
#' @param eps Threshold for convergence.
#' @param threads Number of threads
#' @param verbose Verbosity.
#'
#' @return A matrix containing fractions of each cluster.
#'
#' @export
get_clust_f <- function(x, 
                        droplets = NULL, 
                        lrate = 0.1, 
                        eps = 0.001, 
                        threads = 1, 
                        verbose = TRUE){
    td <- droplet_data(x)
    if (! "Cluster" %in% colnames(td)){
        stop("Assign clusters before running get_clust_f")
    }
    if (is.null(droplets)){
        droplets <- rownames(td)
    }
    droplets <- intersect(droplets, rownames(td))

    clusts <- sort(as.numeric(unique(td$Cluster)))
    names(clusts) <- as.character(clusts)

    genes <- rownames(x@gene_data)[x@gene_data$exprsd]
    A <- x@model$params$Alpha[genes,clusts]
    colnames(A) <- names(clusts)
    rc <- raw_counts(x)[genes, droplets]
    Pi <- x@model$params$Pi[clusts]
    llk <- x@model$llk[droplets,clusts]
    w_init <- get_z(llk, Pi)
    colnames(w_init) <- names(clusts)

    if (verbose){
        message("getting proportions for ", ncol(rc), " droplets ", 
                " with ", nrow(rc), " genes")
    }
    w <- get_clust_dat(rc, A, w_init, lrate = lrate, eps = eps, 
                       threads = threads, verbose = verbose)
    colnames(w) <- names(clusts)
    rownames(w) <- colnames(rc)

    return(w)
}

