
#' Return design matrix of labs vector
#'
#' @param labs Vector of labels. Ordered by factor and corresponds to 
#'  columns of returned matrix.
#'
#' @return A droplet by cluster matrix
z_table <- function(labs){
    labs <- factor(labs)
    nl <- nlevels(labs)
    z <- matrix(0, nrow = length(labs), ncol = nl)
    for (i in 1:nl){
        z[unclass(labs) == i,i] <- 1
    }
    rownames(z) <- names(labs)
    colnames(z) <- levels(labs)
    return(z)
}

#' Run PCA
#'
#' Get PCs from normalized counts.
#' 
#' @param x An SCE object.
#' @param n_pcs Number of PCs to calculate.
#'
#' @return An SCE object with PCs
#'
#' @importFrom irlba prcomp_irlba
#' @importFrom stats prcomp
#' @importFrom Matrix t
run_pca <- function(x, n_pcs = 30){
    m <- nrow(x@norm[x@vg,])
    n <- ncol(x@norm[x@vg,])
    if (n_pcs >= 0.5 * min(m, n)){
        pcs <- prcomp(t(x@norm[x@vg,]), center = TRUE, scale. = TRUE)
    } else {
        pcs <- prcomp_irlba(t(x@norm[x@vg,]), n = n_pcs, center = TRUE, scale. = TRUE)
    }
    pcs <- pcs$x
    rownames(pcs) <- colnames(x@norm[x@vg,])
    x@pcs <- pcs[,1:n_pcs,drop=FALSE]
    return(x)
}

#' Get PCs
#'
#' Run PCA and get top \code{n_pcs} PCs for the test set. 
#' The PCs are used as the features for the initial k-means clustering. 
#' The counts data for the test set are count-normalized to the median 
#' and log transformed. Then the top \code{n_var_genes} 
#' variable genes are calculated using the function 
#' \code{\link{get_var_genes}}. PCA is run on the normalized count data 
#' for these variable genes only.
#'
#' @param x An SCE object.
#' @param n_var_genes Number of top variable genes to use for PCA.
#' @param lss The span parameter of the loess regression, the parameter 
#'  for the function \code{\link[stats]{loess}}. The loess regression is
#'  used to regress out the effect of mean expression on variance.
#' @param n_pcs Number of PCs to return.
#'
#' @return An SCE object with PCs
#'
#' @export
get_pcs <- function(x, 
                    n_var_genes = 2000, 
                    lss = 0.3, 
                    n_pcs = 30){
    x <- get_var_genes(x, n_genes = n_var_genes, lss = lss)
    x <- normalize_data(x)
    x <- run_pca(x, n_pcs = n_pcs)
    return(x)
}

#' Re-assign labels by condensing to 1 to k
#' 
#' @param labs Vector of labels.
#'
#' @return Vector of labels in which the numeric values are contiguous.
condense_labs <- function(labs){
    cts <- setdiff(sort(unique(labs)), "1")
    to <- c("1", as.character(seq(2, length(cts) + 1)))
    from <- c("1", cts)
    names(to) <- from
    nms <- names(labs)
    labs <- to[as.character(labs)]
    names(labs) <- nms
    return(labs)
}

#' Merge small clusters
#' 
#' Merge clusters given in labels until all clusters have a size of at 
#' least \code{min_size}. The \code{dat} matrix contains the data used 
#' to generate the cluster means.
#' Data points are assigned to the next closest cluster.
#'
#' @param labs A vector of labels, with numeric values corresponsing 
#'  to clusters
#' @param dat A feature by observation data matrix used to get cluster means.
#' @param min_size The minimum number of members that belong to a cluster.
#' @param verbose Verbosity.
#'
#' @importFrom stats dist
#' 
#' @return A vector containing new labels
merge_size <- function(labs, dat, min_size = 10, verbose = TRUE){
    labs <- as.character(labs)
    freq <- table(labs)
    nm <- 0
    while(min(freq) < min_size){
        clusts <- names(freq)
        small_clust <- names(freq)[which.min(freq)]

        # Find closest clusters
        means <- sapply(clusts, function(i) colMeans(dat[labs == i,,drop=FALSE]))
        dists <- as.matrix(dist(t(means), upper = TRUE, diag = TRUE))
        others <- colnames(dists)[colnames(dists) != small_clust]
        dists <- dists[others,,drop=FALSE]
        closest <- rownames(dists)[which.min(dists[,small_clust])]
        labs[labs == small_clust] <- closest
        freq <- table(labs)
        nm <- nm + 1
    }
    if (verbose){
        if (nm > 0){
            message("merged ", nm, 
                    " small ", 
                    ngettext(nm, "cluster", "clusters"))
        }
    }

    labs <- condense_labs(labs)
    return(labs)
}

#' Initialize clusters
#' 
#' Initialize the parameters of the Dirichlet-mulitnomial 
#' mixture model. First 
#' run k-means on the principal components of the test set to 
#' initialize the cluster memberships. Then use these memberships 
#' to estimate alpha and pi, the parameters of the mixture model. 
#' These clusters are pruned by removing those that are close in 
#' to the background distribution using a likelihood-based 
#' distance metric. The parameter \code{k_init} specifies 
#' initial number of clusters k to set for k-means.
#'
#' @param x An SCE object.
#' @param k_init The number of clusters to initialize k-means.
#' @param iter.max_init The maximum number of k-means interations 
#'  for the initialization.
#' @param nstart_init The number of starts to use in k-means for 
#'  for the initialization.
#' @param min_size_init The minimum number of droplets that must belong 
#'  to an initialized cluster.
#' @param seedn The seed for random k-means initialization. 
#'  It is set to 1 by default. If you desire truly random initializations
#'  across runs, set to NULL or different values for each run.
#' @param psc Pseudocount to add to all alpha parameter values to avoid 0.
#' @param verbose Verbosity.
#'
#' @return An SCE object
#'
#' @importFrom Matrix Diagonal colSums rowMeans rowSums t
#' @importFrom stats kmeans
#'
#' @export
init <- function(x, 
                 k_init = 30, 
                 iter.max_init = 15, 
                 nstart_init = 30, 
                 min_size_init = 10, 
                 seedn = 1, 
                 psc = 1e-10, 
                 verbose = TRUE){ 
    genes.use <- rownames(x@gene_data)[x@gene_data$exprsd]

    sizes <- colSums(x@counts[genes.use,])

    if (length(x@pcs) == 0){
        stop("Calculate PCs before k-means initialization")
    }

    pcs_s <- scale(x@pcs)

    labs <- rep(0, ncol(x@counts[genes.use,]))
    names(labs) <- colnames(x@counts[genes.use,])
    labs[x@bg_set] <- 1

    all_ret <- list()
    set.seed(seedn)

    for (k in k_init){
        kc <- as.character(k)
        if (verbose){
            message("initializing parameters for k_init = ", k_init)
        }

        old_opt <- options(warn = -1)
        km <- kmeans(x = pcs_s, 
                     centers = k, 
                     iter.max = iter.max_init, 
                     nstart = nstart_init)
        options(old_opt)
        kclust <- km$cluster # integer vector of 1 to k values, corresponding to column in pcs_s

        # Merge small clusters from k-means
        kclust <- merge_size(kclust, pcs_s, min_size_init, verbose)
        kclust <- as.numeric(kclust)
        names(kclust) <- rownames(pcs_s)

        labs_k <- labs
        labs_k[names(kclust)] <- kclust + 1

        Z <- z_table(labs_k)
        Alpha <- get_alpha(x@counts[genes.use,], Z)
        Pi <- get_pi(Z)
        params <- list("Alpha" = Alpha, "Pi" = Pi)
        x@kruns[[kc]] <- list()
        x@kruns[[kc]] <- list("params" = params, "Z" = Z)
        # x@kruns[[kc]][[1]] <- list("params" = params)
    }
    return(x)
}

