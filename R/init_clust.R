
#' @export
init_alpha <- function(x, Z, psc = 1e-10, max_iter = 500, eps = 1e-4){
    if (nrow(Z) != ncol(x)){
        stop("Number of columns in x must match number of rows in Z")
    }
    clusts <- 1:ncol(Z)
    clust_max <- apply(Z, 1, which.max)

    # test a0
    sizes <- colSums(x)
    alphas <- get_alpha(t(x), Z, max_loo = max_iter, eps = eps)
    colnames(alphas) <- clusts
    return(alphas)
}

init_pi <- function(Z){
    Pi <- colSums(Z)
    Pi <- Pi / sum(Pi)
    return(Pi)
}

#' @export
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


#' Get PCs
#' @importFrom irlba prcomp_irlba
#' @importFrom Matrix t
#' @export
get_pcs <- function(x, n_pcs = 30){
    #countsv <- t(x@counts[x@vg,])
    #countsv@x <- log10(countsv@x + 1)
    #prcret <- prcomp_irlba(countsv[x@test_set,], 
    #                       n = n_pcs, 
    #                       center = TRUE, 
    #                       scale. = FALSE)
    #x@pcs <- as.matrix(countsv %*% prcret$rotation)

    pcs <- prcomp_irlba(t(x@norm[x@vg,]), n = n_pcs, center = TRUE, scale. = TRUE)
    pcs <- pcs$x
    rownames(pcs) <- colnames(x@norm)
    x@pcs <- pcs
    return(x)
}

#' @export
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
#' @params labs A vector of labels.
#' @param dat A matrix containing the center means of the groups in each column.
#' @params min_size The minimum number of members that belong to a cluster.
#' @params verbose Verbosity.
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
#' @params x An SCE object.
#' @params n_pcs The number of PCs to use for the k-means.
#' @params k_init The number of clusters to initialize k-means.
#' @params iter.max The maximum number of k-means interations.
#' @params nstart The number of starts to use in k-means.
#' @params min_size The minimum number of members that belong to a cluster.
#' @params fltr The filter threshold that controls the minimum distance to 
#'  the background distribution that a cluster must have. Remove those 
#'  centers with a distance less than this value.
#' @params seedn The seed for random k-means initialization. 
#'  It is set to 1 by default. If you desire random initializations, 
#'  set to NULL, or different values across runs.
#' @params psc Pseudocount to add to all alpha parameter values to avoid 0.
#' @params verbose Verbosity.
#'
#' @return An SCE object
#' @importFrom Matrix Diagonal colSums rowMeans rowSums t
#' @export
init <- function(x, 
                 n_pcs = 30, 
                 k_init = 30, 
                 iter.max = 15, 
                 nstart = 30, 
                 min_size = 10, 
                 fltr = 0.1, 
                 seedn = 1, 
                 psc = 1e-10, 
                 verbose = TRUE){ 
    genes.use <- rownames(x@gene_data)[x@gene_data$exprsd]

    counts <- x@counts[genes.use,]
    countst <- t(counts)
    sizes <- colSums(counts)

    # Get PCs
    if (verbose){
        message("Getting PCs")
    }
    x <- get_var_genes(x)
    x <- normalize_data(x)
    x <- get_pcs(x, n_pcs = n_pcs)
    pcs_s <- scale(x@pcs)

    labs <- rep(0, ncol(counts))
    names(labs) <- colnames(counts)
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
                     iter.max = iter.max, 
                     nstart = nstart)
        options(old_opt)
        kclust <- km$cluster # integer vector of 1 to k values, corresponding to column in pcs_s

        # Merge small clusters from k-means
        kclust <- merge_size(kclust, pcs_s, min_size, verbose)
        kclust <- as.numeric(kclust)
        names(kclust) <- rownames(pcs_s)

        labs_k <- labs
        labs_k[names(kclust)] <- kclust + 1

        Z <- z_table(labs_k)
        Alpha <- get_alpha(countst, Z)
        Pi <- init_pi(Z)
        params <- list("Alpha" = Alpha, "Pi" = Pi)
        #x@kruns[[kc]] <- list()
        #x@kruns[[kc]] <- list("params" = params, "Z" = Z)
        x@init[[kc]] <- list()
        x@init[[kc]][[1]] <- list("params" = params, "Z" = Z)
        x <- get_dist(x, verbose = verbose)
        x <- rm_close(x, k_init = k, fltr = fltr, verbose = verbose)
    }
    return(x)
}

rm_close <- function(x, 
                    k_init = NULL, 
                    fltr = 0.1, 
                    verbose = TRUE){ 

    genes.use <- rownames(x@gene_data)[x@gene_data$exprsd]
    counts <- x@counts[genes.use,]

    if (is.na(fltr)){
        stop("fltr value is set to NA, change to a value between 0 and 1")
    }

    all_ret <- list()
    if (is.null(k_init)){
        k_init <- names(x@init)
    }

    for (k in k_init){
        if (verbose){
            message("removing clusters for k_init = ", k_init, 
                    " using filter = ", fltr)
        }
        kc <- as.character(k)

        i <- length(x@init[[kc]])
        ic <- x@init[[kc]][[i]]
        params <- ic$params
        Alpha <- params$Alpha
        Pi <- params$Pi
        Z <- ic$Z
        ds <- ic$Dist

        skeep <- (ds >= fltr)
        skeep[1] <- TRUE

        if (sum(!skeep) > 0){
            Alpha <- Alpha[,skeep,drop=FALSE]
            Pi <- Pi[skeep]
            Pi <- Pi / sum(Pi)
            Z <- Z[,skeep,drop=FALSE]
            Z <- t(apply(Z, 1, function(i){
                         li <- length(i)
                         if (sum(i) == 0){
                             ret <- rep(0, li)
                         } else {
                             ret <- i / sum(i)
                         }
                         return(ret)}))
        }

        params <- list("Alpha" = Alpha, "Pi" = Pi)
        ret <- list("params" = params, "Z" = Z, "Dist" = ds[skeep])

        if (verbose){
            torm <- sum(!skeep)
            message("removed ", torm, 
                    ngettext(torm, " cluster", " clusters"),
                    " with a distance less than ", 
                    fltr, " to the background distribution")
        }

        x@init[[kc]][[i + 1]] <- ret

        if (verbose){
            nclusters <- ncol(ret$params$Alpha) - 1
            message("using a final k of 1 background and ", 
                    nclusters, " cell type clusters")
        }
    }

    return(x)
}

