
#' @export
get_alpha <- function(x, labs, psc = 1e-16){
    if (length(labs) != ncol(x)){
        stop("Number of columns in x must match length of labs")
    }
    labs <- factor(labs)
    labsu <- unclass(labs)
    clusts <- levels(labs)
    alphas <- sapply(clusts, function(j){
                     clust_mem <- labs == as.character(j)
                     rowSums(x[,clust_mem,drop=FALSE]) + psc
                           })
    z = lapply(1:ncol(x), function(i) {
               ret = rep(0, length(clusts))
               ret[labsu[i]] = 1
               return(ret)})
    z = do.call(rbind, z)
    alphas = dm_loo(t(x), alphas, z)
    colnames(alphas) <- clusts
    return(alphas)
}

get_pi <- function(counts, Alpha=NULL, llks = NULL, sizes = NULL){
    if (is.null(llks)){
        llks <- get_llk(counts, Alpha, sizes)
    }
    llksf <- t(apply(llks, 1, fraction_log))
    pis <- colSums(llksf)
    pis <- pis / sum(pis)
    return(pis)
}

init_z <- function(labs){
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
#' 
#' @export
merge_size <- function(labs, dat, min_size = 10, verbose = TRUE){
    labs <- as.character(labs)
    freq <- table(labs)
    nm <- 0
    while(min(freq) < min_size){
        clusts <- names(freq)
        means <- sapply(clusts, function(i) colMeans(dat[labs == i,,drop=FALSE]))
        dists <- as.matrix(dist(t(means), upper = TRUE, diag = TRUE))
        small_clust <- names(freq)[which.min(freq)]
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

#' Initialize clusters using hierarchical clustering
#' k is a positive integer or vector of positive integers
#' @importFrom Matrix Diagonal colSums rowMeans rowSums
#' @export
init <- function(x, 
                 n_pcs = 30, 
                 K = 15, 
                 iter.max = 30, 
                 nstart = 30, 
                 min_size = 10, 
                 thresh = 0.1, 
                 seedn = 1, 
                 psc = 1e-16, 
                 verbose = TRUE){ 
    if (verbose) message("initializing parameters")

    genes.use <- rownames(x@gene_data)[x@gene_data$exprsd]

    counts <- x@counts[genes.use,]
    sizes <- colSums(counts)

    # Get PCs
    x <- get_var_genes(x)
    x <- normalize_data(x)
    x <- get_pcs(x, n_pcs = n_pcs)
    pcs_s <- scale(x@pcs)

    labs <- rep(0, ncol(counts))
    names(labs) <- colnames(counts)
    labs[x@bg_set] <- 1

    params <- list()
    set.seed(seedn)

    for (k in K){
        options(warn = -1)
        km <- kmeans(x = pcs_s, 
                     centers = k, 
                     iter.max = iter.max, 
                     nstart = nstart)
        options(warn = 0)
        kclust <- km$cluster

        # Merge small clusters from k-means
        kclust <- merge_size(kclust, pcs_s, min_size, verbose = verbose)
        kclust <- as.numeric(kclust)
        names(kclust) <- rownames(pcs_s)
        
        labs_k <- labs
        labs_k[names(kclust)] <- kclust + 1
        alphas <- get_alpha(counts, labs_k, psc = psc)
        Z <- init_z(labs_k)

        # Merge similar cell type clusters
        # alphas <- rm_debris_clust(counts[,x@test_set], alphas, sizes = sizes[x@test_set], thresh = thresh)
        ret <- rm_debris_clust(counts, alphas, Z, NULL, thresh)
        alphas <- ret$Alpha
        pis <- get_pi(counts, alphas, sizes = sizes)

        kc <- as.character(k)
        params[[kc]] <- list("Alpha" = alphas, "Pi" = pis, "Z" = Z) 
    }

    x@init <- params
    return(x)
}

