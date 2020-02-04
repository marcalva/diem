
#' @export
get_alpha <- function(x, labs, psc = 1e-16){
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
merge_size <- function(labs, dat, min_size = 30){
    labs <- as.character(labs)
    freq <- table(labs)
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
                 K = 10, 
                 iter.max = 30, 
                 nstart = 30, 
                 min_size = 10, 
                 thresh = 0.05, 
                 seedn = 1, 
                 top_pct = 100, 
                 psc = 1e-16, 
                 verbose = TRUE){ 
    if (verbose) message("initializing parameters")

    # Subset to top percent of test droplets
    sizes <- x@droplet_data[x@test_set, "total_counts"]
    n_keep <- round(length(x@test_set) * (top_pct / 100))
    o <- order(sizes, decreasing=TRUE)
    keep <- x@test_set[o[1:n_keep]]

    # Get PCs
    x <- get_var_genes(x, droplets.use = keep)
    x <- normalize_data(x, droplets.use = keep)
    x <- get_pcs(x, n_pcs = n_pcs)

    labs <- rep(0, ncol(x@counts))
    names(labs) <- colnames(x@counts)
    labs[x@bg_set] <- 1

    pcs <- x@pcs[keep,,drop=FALSE]
    pcs_s <- scale(pcs)
    dpcs <- dist(pcs_s)

    params <- list()
    set.seed(seedn)
    options(warn = -1)
    #hc <- hclust(dist(pcs_s), method = "ward.D2")

    for (k in K){
        km <- kmeans(x = pcs_s, 
                     centers = k, 
                     iter.max = iter.max, 
                     nstart = nstart)
        kclust <- km$cluster
        wss <- km$betweenss

        # Merge small clusters from k-means
        kclust <- merge_size(kclust, pcs_s, min_size)
        kclust <- as.numeric(kclust)
        names(kclust) <- rownames(pcs_s)
        
        # kclust <- cutree(hc, k = k)
        labs_k <- labs
        labs_k[names(kclust)] <- kclust + 1
        labs_k <- factor(labs_k)
        alphas <- get_alpha(x@counts, labs_k, psc = psc)

        # Merge similar cell type clusters
        to_merge <- c()
        pcts <- c()
        cts <- setdiff(colnames(alphas), "1")
        for (ct in cts){
            ctc <- x@counts[,labs_k == ct, drop=FALSE]
            ctc_sizes <- colSums(ctc)
            p1 <- LlkDirMultSparse(ctc, sizes = ctc_sizes, alpha = alphas[,"1"])
            p2 <- LlkDirMultSparse(ctc, sizes = ctc_sizes, alpha = alphas[,ct])
            pg <- sum(p1 > p2) / length(p2)
            pcts[ct] <- pg
            if (pg > thresh) to_merge <- union(to_merge, ct)
        }

        labs_k[labs_k %in% to_merge] <- "1"
        labs_k <- condense_labs(labs_k)
        alphas <- get_alpha(x@counts, labs_k, psc = psc)

        kc <- as.character(k)
        params[[kc]] <- list("labels" = labs_k, 
                             "Alpha" = alphas, 
                             "SSE" = wss)
    }
    options(warn = 0)

    x@init <- params
    return(x)
}

rdirm <- function(n, size, a){
    nt <- length(a) * n
    xi <- rgamma(n = nt, shape = a, scale = 1)
    datf <- matrix(xi, nrow = length(a))
    datf <- apply(datf, 2, function(i) i / sum(i))
    ret <- sapply(1:n, function(i) rmultinom(1, size = size[i], prob = datf[,i]))
    return(ret)
}


