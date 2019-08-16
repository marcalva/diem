
#' Call clean droplets after running EM
#'
#' Call targets from droplets if the log-likelihood membership 
#' probability (the posterior probability) is higher than \code{pp_thresh}.
#'
#' @param x An SCE object.
#' @param pp_thresh Numeric threshold, where clean droplets must have a 
#'  posterior probability of at least \code{pp_thresh}.
#' @param min_genes Numeric threshold, where clean droplets must have at least 
#'  \code{min_genes} genes detected.
#'
#' @return An SCE object.
#' @export
call_targets <- function(x, pp_thresh=0.95, min_genes=200){

    if (!"CleanProb" %in% colnames(x@droplet_data)) stop("Run DIEM before calling targets")

    calls <- rep("Debris", nrow(x@droplet_data))
    calls[x@droplet_data$CleanProb > pp_thresh & x@droplet_data$n_genes >= min_genes] <- "Clean"
    calls <- as.factor(calls)
    x@droplet_data[,"Call"] <- calls

    return(x)
}

#' @export
score_debris <- function(x, soft=TRUE){

    if (!"CleanProb" %in% colnames(x@droplet_data)) stop("Run DIEM before calling targets")

    droplets.use <- rownames(x@pcs)
    calls <- 1 - x@droplet_data[droplets.use, "CleanProb"]
    scores <- rep(0, length(droplets.use))
    for (i in 1:nrow(x@knn)){
        from <- x@knn[i,"from"]
        scores[from] <- scores[from] + (calls[from] * x@knn[i,"weight"])
    }
    x@droplet_data[droplets.use,"DebrisScore"] <- scores
    return(x)
}

#' Return column IDs of clean droplets
#'
#' @param x An SCE object.
#'
#' @return A character vector with the called droplet IDs.
#' @export
get_clean_ids <- function(x){
    if (!"Call" %in% colnames(x@droplet_data)) stop("Call targets before calling get_clean_ids")
    return(rownames(x@droplet_data)[x@droplet_data$Call == "Clean"])
}

#' Get percent of reads aligning to given genes.
#'
#' @param x An SCE object.
#' @param genes Genes to calculate percentage of in counts.
#' @param name Column name to place in dropl_info.
#'
#' @return An SCE object.
#' @export
#' @examples
#' mm_seur <- get_gene_pct(x=mb_sce, genes="Malat1", name="pct.malat1")
#' mt_genes <- grep(pattern="^mt-", x=rownames(mb_sce@gene_data), value=TRUE, ignore.case=TRUE)
#' mm_seur <- get_gene_pct(x=mb_sce, genes=mt_genes, name="pct.mt")
get_gene_pct <- function(x, genes, name){
    expr <- x@counts[genes,,drop=FALSE]
    if (length(expr) == 0){
        stop("None of genes found in counts.")
    }
    gene_pct <- 100 * Matrix::colSums(expr) / Matrix::colSums(x@counts)
    x@droplet_data[names(gene_pct),name] <- gene_pct
    return(x)
}

#' Merge and re-cluster multiple DIEM objects
#'
#' @param x A list of SCE objects.
#'
#' @return An SCE object.
#' @export
merge_diem <- function(x, use_var_genes=FALSE, eps=1e-8, max_iter=1e3, hcut=0.05, sep="_", name="Merged", verbose=TRUE){
    n_samples <- length(x)

    # Make sure SCE objects have names
    x <- sapply(1:n_samples, function(i){
                sce.this <- x[[i]]
                if (is.null(sce.this@name) || sce.this@name == "SCE"){
                    sce.this@name <- paste0("S", as.character(i))
                }
                return(sce.this)
                          })

    sample_ids <- sapply(x, function(i) i@name)

    # get the union of variable genes
    if (use_var_genes){
        genes.use <- lapply(x, function(i) i@vg)
        genes.use <- Reduce(union, genes.use)
    } else {
        genes.use <- lapply(x, function(i) rownames(i@gene_info)[i@gene_info$exprsd])
        genes.use <- Reduce(intersect, genes.use)
    }

    # get clean droplets
    clean_drop.l <- lapply(1:n_samples, function(i){
                           sce.this <- x[[i]]
                           cds <- get_clean_ids(sce.this)
                           cds <- paste(sample_ids[[i]], cds, sep=sep)
                           return(cds)
                          })

    # get counts
    counts.l <- lapply(1:n_samples, function(i){
                       sce.this <- x[[i]]
                       colnames(sce.this@counts) <- paste(sample_ids[[i]], colnames(sce.this@counts), sep=sep)
                       return(Matrix::t(sce.this@counts[,clean_drop.l[[i]]]))
                          })
    counts.all <- Reduce(rbind, counts.l)
    counts.use <- counts.all[,genes.use]

    # Create SCE
    sce <- create_SCE(Matrix::t(counts.all), name=name)
    sce@dropl_info[,"orig.ident"] <- unlist(sapply(1:n_samples, function(i) rep(sample_ids[[i]], length(clean_drop.l[[i]]))))

    # Get clusters
    clusters.l <- lapply(1:n_samples, function(i){
                         sce.this <- x[[i]]
                         rownames(sce.this@dropl_info) <- paste(sample_ids[[i]], rownames(sce.this@dropl_info), sep=sep)
                         clusts <- sce.this@dropl_info[clean_drop.l[[i]], "mmm_cluster"]
                         clusts <- paste(sample_ids[[i]], clusts, sep=sep)
                         names(clusts) <- clean_drop.l[[i]]
                         return(clusts)
                          })
    clusters.all <- unlist(clusters.l)
    clusters.all <- clusters.all[rownames(counts.use)]

    # Calculate means of all variable genes for each cluster, for all samples
    mu.all <- init_mmm_group(counts.use, clusters.all[rownames(counts.use)])$Mu

    # Merge clusters using hierarchical clustering. Assign new clusters to each sample
    all_groups <- clusters.all
    dist_mu <- as.dist(1 - cor(log1p(1e4*mu.all)))
    mu_clust <- hclust(dist_mu, method="average")
    ct <- cutree(mu_clust, h=hcut)
    # Only merge if ct is > 1
    new_groups <- all_groups
    for (i in names(ct)){
        new_groups[all_groups == i] <- ct[[i]]
    }
    all_groups <- new_groups
    k <- length(unique((all_groups)))
    if (verbose) cat(paste0("Using k=", as.character(k), " groups.\n"))
    if (verbose) print(table(all_groups))
    if (k==1) stop("Can't cluster with k=1 groups.\n")

    # Run unsupervised EM clustering
    mn_params <- init_mmm_group(counts.use, all_groups)

    # Remove droplets with no genes
    no_counts <- Matrix::rowSums(counts.use) == 0
    if (sum(no_counts) > 0){
        counts.use <- counts.use[!no_counts,]
    }

    if (verbose){
        cat(paste0("Running EM\n"))
        cat(paste0("    Classifying ", as.character(nrow(counts.use)), " droplets.\n"))
        cat(paste0("    Using ",as.character(ncol(counts.use)), " genes and ", as.character(nrow(counts.use)), " droplets.\n"))
    }

    # Run EM
    emo <- em(counts.use, k=k, max_iter=max_iter, eps=eps, mn_params=mn_params, verbose=verbose)
    cluster_ids <- paste0("C", as.character(1:ncol(emo$Mu)))
    sce@diem@cluster_ids <- cluster_ids
    sce@diem@emo <- emo
    sce@diem@PP <- t(apply(emo$Z, 1, fraction_log))
    colnames(sce@diem@PP) <- cluster_ids

    sce <- call_targets(sce)
    sce <- fill_dropl_info(sce)

    if (verbose){
        cat("Finished EM\n")
    }

    return(sce)
}


