
#' Assign clusters to droplets from EM
#'
#' Assign cluster to droplets in the test_data from a DIEM k_init run.
#'
#' @param x An SCE object.
#' @param k_init The k_init run to call droplets from. This must 
#'  specify a single value. If left unspecified, will call droplets only if 
#'  there is one k_init run.
#'
#' @return An SCE object.
#'
#' @export
assign_clusters <- function(x, k_init = NULL){
    k_init <- check_k_init(x, k_init)

    emo <- x@kruns[[k_init]]
    if (length(emo$llk) <= 1){
        stop("run EM to assign clusters")
    }

    Z <- get_z(emo$llk, emo$params$Pi)
    Z <- Z[rownames(x@test_data),]

    clust_max <- apply(Z, 1, which.max)
    clust_prob <- apply(Z, 1, function(i) i[which.max(i)])
    
    x@test_data[,"Cluster"] <- clust_max
    x@test_data[,"ClusterProb"] <- clust_prob

    return(x)
}

#' Call debris clusters using the debris score
#'
#' Specify which clusters from DIEM are debris clusters. This uses the 
#' debris score of droplets from \code{\link{estimate_dbr_score}}. The 
#' average debris score of clusters assigned to the fixed debris cluster 
#' (cluster 1) is calculated. Then, for each cluster, the fraction of 
#' droplets with a debris score above this mean value is calculated. 
#' In addition to the fixed debris cluster, any clusters with a 
#' fraction of debris droplets above \code{thresh} (0.5 by default) 
#' are assigned as debris clusters.
#'
#' @param x An SCE object.
#' @param thresh Clusters with a fraction of debris droplets above this 
#'  value are set to debris clusters.
#' @param top_n Number of differentially expressed genes per cluster 
#'  to use for calculating debris score.
#' @param deb_clust Manually specify clusters which should be set as debris 
#'  clusters.
#' @param k_init The k_init run to call droplets from. This must 
#'  specify a single value. If left unspecified, will call droplets only if 
#'  there is one k_init run.
#' @param verbose verbosity
#'
#' @return An SCE object.
call_debris_clusters <- function(x, 
                                 thresh = 0.25, 
                                 top_n = 20, 
                                 deb_clust = NULL, 
                                 k_init = NULL, 
                                 verbose = TRUE){
    k_init <- check_k_init(x, k_init)

    if (!is.null(deb_clust)){
        clusters <- unique(x@test_data$Cluster)
        deb_clust <- c(1, deb_clust)
        deb_clust <- unique(intersect(clusters, deb_clust))
        x@debris_clusters <- deb_clust
        return(x)
    }

    if (is.null(thresh)){
        return(x)
    }

    emo <- x@kruns[[k_init]]
    if (length(emo$llk) <= 1){
        stop("run EM to assign clusters")
    }

    x <- estimate_dbr_score(x, k_init = k_init)
    sm <- summarize_clusters(x)

    db_thresh <- sm["Debris", "avg_dbr_score"]

    pct_debris_drop <- tapply(x@test_data[,"score.debris"], 
                              x@test_data[,"Cluster"], 
                              function(i){
                                  sum(i >= db_thresh)/length(i)
                              })

    deb_clust <- c(1, names(pct_debris_drop)[pct_debris_drop > thresh])
    deb_clust <- unique(deb_clust)
    x@debris_clusters <- as.integer(deb_clust)

    if (verbose){
        message("found ", length(deb_clust), " debris clusters: ", 
                paste(deb_clust, collapse = " "))
    }

    return(x)
}

#' Call clean droplets after running EM using percent debris
#'
#' Call cells or nuclei from an SCE object. Adds meta data to 
#' the test data slot, which can be retreived
#' using the function \code{\link{droplet_data}}. This adds 
#' column "Call" which takes either "Debris" or "Clean". Droplets
#' with a debris scores less than \code{thresh_score} and at least 
#' \code{min_genes} genes detected pass filtering.
#'
#' @param x An SCE object.
#' @param thresh_score The maximum debris score a droplet must have 
#'  to be classified as a cell/nucleus. This takes values between 
#'  0 and 1.
#' @param min_genes The minimum number of genes a droplet must have 
#'  to be classified as a cell/nucleus.
#' @param k_init The k_init run to call droplets from. This must 
#'  specify a single value. If left unspecified, will call droplets only if 
#'  there is one k_init run.
#' @param verbose verbosity
#'
#' @return An SCE object.
#'
#' @export
call_targets <- function(x, 
                         thresh_score = 0.5, 
                         min_genes = 200, 
                         k_init = NULL, 
                         verbose = TRUE){

    k_init <- check_k_init(x, k_init, return_all = FALSE)

    if (verbose){
        message("calling targets from ", sQuote("k_init"), "=", k_init)
    }

    k_init <- as.character(k_init)

    prob_keep <- x@test_data[, "score.debris"] < thresh_score
    gene_keep <- x@test_data[, "n_genes"] >= min_genes
    
    x@test_data[, "Call"] <- "Debris"
    x@test_data[prob_keep & gene_keep,"Call"] <- "Clean"
    x@test_data[, "Call"] <- as.factor(x@test_data[,"Call"])

    if (verbose){
        n_clean <- sum(x@test_data[,"Call"] == "Clean")
        n_rm <- length(x@test_set) - n_clean
        message("removed ", n_rm, " debris droplets from the test set.")
        message("kept ", n_clean, " clean droplets with debris score < ", 
                thresh_score, " and at least ", min_genes, " genes detected")
    }

    return(x)
}

#' Return IDs of clean droplets
#'
#' Return the barcodes that pass DIEM filtering. The function 
#' \code{\link{call_targets}} must be run beforehand.
#'
#' @param x An SCE object.
#'
#' @return A character vector with the called droplet IDs.
#'
#' @export
get_clean_ids <- function(x){
    if (!"Call" %in% colnames(x@test_data)) 
        stop("call targets before calling get_clean_ids")

    clean <- x@test_data$Call == "Clean"
    ids <- rownames(x@test_data)[clean]
    return(ids)
}

#' Return IDs of removed droplets
#'
#' Return the droplet IDs in the test set that have been removed 
#' by DIEM. 
#' The function \code{\link{call_targets}} must be run beforehand.
#'
#' @param x An SCE object.
#' @param min_genes The minimum number of genes detected for a 
#'  droplet to be output. Useful if a threshold was used to call 
#'  targets, as droplets below would never be considered. By 
#'  default is set at 200.
#'
#' @return A character vector with the called droplet IDs.
#' 
#' @export
get_removed_ids <- function(x, min_genes = 200){
    if (!"Call" %in% colnames(x@test_data)) 
        stop("call targets before calling get_clean_ids")

    if (length(x@test_set) == 0) stop("No test set droplets")
    
    ck <- x@test_data$Call == "Debris"
    gk <- x@test_data$n_genes >= min_genes
    debris <- rownames(x@test_data)[ck & gk]
    removed <- intersect(x@test_set, debris)
    return(removed)
}

