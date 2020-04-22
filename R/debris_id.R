
#' Assign clusters to droplets from EM
#'
#' Assign cluster for each droplet in the test set, along with its 
#' cluster probability. These can be retrieved using the 
#' \code{\link{droplet_data}} function.
#'
#' @param x An SCE object.
#'
#' @return An SCE object.
#'
#' @export
#' @examples
#' \donttest{
#'
#' sce <- assign_clusters(sce, type="test")
#' drop_data <- droplet_data(sce)
#' table(drop_data[,"Cluster"])
#' summary(drop_data[,"ClusterProb"])
#'
#' }
assign_clusters <- function(x){

    emo <- x@model
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
#'  to be classified as a cell/nucleus. This should be a value between 
#'  0 and 1. Set to NULL to turn off filtering by threshold.
#' @param clusters Remove droplets assigned to these clusters.
#' @param min_genes The minimum number of genes a droplet must have 
#'  to be classified as a cell/nucleus.
#' @param droplets A character vector that specifies target droplets 
#'  manually. Only droplets that intersect with test set will be 
#'  used. If given, ignores the other parameteres. Default is NULL.
#' @param verbose verbosity
#'
#' @return An SCE object.
#'
#' @export
#' @examples
#' \donttest{
#'
#' # Default filtering
#' sce <- call_targets(sce)
#' drop_data <- droplet_data(sce)
#' table(drop_data[,"Call"])
#'
#' # More stringent filtering
#' sce <- call_targets(sce, thresh = 0.3)
#' drop_data <- droplet_data(sce)
#' table(drop_data[,"Call"])
#'
#' # Less stringent filtering
#' sce <- call_targets(sce, thresh = 0.7)
#' drop_data <- droplet_data(sce)
#' table(drop_data[,"Call"])
#'
#' }
call_targets <- function(x, 
                         thresh_score = 0.5, 
                         clusters = NULL,
                         min_genes = 200, 
                         droplets = NULL, 
                         verbose = TRUE){

    if (verbose){
        message("calling targets")
    }

    if (!is.null(droplets)){
        d <- intersect(droplets, rownames(x@test_data))
        if (length(d) == 0){
            stop("No droplets intersect with the test set.")
        }
        x@test_data[, "Call"] <- "Debris"
        x@test_data[d, "Call"] <- "Clean"
        x@test_data[, "Call"] <- as.factor(x@test_data[,"Call"])
        return(x)
    }

    if (is.null(thresh_score))
        thresh_score <- Inf

    if (is.null(clusters))
        clusters <- ""

    if (identical(clusters, "debris")){
        sm <- summarize_clusters(x)
        clusters <- sm[ sm[,"Type"] == "Debris", "Cluster"]
    }
    if (identical(clusters,"clean")){
        sm <- summarize_clusters(x)
        clusters <- sm[ sm[,"Type"] == "Clean", "Cluster"]
    }

    cluster_keep <- !(x@test_data[, "Cluster"] %in% clusters)
    prob_keep <- x@test_data[, "score.debris"] < thresh_score
    gene_keep <- x@test_data[, "n_genes"] >= min_genes
    keep <- cluster_keep & prob_keep & gene_keep
    
    x@test_data[, "Call"] <- "Debris"
    x@test_data[keep, "Call"] <- "Clean"
    x@test_data[, "Call"] <- as.factor(x@test_data[,"Call"])

    if (verbose){
        n_clean <- sum(x@test_data[,"Call"] == "Clean")
        n_rm <- length(x@test_set) - n_clean
        if (identical(clusters, "")){
            message("removed ", n_rm, " debris droplets from the test set.")
            message("kept ", n_clean, " clean droplets with debris score < ", 
                    thresh_score, " and at least ", min_genes, " genes detected")
        } else {
            lclust <- length(clusters)
            cmessage <- ngettext(lclust, "cluster ", "clusters ")
            appe <- paste0("from ", cmessage, paste(clusters, collapse = " "))
            message("removed ", n_rm, " debris droplets ", appe, " from the test set.")
            message("kept ", n_clean, " clean droplets with debris score < ", 
                    thresh_score, " and at least ", min_genes, " genes detected")
        }
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

