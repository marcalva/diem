
#' Call clean droplets after running EM
#'
#' Call cells or nuclei from an SCE object. Adds meta data to 
#' droplet data.
#' The posterior probability is defined as 1 - Prob(Debris), where 
#' Prob(Debris) is the posterior probability that a droplet is 
#' debris.
#'
#' @param x An SCE object.
#' @param pp_thresh The minimum posterior probability of a droplet to be 
#'  classified as a cell/nucleus.
#' @param min_genes The minimum number of genes a droplet must have 
#'  to be classified as a cell/nucleus.
#' @param verbose verbosity
#'
#' @return An SCE object.
#' @export
call_targets <- function(x, 
                         pp_thresh = 0.95, 
                         min_genes = 200, 
                         k_init = NULL, 
                         verbose = TRUE){

    if (length(x@emo) == 0) stop("Run ", sQuote("run_em"), " before calling targets")

    if (is.null(k_init)){
        if (length(x@emo) > 1){
            stop(sQuote("k_init"), " must be specified to a value from: ", 
                 paste(names(x@emo), collapse = " "))
        } else {
            k_init <- names(x@emo)[1]
        }
    } else {
        if (length(k_init) != 1){
            stop(sQuote("k_init"), " must be a single value")
        }
        if ( ! k_init %in% names(x@emo) ){
            stop(sQuote("k_init"), " must an initialized k value taken from: ", 
                 paste(names(x@emo), collapse = " "))
        }
    }
    if (verbose){
        message("calling targets from ", sQuote("k_init"), "=", k_init)
    }

    k_init <- as.character(k_init)

    emo <- x@emo[[k_init]]
    probs <- 1 - emo$Z[,1]
    calls <- rep("Debris", nrow(x@droplet_data))
    calls[probs >= pp_thresh & x@droplet_data$n_genes >= min_genes] <- "Clean"
    calls <- as.factor(calls)
    x@droplet_data[,"Call"] <- calls

    # Add debris log odds
    llk <- emo$llk[x@test_set,,drop=FALSE]
    x@droplet_data[x@test_set, "DebrisLlk"] <- llk[,1]
    Pi <- emo$params$Pi

    llks_pi <- t(apply(llk, 1, function(j) j + log(Pi)))
    debris_prob <- llks_pi[,1]
    clean_prob <- apply(llks_pi[,2:ncol(llks_pi),drop=FALSE], 1, sum_log)
    x@droplet_data[x@test_set, "DebrisLogOdds"] <- debris_prob - clean_prob

    x@droplet_data[,"Cluster"] <- emo$Cluster
    x@droplet_data[,"ClusterProb"] <- emo$ClusterProb
    x@droplet_data[,"DebrisProb"] <- emo$Z[,1]

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
#' @export
get_clean_ids <- function(x){
    if (!"Call" %in% colnames(x@droplet_data)) 
        stop("Call targets before calling get_clean_ids")

    clean <- x@droplet_data$Call == "Clean"
    ids <- rownames(x@droplet_data)[clean]
    return(ids)
}

#' Return IDs of removed droplets
#'
#' Return the droplet IDs in the test set that have been removed 
#' by DIEM. 
#' The function \code{\link{call_targets}} must be run beforehand.
#'
#' @param x An SCE object.
#'
#' @return A character vector with the called droplet IDs.
#' @export
get_removed_ids <- function(x){
    if (!"Call" %in% colnames(x@droplet_data)) 
        stop("Call targets before calling get_clean_ids")

    if (length(x@test_set) == 0) stop("No test set droplets")
    debris <- rownames(x@droplet_data)[x@droplet_data$Call == "Debris"]
    removed <- intersect(x@test_set, debris)
    return(removed)
}

#' Get percent of reads align to given gene(s)
#'
#' Add a data column for each droplet giving the percentage of raw reads/UMIs 
#' that align to genes given in \code{genes}. The column name is specified by 
#' \code{name}.
#'
#' @param x An SCE object.
#' @param genes Genes to calculate percentage of in counts.
#' @param name Column name to place in dropl_info.
#'
#' @return An SCE object.
#' @importFrom Matrix colSums
#' @export
#' @examples
#' # Add MT%
#' mt_genes <- grep(pattern="^mt-", x=rownames(mb_small@gene_data), ignore.case=TRUE, value=TRUE)
#' mb_small <- get_gene_pct(x = mb_small, genes=mt_genes, name="pct.mt")
#' # Add MALAT1
#' genes <- grep(pattern="^malat1$", x=rownames(mb_small@gene_data), ignore.case=TRUE, value=TRUE)
#' mb_small <- get_gene_pct(x = mb_small, genes=genes, name="MALAT1")
get_gene_pct <- function(x, genes, name){
    expr <- x@counts[genes,,drop=FALSE]
    if (length(expr) == 0){
        stop("None of the given genes found.")
    }
    gene_pct <- 100 * colSums(expr) / colSums(x@counts)
    x@droplet_data[names(gene_pct),name] <- gene_pct
    return(x)
}

