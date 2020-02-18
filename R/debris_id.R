
#' Call clean droplets after running EM
#'
#' Call cells or nuclei from an SCE object. Adds meta data to 
#' the droplet data slot. This data frame can be retreived
#' using the function \code{\link{droplet_data}}. This adds 
#' column "Call" which takes either "Debris" or "Clean", "DebrisLlk" 
#' which is the log likelihood of the debris distribution, 
#' "DebrisLogOdds" which is the log odds of the debris, "DebrisProb" 
#' which is the probability of the debris, "ClusterPob" which is 
#' the probability of the cell type assignment, and "Cluster" which 
#' is an integer that specifies debris (1) or cell type (>1). The 
#' function calls droplets that have a probability of being a cell 
#' type greater of at least \code{pp_thresh} and at least 
#' \code{min_genes} detected as "Clean." 
#'
#' @param x An SCE object.
#' @param pp_thresh The minimum posterior probability of a droplet to be 
#'  classified as a cell/nucleus.
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
                         pp_thresh = 0.95, 
                         min_genes = 200, 
                         k_init = NULL, 
                         verbose = TRUE){

    k_init <- check_k_init(x, k_init, return_all = FALSE)

    if (verbose){
        message("calling targets from ", sQuote("k_init"), "=", k_init)
    }

    k_init <- as.character(k_init)

    emo <- x@kruns[[k_init]]
    Z <- emo$Z
    if (length(emo$llk) <= 1){
        stop("run EM to estimate parameters and likelihoods before calling targets")
    }

    clust_max <- apply(Z, 1, which.max)
    clust_prob <- apply(Z, 1, function(i) i[which.max(i)])

    # Place calls
    calls <- rep("Debris", nrow(x@droplet_data))
    prob_keep <- (1 - Z[,1]) >= pp_thresh
    gene_keep <- x@droplet_data$n_genes >= min_genes
    calls[prob_keep & gene_keep] <- "Clean"
    x@droplet_data[,"Call"] <- as.factor(calls)

    # Add debris log odds
    llk <- emo$llk[x@test_set,,drop=FALSE]
    Pi <- emo$params$Pi
    llk_pi <- t(apply(llk, 1, function(j) j + log(Pi)))
    debris_prob <- llk_pi[,1]
    clean_prob <- apply(llk_pi[,2:ncol(llk_pi),drop=FALSE], 1, sum_log)
    x@droplet_data[x@test_set, "DebrisLlk"] <- debris_prob
    x@droplet_data[x@test_set, "DebrisLogOdds"] <- debris_prob - clean_prob

    x@droplet_data[,"DebrisProb"] <- Z[,1]
    x@droplet_data[,"ClusterProb"] <- clust_prob
    x@droplet_data[,"Cluster"] <- clust_max

    if (verbose){
        n_clean <- sum(x@droplet_data[,"Call"] == "Clean")
        n_rm <- length(x@test_set) - n_clean
        message(paste0("removed ", 
                       as.character(n_rm), 
                       " debris droplets from the test set."))
        message(paste0("kept ", 
                       as.character(n_clean), 
                       " clean droplets with probability >= ", 
                       pp_thresh, " and at least ", min_genes, 
                       " genes detected"))
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
#' @param min_genes The minimum number of genes detected for a 
#'  droplet to be output. Useful if a threshold was used to call 
#'  targets, as droplets below would never be considered. By 
#'  default is set at 200.
#'
#' @return A character vector with the called droplet IDs.
#' 
#' @export
get_removed_ids <- function(x, min_genes = 200){
    if (!"Call" %in% colnames(x@droplet_data)) 
        stop("Call targets before calling get_clean_ids")

    if (length(x@test_set) == 0) stop("No test set droplets")
    
    ck <- x@droplet_data$Call == "Debris"
    gk <- x@droplet_data$n_genes >= min_genes
    debris <- rownames(x@droplet_data)[ck & gk]
    removed <- intersect(x@test_set, debris)
    return(removed)
}

#' Get percent of reads aligned to given gene(s)
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
#' 
#' @importFrom Matrix colSums
#' 
#' @export
#' 
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

