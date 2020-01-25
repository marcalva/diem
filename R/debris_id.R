
#' Get betas from OLS
#' @importFrom stats lm
get_betas_lm <- function(y, X){
    ret <- lm(y ~ 0 + X)
    return(ret$coefficients)
}

#' Get betas from NNLS
#' @importFrom nnls nnls
get_betas_nnls <- function(y, X){
    ret <- nnls(X, y)
    return(ret$x)
}

#' Remove debris reads from droplets
#'
#' This function takes the raw counts matrix and removes the counts that 
#' are estimated to originate from debris. For each droplet, the percent 
#' of debris and cell type are estimated by non-negative least squares. 
#' Then, a debris droplet is calculated with same total read count as 
#' the droplet and multiplied by the debris proportion. The non-integer 
#' numbers are rounded down. Finally, this debris vector is subtracted 
#' from the droplet counts, keeping reads counts to non-negative.
#' 
#' @return An SCE object
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom methods as
#' @export
remove_debris <- function(x, verbose = FALSE){   
    if (verbose) message("Filtering reads")
    dg <- x@assignments == "Debris"
    counts <- x@counts
    total_counts <- colSums(counts)
    drops <- x@test_set
    genes <- rownames(subset(x@gene_data, exprsd == TRUE))
    genes_all <- seq_along(rownames(x@counts)); names(genes_all) <- rownames(x@counts); 
    gene_probs <- x@vemo$params$Beta
    gene_probs <- sweep(gene_probs, 2, colSums(gene_probs), "/")
    coefs <- matrix(nrow = length(drops), ncol = ncol(gene_probs))
    colnames(coefs) <- colnames(gene_probs); rownames(coefs) <- drops
    counts <- counts[genes,drops]
    total_counts <- total_counts[drops]
    gene_probs <- gene_probs[genes,,drop=FALSE]
    countsp <- divide_by_colsum(counts)

    coefs <- sapply(drops, function(i) get_betas_nnls(countsp[,i], gene_probs))
    coefs <- t(coefs)
    dc <- lapply(1:length(drops), function(j){
                  to_rm <- (gene_probs[,dg,drop=FALSE] * total_counts[j]) %*% t(coefs[j,dg,drop=FALSE])
                  to_rm <- floor(to_rm)
                  genes_gt0 <- rownames(to_rm)[to_rm[,1] > 0]
                  i <- genes_all[genes_gt0]
                  if (length(i) == 0) return(NULL)
                  datf <- data.frame("i" = i, "j" = j, "x" = to_rm[genes_gt0, 1])
                  return(datf)})
    dc <- do.call(rbind, dc)
    dc <- Matrix::sparseMatrix(i = dc[,"i"], 
                               j = dc[,"j"], 
                               x = dc[,"x"], 
                               dims = c(nrow(x@counts), length(drops)) )
    colnames(dc) <- drops
    cf <- x@counts[,drops,drop=FALSE] - dc
    cf[cf < 0] <- 0
    x@counts_filt <- as(cf, "CsparseMatrix")
    x@coefs <- coefs
    if (verbose) message("Done")
    return(x)
}

#' Call clean droplets after running EM
#'
#' Call cells or nuclei from an SCE object. EM must be run beforehand.
#' The posterior probability is defined as 1 - Prob(Debris), where 
#' Prob(Debris) is the posterior probability that a droplet is 
#' debris.
#'
#' @param x An SCE object.
#' @param pp_thresh The minimum posterior probability of a droplet to be 
#'  classified as a cell/nucleus.
#' @param min_genes The minimum number of genes a droplet must have 
#'  to be classified as a cell/nucleus.
#'
#' @return An SCE object.
#' @export
call_targets <- function(x, pp_thresh = 0.95, min_genes = 200){

    if (!"CleanProb" %in% colnames(x@droplet_data)) stop("Run DIEM before calling targets")



    calls <- rep("Debris", nrow(x@droplet_data))
    calls[x@droplet_data$CleanProb >= pp_thresh & x@droplet_data$n_genes >= min_genes] <- "Clean"
    calls <- as.factor(calls)
    x@droplet_data[,"Call"] <- calls

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

