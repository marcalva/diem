
#' Divide elements of a column by the column's sum in a sparse matrix
#'
#' @param x A Sparse Matrix.
#'
#' @return The Sparse Matrix x with columns summing to 1.
#' @importFrom Matrix Diagonal colSums
divide_by_colsum <- function(x){
    if (!inherits(x, what="Matrix")) stop("Argument must be of class Matrix.")
    cs <- colSums(x = x)
    d <- Diagonal(x = 1/cs)
    x <- x %*% d
    return(x)
}

#' Normalize counts of a sparse matrix
#'
#' @param counts A Sparse Matrix.
#' @param scale_factor A numeric scaling factor to 
#'  multiply values after dividing counts matrix by total counts.
#' @param logt A logical indicating whether to log transform after 
#'  normalizing for columns to sum to 1 and multiplying by scaling factor.
#'
#' @return Sparse Matrix
norm_counts <- function(counts, scale_factor = 1e4, logt = TRUE){
    if (!inherits(counts, what="Matrix")) stop("Argument must be of class Matrix.")
    counts <- divide_by_colsum(counts)
    counts <- counts * scale_factor
    if (logt) counts <- log1p(counts)
    return(counts)
}

#' Normalize raw counts.
#'
#' Normalization of raw counts in an SCE object. Normalization is performed 
#' for the initialiation of the EM. The initialization involves 
#' clustering the high-count high-confidence droplets to approximately 
#' identify the cell types present in the data. To best identify 
#' cell types, clustering is done on the normalized counts.
#' 
#' By default this only normalizes droplets in the cluster set, as 
#' only these droplets that are used for the intialization and are 
#' the only droplets that require count normalization. Unless specified 
#' with \code{genes.use}, only variable genes are included in the 
#' normalization. 
#' The data is normalized by dividing counts by the total counts per droplet. 
#' Then, the counts are multiplied by a scaling factor, given by 
#' \code{sf} (the median of total counts by default). Finally, the data is 
#' log transformed after adding a constant value of 1.
#'
#' @param x An SCE object.
#' @param droplets.use A character vector of droplet IDs to subset the counts data. 
#'  Normalization will only be run on these droplets.
#' @param genes.use A character vector of gene names to subset the counts data. 
#'  Normalization will only be run for these genes.
#' @param use_var A logical indicating whether to subset the data to include
#'  only variable genes. This overrides \code{genes.use}.
#'  The default is TRUE as it may better identify cell types.
#' @param sf Either a numeric scaling factor to multiply counts after 
#'  division by column sums, or "median" indicating to multiply by the median 
#'  number of total read/UMI counts in droplets (default).
#' @param logt A logical specifying whether to log(x+1) transform counts after 
#'  size normalization. Default is TRUE.
#'
#' @return An SCE object
#' @importFrom Matrix rowSums colSums
#' @importFrom stats median
normalize_data <- function(x, 
                           droplets.use = NULL, 
                           genes.use = NULL, 
                           use_var = TRUE, 
                           sf = "median", 
                           logt = TRUE){
    if (is.null(droplets.use)) droplets.use <- x@cluster_set
    if (length(droplets.use) == 0) stop("0 droplets specified in 'droplets.use'.")
    if (is.null(genes.use)){
        if (sum(x@gene_data[,"exprsd"]) == 0){
            genes.use <- rownames(x@counts)[rowSums(x@counts[,droplets.use]) > 0]
        } else {
            genes.use <- rownames(x@gene_data)[x@gene_data[,"exprsd"]]
        }
    }

    if (use_var){ 
        if (length(x@vg) == 0) stop("0 variable genes found.")
        genes.use <- x@vg
    }

    expr <- x@counts[genes.use, droplets.use]
    if (sf == "median") sf <- median(colSums(expr))
    x@norm <- norm_counts(expr, scale_factor=sf, logt=logt)

    return(x)
}

