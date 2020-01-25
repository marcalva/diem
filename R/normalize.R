
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

