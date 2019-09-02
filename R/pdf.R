
#' Get log multinomial density of columns in a sparse matrix. 
#'
#' @param X A sparseMatrix that is a sample by feature matrix of counts.
#' @param prob Probability parameter of multinomial. Must be same size as 
#'  number of columns in X.
#'
#' @import Matrix
#' @export
dmultinom_sparse <- function(X, prob){
    if (length(prob) != ncol(X)){
        stop("Length of prob and rows in X must be the same.")
    }
    X <- as(X, "dgCMatrix")
    Xlg <- X
    Xlg@x <- lgamma(Xlg@x + 1) # Get log gamma

    m <- lgamma(Matrix::rowSums(X) + 1)
    xls <- Matrix::rowSums(Xlg)
    px <- as.numeric(X %*% log(prob))

    Llks <- m - xls + px
    return(Llks)
}

