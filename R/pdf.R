
#' Get log multinomial density of columns in a sparse matrix. 
#'
#' @param X A sparseMatrix, observations in rows, variables in columns
#' @param prob Probability parameter of multinomial. Length must be the
#'  same as the number of columns in X.
#'
#' @importFrom Matrix rowSums
#' @importFrom methods as
#' @importMethodsFrom Matrix %*%
dmultinom_sparse <- function(X, prob){
    if (length(prob) != ncol(X)){
        stop("Length of prob and rows in X must be the same.")
    }
    X <- as(X, "dgCMatrix")
    Xlg <- X
    Xlg@x <- lgamma(Xlg@x + 1) # Get log gamma

    m <- lgamma(rowSums(X) + 1)
    xls <- rowSums(Xlg)
    px <- as.numeric(X %*% log(prob))

    Llks <- m - xls + px
    return(Llks)
}

