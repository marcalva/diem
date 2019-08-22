
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

#' Log density of Dirichlet-multinomial
#' X is a sample by gene matrix
#' @import Matrix
#' @export
ddm_sparse <- function(X, alpha){
    if (length(alpha) != ncol(X)) stop("length of alpha should be same as number of columns of X.")

    rsX <- Matrix::rowSums(X)
    a0 <- sum(alpha)

    # RET <- lgamma(rsX+1)
    # RET <- RET + lgamma(a0)
    # RET <- RET - lgamma(rsX + a0)
    # RET <- RET - sum(lgamma(alpha))
    # ret <- as.double(0)
    ret <- sum(lgamma(rsX+1))
    ret <- ret + (nrow(X)*lgamma(a0))
    ret <- ret - sum(lgamma(rsX + a0))

    # Traverse by gene
    nz <- nrow(X)
    for (i in seq(ncol(X))){
        ix1 <- X@p[i] + 1
        ix2 <- X@p[i+1] + 1
        n_ix <- ix2-ix1
        ix_seq <- seq(ix1, length.out=n_ix)
        ret <- ret + ((nz - n_ix)*lgamma(alpha[i]))
        ret <- ret + sum(lgamma(X@x[ix_seq] + alpha[i]))
        ret <- ret - sum(lgamma(X@x[ix_seq] + 1))
        ret <- ret - (nz*lgamma(alpha[i]))
    }

    return(as.numeric(ret))
}


#' Log density of Dirichlet-multinomial
#' X is a sample by gene matrix
#' @import Matrix
#' @export
ddm_sparseI <- function(X, alpha){
    # X <- as(X, "dgTMatrix")
    rsX <- Matrix::rowSums(X)
    a0 <- sum(alpha)

    RET <- lgamma(rsX+1)
    RET <- RET + lgamma(a0)
    RET <- RET - lgamma(rsX + a0)
    RET <- RET - sum(lgamma(alpha))

    # Traverse by gene
    nz <- nrow(X)
    for (j in seq(ncol(X))){
        ix1 <- X@p[j] + 1
        ix2 <- X@p[j+1] + 1
        n_ix <- ix2-ix1
        ix_seq <- seq(ix1, length.out=n_ix)
        i <- X@i[ix_seq] + 1
        to_add <- rep(lgamma(alpha[j]), nz)
        to_add[i] <- lgamma(X@x[ix_seq] + alpha[j])
        RET <- RET + to_add
        RET[i] <- RET[i] - lgamma(X@x[ix_seq] + 1)
    }
    return(RET)
}

#' Log density of Dirichlet-multinomial
#' X is a sample by gene matrix
#' @import Matrix
#' @export
ddm_sparseIM <- function(X, alpha){
    # X <- as(X, "dgTMatrix")
    if (nrow(X) != nrow(alpha) || ncol(X) != ncol(alpha)) stop("Dimensions of X and alpha should be the same.")
    rsX <- Matrix::rowSums(X)
    a0 <- rowSums(alpha)

    RET <- lgamma(rsX+1)
    RET <- RET + lgamma(a0)
    RET <- RET - lgamma(rsX + a0)
    RET <- RET - rowSums(lgamma(alpha))

    # Traverse by gene
    nz <- nrow(X)
    for (j in seq(ncol(X))){
        ix1 <- X@p[j] + 1
        ix2 <- X@p[j+1] + 1
        n_ix <- ix2-ix1
        ix_seq <- seq(ix1, length.out=n_ix)
        i <- X@i[ix_seq] + 1
        to_add <- lgamma(alpha[,j])
        to_add[i] <- lgamma(X@x[ix_seq] + alpha[i,j])
        RET <- RET + to_add
        RET[i] <- RET[i] - lgamma(X@x[ix_seq] + 1)
    }
    return(RET)
}

