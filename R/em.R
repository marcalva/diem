
#' Run EM to estimate parameters of mixture model
#' 
#' Estimate the parameters of the mixture model, 
#' and estimate the posterior probability a droplet belongs to each 
#' of the clusters.
#' 
#' @param x An SCE object.
#' @param A Initial alpha parameter.
#' @param P Initial pi parameter.
#' @param Z Initial latent variable Z values.
#' @param model The mixture model to assume. Can be either "DM" for 
#'  a Dirichlet-multinomial or "mltn" for a multinomial.
#' @param alpha_prior Add a non-informative prior by adding a count of
#'  \code{alpha_prior} to all genes in each cluster. Only valid for 
#'  the multinomial model.
#' @param pi_prior Add a non-informative prior by adding a count of 
#'  \code{pi_prior} to the each cluster's membership.
#' @param max_iter Maximum number of iterations for the EM estimation 
#'  of the mixture model.
#' @param eps The delta threshold for when to call convergence for 
#'  the EM estimation of the mixture model. The EM 
#'  stops when delta falls below this value. We define delta as the 
#'  average change in posterior probability. By default this is set to 
#'  1e4, so that the EM converges when less than 1 in 10,000 labels 
#'  change on average.
#' @param psc Pseudocount to add to avoid collapsing likelihood to 0.
#' @param threads Number of threads for parallel execution. Default is 1.
#' @param verbose Verbosity.
#'
#' @return An SCE object.
#'
#' @importFrom Matrix colSums
#' @export
#'
run_em <- function(x, 
                   A = NULL, 
                   P = NULL, 
                   Z = NULL,
                   model = "mltn", 
                   alpha_prior = 0, 
                   pi_prior = 0, 
                   max_iter = 1e3, 
                   eps = 1e-4, 
                   psc = 1e-10, 
                   threads = 1, 
                   verbose = TRUE){

    genes.use <- rownames(x@gene_data)[x@gene_data$exprsd]
    droplets.use <- colnames(x@counts)

    if (length(genes.use) == 0 | length(droplets.use) == 0) stop("specify test set and filter genes before running EM.")
    if (length(x@model) == 0) stop("initialize parameters before running run_em")

    if (is.null(A)){
        A <- x@model$params$Alpha
    }

    if (is.null(P)){
        P <- x@model$params$Pi
    }

    if (is.null(Z)){
        Z <- x@model$Z
    }

    counts <- x@counts[genes.use,]

    N <- length(droplets.use)
    G <- length(genes.use)

    fixed <- which(colnames(x@counts) %in% x@bg_set)

    if (model == "mltn"){
        ret <- mult_em(counts, A, P, Z, fixed, alpha_prior = alpha_prior, 
                       pi_prior = pi_prior, threads = threads, max_iter = max_iter, 
                       eps = eps, psc = psc, display_progress = verbose)
    } else if (model == "DM"){
        ret <- dirmult_em(counts, A, P, Z, fixed, alpha_prior = alpha_prior, 
                          pi_prior = pi_prior, threads = threads, max_iter = max_iter, 
                          eps = eps, psc = psc, display_progress = verbose)
    }

    rownames(ret[[1]]) <- rownames(A)
    colnames(ret[[1]]) <- colnames(A)
    ret[[2]] <- ret[[2]][,1]
    names(ret[[2]]) <- colnames(A)
    rownames(ret[[3]]) <- colnames(counts)
    colnames(ret[[3]]) <- colnames(A)
    rownames(ret[[4]]) <- colnames(counts)
    colnames(ret[[4]]) <- colnames(A)

    x@model$params$Alpha <- ret[[1]]
    x@model$params$Pi <- ret[[2]]
    x@model$Z <- ret[[3]]
    x@model$llk <- ret[[4]][x@test_set,]

    return(x)
}

