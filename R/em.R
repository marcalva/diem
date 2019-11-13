
#' Initialize EM parameters given group assignments
#'
#' @param counts an observation by variable matrix of non-negative
#'  integer counts.
#' @param groups A character vector or factor, the same length 
#'  as the number of columns in counts. The assignment of the 
#'  index in groups corresponds to the group assignment of the 
#'  column index in counts. Thus, the ordering of groups must 
#'  be consistent with the column ordering in counts.
#' @param psc Pseudocount to add to the resulting mean for 
#'  each multinomial mean parameter.
#'
#' @return List
#' \describe{
#'  \item{mu0}{A observation by k matrix of multinomial means.}
#'  \item{mc0}{A numeric vector of length k of micture coefficients.}
#' }
#' @importFrom Matrix colSums
init_param <- function(counts, groups, psc=1e-4){

    if (length(groups) != nrow(counts)){
        stop("Length of groups must be same as number of rows in counts.")
    }

    groups <- droplevels(as.factor(groups))
    ks <- levels(groups)
    k <- nlevels(groups)

    mu0 <- matrix(nrow=ncol(counts), ncol=k)
    rownames(mu0) <- colnames(counts)
    colnames(mu0) <- ks

    for (clust in ks){
        k_barcodes <- rownames(counts)[groups == clust]
        k_mu <- colSums(counts[k_barcodes,,drop=FALSE]) + psc
        mu0[,clust] <- k_mu/sum(k_mu)
    }

    mc0 <- table(groups)[ks]
    mc0 <- mc0/sum(mc0)
    ret <- list(Mu=mu0, Mc=mc0)
    return(ret)
}

#' fraction of logs
#'
#' @param x numeric vector
fraction_log <- function(x){
    x_c = x - max(x)
    x_c = exp(x_c)
    frac = x_c / sum(x_c);
    return(frac);
}

#' sum of logs
#'
#' @param x numeric vector
sum_log <- function(x){
    max_x <- max(x)
    x_c = x - max_x
    x_sum <- log(sum(exp(x_c))) + max_x
    return(x_sum);
}

#' Compute density of multinomial mixture for a matrix
#'
#' Given an n x p matrix of count data \code{x}, as well as a 
#' p x k matrix of probabilities \code{p}, return the probabilities 
#' of each sample for the k groups in a n x k matrix. If \code{x} 
#' has any non-integer values, they are rounded to the nearest integer.
#'
#' @param x A sparseMatrix that is a sample by feature matrix of counts.
#' @param p Matrix with probability parameters of multinomials for each group.
#'  Number of rows must be same as number of columns in X
#' @param mc Numeric vector with mixture coefficients. Length must be same as 
#'  number of columns in \code{p}.
#' @param labels Numeric vector of same length as number of rows in x. Fixes 
#'  the group probabilities of the integer in this vector element to 1. In 
#'  other words, the latent variable for these samples are treated as known.
#'
#' @return a numeric matrix with n samples by k groups containing log multinomial 
#'  densities under the given parameters of the MMM.
dmmn <- function(x, p, mc, labels=NULL){
    if (ncol(x) != nrow(p)) stop("Number of columns in x must be the same as the number of rows in p.")
    if (length(mc) != ncol(p)) stop("Length of mc must be the same as the number of columns in p.")
    if ( any(round(colSums(p), 1) != 1) ) stop("Probability columns in p must sum to 1.")
    if (round(sum(mc), 1) != 1) stop("Mixture probabilities in mc must sum to 1.")

    n <- nrow(x)
    k <- ncol(p)

    if (is.null(labels)) labels <- numeric(n)
    if (length(labels) != n){
        stop("Size of labels must be the same as the number of rows in x.")
    }
    if ( any(is.na(labels)) ){
        stop("Labels must not contain any NA values.")
    }
    ks <- sort(unique(labels))
    if (max(ks) > k){
        stop("Labels must have a max integer value no more than the number of columns in p.")
    }
    ks <- ks[ks != 0]

    Llks <- matrix(nrow=n, ncol=k)
    colnames(Llks) <- colnames(p)
    rownames(Llks) <- rownames(x)
    for (ki in 1:k){
        Llks[,ki] <- dmultinom_sparse(X=x, prob=p[,ki])
        mc_v <- rep(x=log(mc[ki]), times=n) # Add log of mixture coefficient
        mc_v[(labels != 0) & (labels != ki)] <- -Inf
        Llks[,ki] <- Llks[,ki] + mc_v
    }

    if (any(apply(Llks, 1, function(i) all(is.infinite(i))))){
        warning("For at least one sample, probabilities returned 0 for all k groups. ", 
                "Setting likelihood to 1 for each group.")
        Llks[apply(Llks, 1, function(i) all(is.infinite(i))), ] <- c(0,0)
    }

    return(Llks)
}

#' Computed expected log likelihood of multinomial mixture
#'
#' Given an n x p matrix of count data \code{x}, as well as a 
#' p x k matrix of probabilities \code{p}, return the responsibilities 
#' of each sample for the k groups in a n x k matrix. If \code{x} 
#' has any non-integer values, they are rounded to the nearest integer.
#'
#' @param x sparseMatrix. A sample by feature matrix of counts.
#' @param p numeric matrix. A feature by group matrix of probabilities.
#' @param mc numeric. Mixture coefficients
#' @param Llks If log likelihood was computed for the mixture, provide the 
#'  matrix here.
#' @param labels Numeric vector of same length as number of rows in x. Fixes 
#'  the group probabilities of the integer in this vector element to 1. In 
#'  other words, the latent variable for these samples are treated as known.
#'
#' @return a numeric matrix with n samples by k groups.
e_step_mn <- function(x, p, mc, Llks=NULL, labels=NULL){
    if (any(p == 0)) 
        stop("Probability of a genes(s) in at least 1 group is 0,",  "
             which would collapse likelihood to 0 and prevent classification. ", 
             "Increase cpm threshold to prevent this.")

    if (is.null(Llks)) Llks <- dmmn(x=x, p=p, mc=mc, labels=labels)

    r <- t(apply(Llks, 1, fraction_log))
    if (any(is.infinite(r))) stop("One or more responsibilities are infinite, probably from dividing by 0.")
    return(r)
}

#' Maximize multinomial parameters in EM
#'
#' Given an n x p matrix of count data, as well as a 
#' n x k matrix of responsibilities, return the maximum 
#' likelihood estimate (MLE) of p for each k group in 
#' a p x k matrix.
#'
#' @param x An n x p count matrix.
#' @param r An n x k membership prob. matrix.
#' @param psc Pseudocount.
#' 
#' @return A list with 
#' \describe{
#'  \item{Mu}{An n x k matrix of multinomial means.}
#'  \item{Mc}{A k-length numeric vector of mixture coefficients.}
#' }
#' 
#' @importFrom Matrix t
#' @importMethodsFrom Matrix %*%
m_step_mn <- function(x, r, psc=1e-4){
    n <- nrow(x)
    k <- ncol(r)
    wm <- as.matrix(t(x) %*% r) + psc
    Mu <- sweep(wm, 2, colSums(wm), "/")
    mc <- colSums(r)
    mc <- mc/sum(mc)
    return(list(Mu=Mu, Mc=mc))
}

#' EM function
#' 
#' Run EM for a multinomial mixture model on a sample x feature matrix. 
#' Take a matrix \code{counts} and classify the samples in each row into 
#' one of \code{k} clusters. The initial parameters of the multinomial 
#' mixture model must be given as a list in the parameter \code{mn_params}.
#'
#' @param counts observation by variable matrix of non-negative 
#'  integer counts.
#' @param k Number of clusters.
#' @param mn_params A list containing 
#'  \describe{
#'    \item{Mu}{A variable by k matrix containing the means of the 
#'    of the k multinomial distributions.}
#'    \item{Mc}{A numeric vector of mixture coefficients.}
#' }
#' @param max_iter A numeric value indivating the maximum 
#'  number of iterations.
#' @param eps The epsilon value, which is the convergence 
#'  threshold of the percent change in the log likelihood.
#' @param psc Pseudocount to add to the multinomial mean 
#'  parameter to avoid the likelihood collapsing to 0.
#' @param labels Numeric vector of same length as number of 
#'  observations in counts. Fixes the group probabilities 
#'  of the integer in this vector element to 1. In 
#'  other words, the latent variable for these samples 
#'  are treated as known.
#' @param verbose verbosity.
#'
#' @return A list with
#' \describe{
#'   \item{Z}{An observation by cluster matrix of log 
#'   log likelihoods. Each element is the log likelihood 
#'   of that data point under the the k multinomial.}
#'   \item{Mu}{A variable by cluster matrix of 
#'   multinomial parameters.}
#'   \item{Mc}{A numeric vector of mixture coefficients.}
#'   \item{loglk}{Data log likelihood.}
#'   \item{converged}{logical indicating whether the 
#'   EM algorithm converged (TRUE) or reached the 
#'   maximum number of iterations.}
#' }
em <- function(counts, 
               k, 
               mn_params, 
               max_iter = 1e2, 
               eps = 1e-8, 
               psc = 1e-4, 
               labels = NULL, 
               verbose = TRUE){
    # Initialize parameters
    loglk <- -Inf
    loglks <- c(loglk)
    iter <- 1

    if (ncol(counts) != nrow(mn_params$Mu)) stop("Number of columns in counts must be the same as the number of features in mn_params.")
    if (nrow(counts) != length(labels)) stop("Number of rows in counts must be the same as the length of labels.")

    Llks <- dmmn(x=counts, p=mn_params$Mu, mc=mn_params$Mc, labels=labels)
    while (iter <= max_iter){
        # E Step
        r <- e_step_mn(x=counts, p=mn_params$Mu, mc=mn_params$Mc, Llks=Llks, labels=labels)
        # M Step
        mn_params <- m_step_mn(x=counts, r=r, psc=psc)
        # Evaluate likelihood
        Llks <- dmmn(x=counts, p=mn_params$Mu, mc=mn_params$Mc, labels=labels)
        loglk <- loglks[[iter]] <- sum(apply(Llks, 1, sum_log))
        # Save output
        emo <- list(Z=Llks, Mu=mn_params$Mu, Mc=mn_params$Mc, loglk=loglk, converged=FALSE)

        if (verbose){
            message("iteration ", iter, "; llk ", loglk)
        }

        if (iter > 1){
            dloglk <- (loglks[[iter]] - loglks[[iter-1]])/abs(loglks[[iter-1]])
            if (dloglk < 0){
                warning("likelihood decreased")
            }
            if (dloglk < eps){
                emo$converged=TRUE
                if (verbose) message("converged!")
                break
            }
        }
        iter <- iter + 1
    }
    if (!emo$converged) warning("EM did not converge")

    emo$PP <- t(apply(emo$Z, 1, fraction_log))

    return(emo)
}

#' Run EM on counts to estimate multinomial mixture model
#' 
#' Run expectation maximization (EM) to estimate the parameters of the 
#' multinomial mixture model. This function takes an SCE 
#' object as input, and returns an SCE object with the EM output. 
#' The number of clusters and their initialized multinomial means 
#' are taken from the initial clustering assignments calculated with 
#' \code{\link{initialize_clusters}}. The EM algorithm is run  
#' by repeatedly updating the membership probabilities and then 
#' estimating the MLE parameters of the multinomial mixture model.
#' The algorithm converges when the percent change in 
#' the log likihood is less than \code{eps}. If the 
#' algorithm doesn't converge by \code{max_iter}, it breaks off.
#' The posterior probability is the calculated by taking the 
#' sum of the likelihood fractions across the cell types.
#'
#' @param x An SCE object.
#' @param eps Numeric threshold. The EM algorithm converges when the 
#'  percent change in log likihood is less than \code{eps}.
#' @param max_iter The maximum number of iterations allowed to run.
#' @param psc Pseudocount to add to multinomial parameters 
#'  to avoid collapsing likelihood to 0.
#' @param verbose Logical indicating verbosity.
#'
#' @return An SCE object.
#' @importFrom igraph gsize
#' @importMethodsFrom Matrix %*%
#' @export
run_em <- function(x, 
                   eps = 1e-8, 
                   max_iter = 1e3, 
                   psc = 1e-4, 
                   verbose = TRUE){

    genes.use <- rownames(x@gene_data)[x@gene_data$exprsd]
    droplets.use <- c(x@test_set, x@bg_set)

    if (length(genes.use) == 0 | length(droplets.use) == 0) stop("Specify test set and filter genes before running EM.")
    if (length(x@ic$clusters) == 0) stop("No clusters found, run initialize_clusters before running EM.")

    counts <- t(x@counts[genes.use,droplets.use])

    # Initialize means of clusters
    ic <- x@ic$clusters
    mn_params <- init_param(counts[names(ic),], ic, psc=psc)

    # Fix labels
    labels <- rep(0, length(droplets.use))
    names(labels) <- droplets.use
    labels[x@bg_set] <- 1
    bg_clusts <- names(x@ic$assignments)[x@ic$assignments == "Debris"]

    k <- nlevels(ic)

    if (verbose){
        message("running EM")
    }

    # Run EM
    emo <- em(counts, 
              k = k, 
              max_iter = max_iter, 
              eps = eps, 
              labels = labels, 
              psc = psc, 
              mn_params = mn_params, 
              verbose = verbose)
    x@emo <- emo

    # Naive Bayes assignment
    a <- x@ic$assignments
    PP_summary <- sapply(levels(a), function(i) rowSums(emo$PP[,names(a)[a == i],drop=FALSE]))
    clust_prob <- apply(emo$PP, 1, function(i) i[which.max(i)])
    clust_max <- apply(emo$PP, 1, function(i) colnames(emo$PP)[which.max(i)])
    clust_max <- as.factor(clust_max)

    # Fill in droplet data
    x@droplet_data[rownames(PP_summary),"CleanProb"] <- PP_summary[,"Clean"]
    x@droplet_data$Cluster <- NULL
    x@droplet_data[names(clust_max), "Cluster"] <- droplevels(clust_max)
    x@droplet_data[names(clust_prob), "ClusterProb"] <- clust_prob

    if (verbose){
        message("finished EM")
    }

    return(x)
}
