
#' Initialize randomly
#'
#' @export
init_random <- function(X, 
                        N, 
                        G, 
                        K, 
                        labs = NULL, 
                        n_start = 1, 
                        Alpha0 = 1e-3, 
                        Beta0 = 1e-3){
    clusts <- 1:K
    if (length(labs) != N) labs <- rep(0, N)
    if (length(labs) != N){
        stop("Size of labs must be the same as the number of rows in x.")
    }

    Xt <- t(X)

    Alpha <- rep(Alpha0, K)
    
    # Get labeled groups
    clusts <- 1:K
    labeled <- setdiff(unique(labs), 0)
    unlabeled <- setdiff(clusts, labeled)
    
    labs_init <- labs
    names(labs_init) <- seq(N)

    lbeta <- sapply(labeled, function(k){
                    ru <- rowSums(Xt[,labs_init == k]) + Beta0
                    return(ru)
                    })
    print(head(lbeta))

    uix <- seq(N)[labs_init == 0]
    params <- lapply(1:n_start, function(repl) {
                         labs_init <- labs
                         names(labs_init) <- seq(N)
                         uix_s <- sample(uix, size = length(unlabeled))
                         labs_init[uix_s] <- unlabeled
                         ubeta <- as.matrix(Xt[,uix_s]) + Beta0
                         colnames(ubeta) <- NULL
                         thisBeta <- cbind(lbeta, ubeta)
                         thisAlpha <- Alpha
                         return(list("Alpha" = thisAlpha, 
                                     "Beta" = thisBeta, 
                                     "R" = NULL))
                         })

    return(params)
}


#' fraction of logs
#'
#' @param x numeric vector
#' @export
fraction_log <- function(x){
    x_c = x - max(x)
    x_c = exp(x_c)
    frac = x_c / sum(x_c);
    return(frac);
}

#' sum of logs
#'
#' @param x numeric vector
#' @export
sum_log <- function(x){
    max_x <- max(x)
    x_c = x - max_x
    x_sum <- log(sum(exp(x_c))) + max_x
    return(x_sum);
}

#' Normalizing constant for Dirichlet
#' @export
dir_c <- function(x, log=TRUE){
    if (log){
        lgamma(sum(x)) - sum(lgamma(x))
    } else {
        gamma(sum(x)) / prod(gamma(x))
    }
}

#' Entropy of Dirichlet
#'
hd <- function(a){
    K <- length(a)
    a0 <- sum(a)
    (digamma(a0) * (a0 - K) )- (digamma(a) %*% (a-1)) - dir_c(a)
}

#' Get expected log rho
#' @export
get_e_rho <- function(Beta){
    apply(Beta, 2, function(b){ digamma(b) - digamma(sum(b)) })
}

#' Get expected log pi
#' @export
get_e_pi <- function(Alpha){
    digamma(Alpha) - digamma(sum(Alpha))
}

#' Get expected Z
#' @export
get_e_z <- function(R){
    R
}

# Store m_num
#' @export
get_mn_num <- function(X, K){
    N <- nrow(X)
    nsum <- rowSums(X)
    mnnum <- lfactorial(nsum)
    mnnum <- matrix(mnnum, nrow = N, ncol = K, byrow=FALSE)
    return(mnnum)
}

# Store m_den
#' @export
get_mn_den <- function(X, K){
    N <- nrow(X)
    Xlf <- X
    Xlf@x <- lfactorial(Xlf@x)
    xlogg <- rowSums(Xlf)
    xlogg <- matrix(xlogg, nrow = N, ncol = K, byrow=FALSE)
    return(xlogg)
}

# Store Xtrho
#' @export
get_Xrho <- function(X, e_rho){
    return(as.matrix(X %*% e_rho))
}

#' Update Z
#' @export
update_z <- function(X, 
                     K, 
                     e_rho, 
                     e_pi, 
                     mn_num, 
                     mn_den, 
                     Xrho, 
                     labs = NULL){
    groups <- 1:K
    N <- nrow(X)
    G <- ncol(X)
    if (length(labs) != N) labs <- rep(0, N)

    if (is.null(rownames(X))) rownames(X) <- seq(N)
    seq_n <- seq(N)
    fixd <- seq_n[labs != 0]
    unfixd <- seq_n[labs == 0]

    Xu <- X[unfixd,,drop=FALSE]
    mn_num_u <- mn_num[unfixd,,drop=FALSE]
    mn_den_u <- mn_den[unfixd,,drop=FALSE]
    # Xurho <- get_Xrho(Xu, e_rho)
    Xurho <- Xrho[unfixd,,drop=FALSE]
    e_pi_u <- matrix(e_pi, nrow = nrow(Xu), ncol = K, byrow = TRUE)

    zu <- as.matrix(mn_num_u - mn_den_u + Xurho + e_pi_u)
    rownames(zu) <- rownames(X)[unfixd]

    # Get E(Z) for fixed observations
    zf <- matrix(-Inf, nrow = length(fixd), ncol = K)
    rownames(zf) <- rownames(X)[fixd]
    rv <- 1:length(fixd)
    cv <- labs[fixd]
    zf[rv + nrow(zf) * (cv - 1)] <- 0

    z <- rbind(zu, zf)
    z <- z[rownames(X),]

    z <- apply(z, 1, fraction_log)
    z <- t(z)

    return(z)
}

#' @export
update_alpha <- function(hyperp, e_z){
    colSums(e_z) + hyperp$Alpha
}

#' @export
update_beta <- function(X, K, hyperp, e_z){
    Beta0_m <- matrix(hyperp$Beta, nrow = length(hyperp$Beta), ncol = K, byrow=FALSE)
    return(Beta0_m + as.matrix(t(X) %*% e_z))
}

#' Compute lower bound
#'
#' @export
lower_bound <- function(X, 
                        K, 
                        params, 
                        hyperp, 
                        e_z, 
                        e_rho, 
                        e_pi, 
                        mn_num, 
                        mn_den, 
                        Xrho){
    groups <- 1:K
    N <- nrow(X)
    G <- ncol(X)

    p_x <- sum(e_z * (mn_num - mn_den + Xrho))

    
    e_pi_m <- matrix(e_pi, nrow = nrow(X), ncol = K, byrow = TRUE)
    p_z <- sum(e_z * e_pi_m)

    p_pi <- dir_c(hyperp$Alpha) + t(hyperp$Alpha - 1) %*% e_pi

    p_rho <- sum(rep(dir_c(hyperp$Beta), K) + (t(hyperp$Beta - 1) %*% e_rho))

    ln_e_z <- log(e_z)
    ln_e_z[is.infinite(ln_e_z)] <- 0
    q_z <- e_z * ln_e_z
    q_z <- sum(q_z)

    q_pi <- dir_c(params$Alpha) + ((params$Alpha - 1) %*% e_pi)

    q_rho <- apply(params$Beta, 2, dir_c) + 
        sapply(1:K, function(k){ (params$Beta[,k]- 1 ) %*% e_rho[,k] })
    q_rho <- sum(q_rho)

    ret <- p_x + p_z + p_pi + p_rho - q_z - q_pi - q_rho

    return(ret)
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
#' @param labs Numeric vector of same length as number of 
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
vem <- function(counts, 
                K, 
                params, 
                hyperp, 
                max_iter = 1e2, 
                eps = 1e-6, 
                labs = NULL, 
                verbose = TRUE){
    N <- nrow(counts)
    G <- ncol(counts)

    if (length(labs) != N) labs <- rep(0, N)
    if (nrow(params$Beta) != G) stop("Number of columns in counts must be the same as the number of features in Beta.")
    if (ncol(params$Beta) != K) stop("Number of columns in Beta must be the number of clusters K.")
    if (length(params$Alpha) != K) stop("Length of Alpha must be the number of clusters K.")
    # if (nrow(params$R) != N) stop("Number of rows in R must be the same as the number of observations in counts.")
    # if (ncol(params$R) != K) stop("Number of columns in R must be the number of clusters K.")
    if (length(labs) != N) stop("Length of labs must be the same as the number of observations in counts.")

    mn_num <- get_mn_num(counts, K)
    mn_den <- get_mn_den(counts, K)
    
    # Start with parameters Alpha and Beta. Iteration updates Z, then alpha and beta
    e_rho <- get_e_rho(params$Beta)
    Xrho <- get_Xrho(counts, e_rho)
    e_pi <- get_e_pi(params$Alpha)
    lb_old <- -Inf
    iter <- 1
    converged <- FALSE
    while (iter <= max_iter){
        if (verbose) message("Iteration ", iter)
        # Update Z
        R <- update_z(X = counts, K = K, e_rho = e_rho, 
                      e_pi = e_pi, mn_num = mn_num, mn_den = mn_den,
                      Xrho = Xrho, labs = labs)
        e_z <- get_e_z(R)

        # Update Alpha and Beta
        Alpha <- update_alpha(hyperp, e_z)
        e_pi <- get_e_pi(Alpha)
        Beta <- update_beta(X = counts, K, hyperp, e_z)
        e_rho <- get_e_rho(Beta)
        Xrho <- get_Xrho(counts, e_rho)

        params <- list("R" = R, "Alpha" = Alpha, "Beta" = Beta)

        lb_new <- lower_bound(X = counts, K = K, params = params, hyperp = hyperp, 
                              e_z = e_z, e_rho = e_rho, e_pi = e_pi, 
                              mn_num = mn_num, mn_den = mn_den, Xrho = Xrho)

        if (iter > 1) delt <- (lb_new - lb_old) / abs(lb_old)
        else delt <- Inf

        if (delt < 0) warning("lower bound decreased")
        if (delt < eps){
            if (verbose) message("converged")
            converged <- TRUE
            break
        }
        lb_old <- lb_new
        iter <- iter + 1
    }
    vmo <- list("params" = params, 
                "lb" = lb_new, 
                "converged" = converged, 
                "n_iter" = iter)
    return(vmo)
}

#' Initialize best from random starts
init_best <- function(X, 
                         K, 
                         hyperp, 
                         n_start = 100, 
                         n_iter = 3, 
                         labs = NULL, 
                         verbose = TRUE){
    N <- nrow(X)
    G <- ncol(X)
    pl <- init_random(X = X, N = N, G = G, K = K, labs = labs, 
                      Alpha = hyperp$Alpha[1], Beta = hyperp$Beta[1], 
                      n_start = n_start)
    print(class(pl[[1]]$Beta))
    vl <- lapply(pl, function(p) {
                 vem(counts = X, K = K, params = p, hyperp = hyperp, 
                     max_iter = n_iter, eps = -Inf, labs = labs)
                 })
    return(vl)
    return(vmol)
    #lbs <- sapply(vmol, function(i) i$lb)
    #i <- which.max(lbs)
    #return(vmol[[i]])
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
#' @param verbose Logical indicating verbosity.
#'
#' @return An SCE object.
#' @importFrom igraph gsize
#' @importMethodsFrom Matrix %*%
#' @export
run_vem <- function(x, 
                    K = 30, 
                    Alpha0 = 1, 
                    Beta0 = 5, 
                    eps = 1e-6, 
                    max_iter = 1e2, 
                    verbose = TRUE){

    genes.use <- rownames(x@gene_data)[x@gene_data$exprsd]
    # droplets.use <- c(x@test_set, x@bg_set)
    droplets.use <- colnames(x@counts)

    if (length(genes.use) == 0 | length(droplets.use) == 0) stop("Specify test set and filter genes before running EM.")

    counts <- t(x@counts[genes.use,droplets.use])

    N <- nrow(counts)
    G <- ncol(counts)

    # Fix labs
    labs <- rep(0, N)
    names(labs) <- droplets.use
    labs[x@bg_set] <- 1

    hyperp <- list("Alpha" = rep(Alpha0, K), 
                   "Beta" = rep(Beta0, G))

    # Get initial parameters
    if (is.null(x@ic)){
        countsv <- counts[,x@vg]
        Gs <- ncol(countsv)
        hyperps <- list("Alpha" = rep(Alpha0, K), "Beta" = rep(Beta0, Gs))
        vmol <- init_best(countsv, K = K, hyperp = hyperps, labs = labs)
        lbs <- sapply(vmol, function(i) i$lb)
        i_max <- which.max(lbs)
        params <- vmol[[i_max]]
    } else {
        params <- list()
        clusts <- 1:K
        Beta <- sapply(clusts, 
                       function(i){
                           rowSums(x@counts[genes.use, x@ic$labels == i, drop=FALSE])})
        params$Beta <- Beta + Beta0
        params$R <- NULL
        params$Alpha <- rep(Alpha0, K)
    }

    if (verbose){
        message("Estimating parameters")
    }

    # Run EM
    vmo <- vem(counts, 
               K = K, 
               params = params,
               hyperp = hyperp,
               max_iter = max_iter, 
               eps = eps, 
               labs = labs, 
               verbose = verbose)

    # Naive Bayes assignment

    clusts <- colSums(vmo$params$R) > 0
    vmo$params$R <- vmo$params$R[,clusts,drop=FALSE]
    vmo$params$Beta <- vmo$params$Beta[,clusts,drop=FALSE]
    vmo$params$Alpha <- vmo$params$Alpha[clusts]
    x@assignments <- x@assignments[clusts]

    a <- x@assignments
    PP_summary <- sapply(levels(a), function(i){
                         rowSums(vmo$params$R[,a == i,drop=FALSE])})
    clust_prob <- apply(vmo$params$R, 1, function(i) i[which.max(i)])
    clust_max <- apply(vmo$params$R, 1, function(i){
                       (1:ncol(vmo$params$R))[which.max(i)]
               })
    #clust_max <- paste0("C", as.character(clust_max))
    #clust_max <- as.factor(clust_max)
    #clust_max <- droplevels(clust_max)
    #cl <- paste0("C", seq(1,nlevels(clust_max)-1))
    #levels(clust_max) <- c("Debris", cl)

    # Fill in droplet data
    x@droplet_data[,"CleanProb"] <- PP_summary[,"Clean"]
    x@droplet_data$Cluster <- NULL
    x@droplet_data$Cluster <- clust_max
    x@droplet_data$ClusterProb <- clust_prob

    x@emo <- vmo

    if (verbose){
        message("finished variational EM")
    }

    return(x)
}

