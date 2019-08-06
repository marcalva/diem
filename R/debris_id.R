
#' Initialize EM parameters for k=2 low and high groups
#'
#' Initialize the means of the two multinomial groups. Labels=2 are debris
#'
#' @return List with mu0 a p x 2 matrix and mc0 a length 2 numeric
#'
#' @importFrom Matrix rowSums
init_tails <- function(counts, labels){

    topn <- colnames(counts)[labels != 2]
    mu_top <- Matrix::rowSums(counts[,topn,drop=FALSE])
    mu_top <- mu_top/sum(mu_top)

    botn <- colnames(counts)[labels == 2]
    mu_bot <- Matrix::rowSums(counts[,botn,drop=FALSE])
    mu_bot <- mu_bot/sum(mu_bot)

    mu0 <- matrix(cbind(mu_top, mu_bot), ncol=2)
    colnames(mu0) <- c("Signal", "Background")
    rownames(mu0) <- rownames(counts)
    mc0 <- c(0.5, 0.5)
    ret <- list(Mu=mu0, Mc=mc0)
    return(ret)
}

#' x is a vector
fraction_log <- function(x){
    x_c = x - max(x)
    x_c = exp(x_c)
    frac = x_c / sum(x_c);
    return(frac);
}

#' x is a vector
sum_log <- function(x){
    max_x <- max(x)
    x_c = x - max_x
    x_sum <- log(sum(exp(x_c))) + max_x
    return(x_sum);
}

#' Get log multinomial density of columns in a sparse matrix. 
#'
#' @param X A sparseMatrix that is a sample by feature matrix of counts.
#' @param prob Probability parameter of multinomial. Must be same size as 
#'  number of columns in X.
#'
#' @import Matrix
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
#' @param labels Numeric vector of same length as number of rows in x. Fixed 
#'  the group probabilities of the integer in this vector element to 1. In 
#'  other words, the latent variable for these samples are treated as known.
#'
#' @return a numeric matrix with n samples by x variables.
dmmn <- function(x, p, mc, labels=NULL){
    if (ncol(x) != nrow(p)) stop("Number of columns in x must be the same as the number of rows in p.")
    if (length(mc) != ncol(p)) stop("Length of mc must be the same as the number of columns in p.")
    if ( any(round(colSums(p), 1) != 1) ) stop("Probability columns in p must sum to 1.")
    if (round(sum(mc), 1) != 1) stop("Mixture probabilities in mc must sum to 1.")
    if ( any(Matrix::rowSums(x) == 0) ) stop("Sample rows of x must have at least 1 count.")

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
        cat("WARNING: log probabilities returned -Inf for both groups.\n")
        cat("         Setting llk to 0 for each group.\n")
    }

    return(Llks)
}

#' Computed log likelihood of multinomial mixture for a matrix
#'
#' Given an n x p matrix of count data \code{x}, as well as a 
#' p x k matrix of probabilities \code{p}, return the probabilities 
#' of each sample for the k groups in a n x k matrix. If \code{x} 
#' has any non-integer values, they are rounded to the nearest integer.
#'
#' @param x sparseMatrix. A sample by feature matrix of counts.
#' @param p Matrix with probability parameters of multinomials for each group.
#'  Number of rows must be same as number of columns in X
#' @param mc Numeric vector with mixture coefficients. Length must be same as 
#'  number of columns in \code{p}.
#' @param Llks If log likelihood was computed for the mixture, provide the 
#'  matrix here.
#' @param labels Numeric vector of same length as number of rows in x. Fixed 
#'  the group probabilities of the integer in this vector element to 1. In 
#'  other words, the latent variable for these samples are treated as known.
#'
#' @return a numeric matrix with n samples by x variables.
dmmn_llk <- function(x, p, mc, Llks=NULL, labels=NULL){
    if (is.null(Llks)) Llks <- dmmn(x=x, p=p, mc=mc, labels=labels)
    Llks[apply(Llks, 1, function(i) all(is.infinite(i))), ] <- c(0,0)
    Llks_sum <- apply(Llks, 1, sum_log)
    llk <- sum(Llks_sum)
    return(llk)
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
#' @param log Logical. If TRUE, return matrix with log probabilities.
#'
#' @return a numeric matrix with n samples by k groups.
e_step_mn <- function(x, p, mc, Llks=NULL, labels=NULL){
    if (any(p == 0)){
        stop("Probability of a genes(s) in at least 1 group is 0, which would collapse likelihood to 0 and prevent classification. Increase cpm threshold to prevent this.")
    }
    if (is.null(Llks)) Llks <- dmmn(x=x, p=p, mc=mc, labels=labels)
    Llks[apply(Llks, 1, function(i) all(is.infinite(i))), ] <- c(0,0)
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
m_step_mn <- function(x, r){
    n <- nrow(x)
    k <- ncol(r)
    wm <- as.matrix(Matrix::t(x) %*% r)
    Mu <- sweep(wm, 2, colSums(wm), "/")
    mc <- colSums(r)
    mc <- mc/sum(mc)
    return(list(Mu=Mu, Mc=mc))
}

#' Run EM on counts to estimate multinomial mixture model
#' 
#' Run expectation maximization (EM) to estimate the parameters of the multinomial 
#' mixture model. This function takes an SCE object as input, and returns an SCE 
#' object with diem information. First genes with no expression are removed. This is 
#' to make sure the entire data likelihood doesn't collapse to 0. However, during the 
#' EM run, a gene may cease to be expressed in either of the groups, collapsing the 
#' likelihood to 0 and causing the run to fail. If this happens, it is necessary to 
#' increase \code{cpm_threshold} in \code{\link{diem}}. Then, droplets with 0 
#' counts in the resulting genes are removed. To initialize EM, the emprical means 
#' are calculated from the droplets that are fixed as debris and those that are 
#' unlabeled for the debris and clean groups, respectively. The parameters of 
#' of the multinomial are re-estimated after each run, and the algorithm converges
#' when the percent change in log likihood is less than \code{eps}. If the 
#' algorithm doesn't converge by \code{max_iter}, it breaks off.
#'
#' @param x An SCE object.
#' @param eps Numeric threshold. The EM algorithm converges when the percent change in log likihood is 
#'  less than \code{eps}.
#' @param max_iter The maximum number of iterations allowed to run.
#' @param verbose Logical indicating verbosity.
#'
#' @return An SCE object with EM output.
#' @export
run_em <- function(x, eps=1e-8, max_iter=1e3, verbose=TRUE){

    # Store EM output
    emo <- list()
    x@diem@converged  <- FALSE

    genes <- rownames(x@gene_info)[x@gene_info$exprsd]
    drops <- colnames(x@counts)
    counts <- Matrix::Matrix(x@counts[genes, drops], sparse=TRUE)
    labels <- x@labels

    # Remove genes with no expression
    no_counts <- Matrix::rowSums(counts) == 0
    if (sum(no_counts) > 0){
        cat(paste0("Removing ", as.character(sum(no_counts)), " genes with no expression.\n"))
        counts <- counts[genes,]
    }

    # Remove droplets with no genes
    no_counts <- Matrix::colSums(counts) == 0
    if (sum(no_counts) > 0){
        cat(paste0("Removing ", as.character(sum(no_counts)), " droplets with no genes detected.\n"))
        counts <- counts[,!no_counts]
        labels <- labels[!no_counts]
    }

    if (verbose){
        cat("Running EM on ", as.character(nrow(counts)), " genes and ", as.character(ncol(counts)), " droplets.\n")
        cat("Classifying ", as.character(sum(labels == 0)), " droplets.\n")
    }

    # Run EM
    mn_params <- init_tails(counts, labels)
    loglk <- -Inf
    loglks <- c(loglk)
    iter <- 1
    counts <- Matrix::t(counts)
    Llks <- dmmn(x=counts, p=mn_params$Mu, mc=mn_params$Mc, labels=labels)
    while (iter < max_iter){
        resp <- e_step_mn(x=counts, p=mn_params$Mu, mc=mn_params$Mc, Llks=Llks, labels=labels)
        mn_params <- m_step_mn(x=counts, r=resp)

        Llks <- dmmn(x=counts, p=mn_params$Mu, mc=mn_params$Mc, labels=labels)
        # Evaluate log likelihood
        if (any(apply(Llks, 1, function(i) all(is.infinite(i))))){
            cat("WARNING: log probabilities returned -Inf for both groups.\n")
        }
        loglk <- loglks[[iter]] <- sum(apply(Llks, 1, sum_log))
        emo[[iter]] <- list(Z=Llks, Mu=mn_params$Mu, Mc=mn_params$Mc, loglk=loglk)

        if (verbose){
            cat(paste0("Iteration ", as.character(iter), "; llk ", as.character(loglk), "\n"))
        }

        if (iter > 1){
            dloglk <- (loglks[[iter]] - loglks[[iter-1]])/abs(loglks[[iter-1]])
            if (dloglk < 0){
                cat("Warning: Likelihood decreased.\n")
            }
            if (dloglk < eps){
                x@diem@converged  <- TRUE
                if (verbose) cat("Converged!\n")
                break
            }
        }
        iter <- iter + 1
    }
    x@diem@PP <- t(apply(Llks, 1, fraction_log))

    x@diem@emo <- emo
    if (verbose){
        cat("Finished EM\n")
    }
    return(x)
}

#' Call clean droplets after running EM
#'
#' Call targets from droplets if the log-likelihood membership 
#' probability (the posterior probability) is higher than \code{pp_thresh}.
#'
#' @param x An SCE object.
#' @param pp_thresh Numeric threshold, where clean droplets must have a 
#'  posterior probability of at least \code{pp_thresh}.
#' @param min_genes Numeric threshold, where clean droplets must have at least 
#'  \code{min_genes} genes detected.
#'
#' @return An SCE object.
#' @export
call_targets <- function(x, pp_thresh=0.95, min_genes=200){

    PP <- x@diem@PP[rownames(x@diem@dropl_info),,drop=FALSE]

    calls <- rep("Debris", nrow(PP))
    names(calls) <- rownames(PP)

    tb <- (PP[,1] > pp_thresh) & (x@dropl_info[,"n_genes"] >= min_genes)
    target_names <- rownames(PP)[tb]
    calls[tb] <- "Clean"

    x@diem@calls <- calls

    return(x)
}

#' Return column IDs of clean droplets
#'
#' @param x An SCE object.
#'
#' @return A character vector with the called droplet IDs.
#' @export
get_clean_ids <- function(x){
    if (length(x@diem@calls) == 0) stop("Run DIEM before calling get_clean_ids")
    return(names(x@diem@calls)[x@diem@calls == "Clean"])
}

#' Add calls, PP, and LLK to dropl_info.
#'
#' @param x An SCE object.
#'
#' @return An SCE object with call info, PP, and log-likilhood added to dropl_info.
#' @export
fill_dropl_info <- function(x){
    # Fill calls
    x@dropl_info[names(x@diem@calls),"Call"] <- x@diem@calls

    # Fill PP
    x@dropl_info[rownames(x@diem@PP),"PP"] <- x@diem@PP[,1,drop=FALSE]

    # Fill droplet llk
    llk <- x@diem@emo[[length(x@diem@emo)]]$Z
    x@dropl_info[rownames(llk), "LLK"] <-  llk[,1]

    return(x)
}

#' Get percent of reads aligning to given genes.
#'
#' @param x An SCE object.
#' @param genes Genes to calculate percentage of in counts.
#' @param name Column name to place in dropl_info.
#'
#' @return An SCE object.
#' @export
#' @examples
#' mm_seur <- get_gene_pct(x=mb_sce, genes="Malat1", name="pct.malat1")
#' mt_genes <- grep(pattern="^mt-", x=rownames(mb_sce), ignore.case=TRUE)
#' mm_seur <- get_gene_pct(x=mb_sce, genes=mt_genes, name="pct.mt")
get_gene_pct <- function(x, genes, name){
    expr <- x@counts[genes,,drop=FALSE]
    if (length(expr) == 0){
        stop("None of genes found in counts.")
    }
    gene_pct <- Matrix::colSums(x@counts[genes,]) / Matrix::colSums(x@counts)
    x@dropl_info[names(gene_pct),name] <- gene_pct
    return(x)
}

