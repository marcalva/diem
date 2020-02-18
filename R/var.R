
#' Get variable genes
#'
#' This function finds the variable genes in an SCE object. 
#' The means and variances of the genes are calculated from the raw counts 
#' of the cluster set. Then, the means and variances are log-transformed 
#' after adding a constant of 1. A loess regression line is fit between 
#' the log counts and log variance, and the only top genes 
#' ranked by residual are variable genes used to initialize the clusters. 
#' The number of genes is specified with \code{n_genes}. 
#' The span for the loess function is given by \code{lss} 
#' (default is 0.3).
#'
#' @param x An SCE object.
#' @param n_genes Number of variable genes to return.
#' @param lss Numeric value of the span parameter of the loess regression.
#' @param droplets.use Vector of droplet IDs to use for getting variable 
#'  genes. Default is to use the test set.
#' @param verbose Verbosity.
#'
#' @return An SCE object
#' @importFrom Matrix rowMeans
#' @importFrom stats loess
#' @export
get_var_genes <- function(x, 
                          n_genes=2000, 
                          lss=0.3, 
                          droplets.use = NULL, 
                          verbose=FALSE){
    if (verbose) message("getting variable genes")

    if (is.null(droplets.use)){
        droplets.use <- x@test_set
    }
    counts <- x@counts[,droplets.use]
    if (sum(x@gene_data$exprsd) == 0) x <- filter_genes(x)
    counts <- counts[x@gene_data$exprsd,]
    exprsd <- rowSums(counts) > 0
    counts <- counts[exprsd,,drop=FALSE]
    gene_means <- rowMeans(counts)
    gene_names <- rownames(counts)

    log_mean <- log10(gene_means + 1)
    log_var <- log10(fast_varCPP(counts, gene_means) + 1)
    fit <- loess(log_var ~ log_mean, span=lss)
    rsd <- log_var - fit$fitted
    topi <- order(rsd, decreasing=TRUE)[1:n_genes]
    vg <- gene_names[topi]
    datf <- data.frame(mean=log_mean, var=log_var,
                       fit=fit$fitted, rsd=rsd)
    rownames(datf) <- gene_names
    x@vg <- vg[!is.na(vg)]
    x@vg_info <- datf
    if (verbose) message("found ", n_genes, " variable genes")
    return(x)
}

