
#' @export
get_var <- function(counts, cpm_thresh=10, ng=2000, lss=0.3, verbose=FALSE){
    cpm <- Matrix::rowSums(counts)
    cpm <- 1e6*cpm/sum(cpm)
    exprd <- cpm > cpm_thresh
    counts <- counts[exprd,]
    gene_means <- Matrix::rowMeans(counts)
    gene_names <- rownames(counts)

    log_mean <- log10(gene_means + 1)
    log_var <- log10(fast_varCPP(counts, gene_means) + 1)
    fit <- loess(log_var ~ log_mean, span=lss)
    rsd <- log_var - fit$fitted
    topi <- order(rsd, decreasing=TRUE)[1:ng]
    vg <- gene_names[topi]
    return(vg)
}


#' Get variable genes
#'
#' Finds top \code{nf} variable genes by running a loess regression 
#' on log(variance) and log(mean). Residuals are used to rank 
#' features, and the top \code{nf} are output.
#'
#' @param x SCE. SCE object.
#' @param min_mean Numeric. Minimum mean of counts to consider.
#' @param nf Numeric. Output top nf variable features.
#' @param lss Numeric. The span paramater for loess function.
#' @param verbose Boolean.
#'
#' @return SCE object with \code{vg} and \code{vg_info} slots filled for variable genes
#' and variable gene info. \code{vg_info} is a data frame with log mean, log var, 
#' fitted log var values, and residuals. \code{vg} is a vector of genes
#' @export
get_var_genes <- function(x, cpm_thresh=10, ng=2000, lss=0.3, verbose=FALSE){
    if (verbose) cat("Getting variable genes\n")
    
    counts <- x@counts[,x@cluster_set]
    exprd <- x@gene_data$exprsd & (Matrix::rowSums(counts) > 0)
    counts <- counts[exprd,]
    gene_means <- Matrix::rowMeans(counts)
    gene_names <- rownames(counts)

    log_mean <- log10(gene_means + 1)
    log_var <- log10(fast_varCPP(counts, gene_means) + 1)
    fit <- loess(log_var ~ log_mean, span=lss)
    rsd <- log_var - fit$fitted
    topi <- order(rsd, decreasing=TRUE)[1:ng]
    vg <- gene_names[topi]
    datf <- data.frame(mean=log_mean, var=log_var,
                       fit=fit$fitted, rsd=rsd)
    rownames(datf) <- gene_names
    x@vg <- vg
    x@vg_info <- datf
    if (verbose) cat("Found variable genes\n")
    return(x)
}


