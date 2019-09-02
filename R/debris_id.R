
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

    if (!"CleanProb" %in% colnames(x@droplet_data)) stop("Run DIEM before calling targets")

    calls <- rep("Debris", nrow(x@droplet_data))
    calls[x@droplet_data$CleanProb >= pp_thresh & x@droplet_data$n_genes >= min_genes] <- "Clean"
    calls[x@droplet_data$CleanProb < pp_thresh & x@droplet_data$n_genes >= min_genes] <- "Removed"
    calls <- as.factor(calls)
    x@droplet_data[,"Call"] <- calls

    return(x)
}

#' Return column IDs of clean droplets
#'
#' @param x An SCE object.
#'
#' @return A character vector with the called droplet IDs.
#' @export
get_clean_ids <- function(x){
    if (!"Call" %in% colnames(x@droplet_data)) stop("Call targets before calling get_clean_ids")
    return(rownames(x@droplet_data)[x@droplet_data$Call == "Clean"])
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
#' mt_genes <- grep(pattern="^mt-", x=rownames(mb_sce@gene_data), value=TRUE, ignore.case=TRUE)
#' mm_seur <- get_gene_pct(x=mb_sce, genes=mt_genes, name="pct.mt")
get_gene_pct <- function(x, genes, name){
    expr <- x@counts[genes,,drop=FALSE]
    if (length(expr) == 0){
        stop("None of genes found in counts.")
    }
    gene_pct <- 100 * Matrix::colSums(expr) / Matrix::colSums(x@counts)
    x@droplet_data[names(gene_pct),name] <- gene_pct
    return(x)
}

