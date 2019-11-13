
#' Set debris and test droplets
#'
#' Specifies the droplets that are in the test set and debris set. 
#' The test set consists of droplets that will be classified by DIEM, 
#' the labels of which are allowed to change during the EM. 
#' The debris set consists of droplets with labels fixed to the 
#' debris group during the EM.
#' 
#' \code{min_counts} and \code{min_genes} specify the lower limit 
#' that droplets in the test set must have. In addition, 
#' the parameter \code{top_n} specifies that at most the top \code{top_n} 
#' should be included in the test set. The default is 10,000 
#' in accordance with the typical maximum number of droplets 
#' generated during a single-cell experiment, but should be changed 
#' accordingly.
#' If the read/UMI count of the \code{top_n} count-ranked droplet is greater 
#' than \code{min_counts}, then the \code{min_counts} is 
#' increased to this limit. The \code{fix_debris} parameter allows 
#' additional flexibility in specifying those droplets that are 
#' included in the debris set. For example, it may be useful to 
#' set droplets above a pre-determined MT% threshold as debris, 
#' regardless of the total number of read/UMI counts present.
#'
#' @param x An SCE object.
#' @param top_n A numeric value giving the total number of test droplets 
#'  and take the counts of the \code{top_n} ranked droplet. Set to \code{NULL} 
#'  to disable and set all droplets above \code{min_counts} and 
#'  \code{min_genes} to the test set.
#' @param min_counts A numeric value that specifies the minimum 
#'  total read/UMI counts a test droplet can have.
#' @param min_genes A numeric value that specifies the minimum 
#'  total genes detected a test droplet can have.
#' @param fix_debris A character vector of droplet IDs that will be assigned 
#'  to the debris set, regardless of its total counts or genes detected.
#' @param verbose verbosity
#'
#' @return An SCE object.
#' @importFrom Matrix colSums
#' @export
set_debris_test_set <- function(x, 
                                top_n = 1e4, 
                                min_counts = 100, 
                                min_genes = 100, 
                                fix_debris = NULL, 
                                verbose = FALSE){
    if (is.null(top_n)){
        top_n <- ncol(x@counts)
    }
    if (top_n < 0){
        stop(sQuote("top_n"), " must be greater than 0")
    }
    if (min_counts < 0 || min_genes < 0){
        stop(sQuote("min_counts"), " and ", sQuote("min_genes"), 
             " must be greater than 0")
    }

    top_n <- min(top_n, ncol(x@counts))

    barcodes <- colnames(x@counts)
    dc <- colSums(x@counts)
    dg <- colSums(x@counts > 0)
    keep <- (dc >= min_counts) & (dg >= min_genes)

    # Keep top top_n
    dc <- dc[keep]
    o <- order(dc, decreasing=TRUE)
    dco <- dc[o]
    min_counts <- max(min_counts, dco[top_n], na.rm=TRUE)
    ts <- names(dco)[dco >= min_counts]
    x@test_set <- setdiff(ts, fix_debris)
    x@bg_set <- setdiff(colnames(x@counts), x@test_set)

    if (length(x@test_set) <= 1) stop("1 or less test droplets found. Check ", 
                                      sQuote("min_counts"), " and ", sQuote("min_genes"))

    if (verbose){
        message("using ", length(x@test_set), " droplets in the test set ", 
                "and ", length(x@bg_set), " droplets in the debris set")
    }

    return(x)
}

#' Set droplets for cluster initialization
#' 
#' This function sets the droplets that will be used for the 
#' initialization. The cell types from clustering these droplets are used to
#' to intialize the parameters of the multinomal for the EM. The top 
#' \code{cluster_n} droplets ranked by either gene or count are fixed as 
#' to the debris group.
#' 
#' @param x An SCE object.
#' @param cluster_n Numeric value specifying the number of droplets to use 
#'  in the cluster set. The top \code{cluster_n} droplets, ranked by 
#'  count or gene, are included in the cluster set.
#' @param order_by Whether to order the droplets by total number of total 
#'  counts or total number of genes detected.
#' @param verbose verbosity
#'
#' @return An SCE object.
#' @importFrom Matrix colSums
#' @export
set_cluster_set <- function(x, 
                            cluster_n = 1000, 
                            order_by = "gene", 
                            verbose = FALSE){
    if (cluster_n <= 0) stop(sQuote("cluster_n"), " must be greater than 0")
    if (order_by != "gene" & order_by != "count") stop(sQuote("order_by"), " must be one of ", 
                                                       dQuote("gene"), " or ", 
                                                       dQuote("count"))
    if (length(x@test_set) == 0)
        stop("No test droplets found. Run ", sQuote("set_debris_test_set"), " before setting cluster droplets.")
    if (order_by == "gene") totals <- colSums(x@counts > 0)
    else totals <- colSums(x@counts)

    cluster_n <- min(cluster_n, length(x@test_set))
    totals <- totals[x@test_set]
    o <- order(totals, decreasing=TRUE)
    totals <- totals[o]
    min_c_counts <- totals[cluster_n]
    x@cluster_set <- names(totals)[totals >= min_c_counts]

    if (verbose){
        message("Using top ", length(x@cluster_set), 
                " droplets ranked by total ", order_by, 
                "s for clustering")
    }

    return(x)
}

#' Filter out lowly expressed genes
#'
#' This function removes genes that are lowly expressed. 
#' The droplets are split into the test set and debris set, and the Counts 
#' Per Million mapped reads (CPM) are calculated for each. Then, only genes 
#' that have a CPM of at least \code{cpm_thresh} in both the test set and 
#' debris set are kept.
#' 
#' @param x An SCE object.
#' @param cpm_thresh The minimum CPM threshold for removing genes.
#' @param verbose verbosity
#'
#' @return An SCE object
#' @importFrom Matrix rowSums
#' @export
filter_genes <- function(x, cpm_thresh = 10, verbose = FALSE){
    if (length(x@test_set) == 0 || length(x@bg_set) == 0)
        stop("No test droplets found. Run ", sQuote("set_debris_test_set"), " before setting cluster droplets.")
    groups <- list(x@test_set, x@bg_set)
    keep_all <- sapply(groups, function(g){
                       expr <- rowSums(x@counts[,g,drop=FALSE])
                       cpm <- 1e6*expr/sum(expr)
                       keep <- cpm >= cpm_thresh
                       return(keep)
                           })
    keep <- apply(keep_all, 1, all)
    if (sum(keep) == 0){
        stop("no genes pass ", sQuote("cpm_thresh"), " threshold")
    }
    x@gene_data[,"exprsd"] <- keep
    if (verbose){
        message("using ", sum(keep), " expressed genes")
    }
    return(x)
}

