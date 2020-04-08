
#' Set debris and test droplets
#'
#' Specifies the droplets that are in the test set and debris set. 
#' The test set consists of droplets that will be classified by DIEM, 
#' the labels of which are allowed to change during the EM. 
#' The debris set consists of droplets with labels fixed to the 
#' debris group during the EM.
#' 
#' \code{min_counts} and \code{min_genes} specify the lower limit 
#' that droplets in the test set must have. The default it to 
#' set all droplets with at least 100 read counts. Optionally, the 
#' \code{top_n} parameter can specify that the top 
#' \code{top_n} are set the test set, and all others are set to the 
#' debris set.
#' If the read/UMI count of the \code{top_n} count-ranked droplet is greater 
#' than \code{min_counts}, then the \code{min_counts} is 
#' increased to this limit. The \code{debris_ids} parameter allows 
#' additional flexibility in specifying those droplets that are 
#' included in the debris set. For example, it may be useful to 
#' set droplets above a pre-determined MT% threshold as debris.
#'
#' @param x An SCE object.
#' @param min_counts A numeric value that specifies the minimum 
#'  total read/UMI counts a test droplet can have.
#' @param min_genes A numeric value that specifies the minimum 
#'  total genes detected a test droplet can have. Set to 0 by default.
#' @param top_n A numeric value giving the total number of test droplets 
#'  and take the counts of the \code{top_n} ranked droplet. Turned off
#'  by default.
#' @param debris_ids A character vector of droplet IDs that will be assigned 
#'  to the debris set, regardless of its total counts or genes detected.
#' @param verbose verbosity
#'
#' @return An SCE object.
#' @importFrom Matrix colSums
#' @export
set_debris_test_set <- function(x, 
                                min_counts = 100, 
                                min_genes = 0, 
                                top_n = NULL, 
                                debris_ids = NULL, 
                                verbose = TRUE){
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
    min_counts <- max(min_counts, dc[o[top_n]], na.rm=TRUE)
    ts <- names(dc)[dc >= min_counts]
    x@test_set <- setdiff(ts, debris_ids)
    x@bg_set <- setdiff(colnames(x@counts), x@test_set)

    if (length(x@test_set) <= 1) stop("1 or less test droplets found. Check ", 
                                      sQuote("min_counts"), " and ", sQuote("min_genes"))

    if (verbose){
        message(length(x@test_set), " droplets in the test set ", 
                "and ", length(x@bg_set), " droplets in the debris set")
    }

    x@test_data <- x@droplet_data[x@test_set,]

    return(x)
}

#' Filter out lowly expressed genes
#'
#' This function removes genes that are lowly expressed. 
#' The droplets are split into the test set and debris set, and the Counts 
#' Per Million mapped reads (CPM) are calculated for each. Then, only genes 
#' that have a CPM greater than \code{cpm_thresh} in across the entire 
#' data set are kept.
#' 
#' @param x An SCE object.
#' @param cpm_thresh CPMs must be greater than this value. This is 0 by 
#'  default to include all expressed genes.
#' @param verbose verbosity
#'
#' @return An SCE object
#' @importFrom Matrix rowSums
#' @export
filter_genes <- function(x, cpm_thresh = 0, verbose = TRUE){
    if (length(x@test_set) == 0 || length(x@bg_set) == 0)
        stop("no test droplets found. Run ", sQuote("set_debris_test_set"), " before setting cluster droplets.")
    groups <- list(c(x@test_set, x@bg_set))
    keep_all <- sapply(groups, function(g){
                       expr <- rowSums(x@counts[,g,drop=FALSE])
                       cpm <- 1e6*expr/sum(expr)
                       keep <- cpm > cpm_thresh
                       return(keep)
                           })
    keep <- apply(as.matrix(keep_all), 1, all)
    expr <- rowSums(x@counts)
    cpm <- 1e6*expr/sum(expr)
    keep <- cpm > cpm_thresh
    if (sum(keep) == 0){
        stop("no genes pass ", sQuote("cpm_thresh"), " threshold")
    }
    x@gene_data[,"exprsd"] <- keep
    if (verbose){
        message("using ", sum(keep), " expressed genes")
    }
    return(x)
}

