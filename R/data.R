
#' Divide elements of a column by the column's sum in a sparse matrix
#'
#' @param x Sparse Matrix
#'
#' @return The Sparse Matrix x with columns summing to 1.
#' @importFrom Matrix Diagonal colSums
#' @export
divide_by_colsum <- function(x){
    cs <- Matrix::colSums(x = x)
    d <- Matrix::Diagonal(x = 1/cs)
    x <- x %*% d
    return(x)
}

#' Normalize counts of a sparse matrix
#'
#' @param counts Sparse Matrix
#' @param scale_factor Numeric. A scaling factor to multiply values after division by total count.
#'  Default is 1.
#' @param logt Boolean. Log transform after normalizing for columns to sum to 1 and multiplying by scaling factor.
#'  Default is TRUE
#'
#' @return Sparse Matrix
#' @export
norm_counts <- function(counts, scale_factor=1, logt=TRUE){
    counts <- divide_by_colsum(counts)
    counts <- counts * scale_factor
    if (logt) counts <- log1p(counts)
    return(counts)
}

#' Normalize raw counts.
#'
#' Counts are normalized by dividing by total counts per droplet. Note that this only normalizes 
#' droplets that are candidates in the test set. Optionally, droplets can be normalized 
#' by log1p transformation, with a counts being multiplied by a scaling factor \code{scale_factor}.
#'
#' @param x SCE.
#' @param genes.use A character vector to subset count data to. Normalization will only be run for these genes.
#' @param scale_factor A numeric scaling factor to multiply counts after division by column sum
#'  Default is 1.
#' @param logt Logical indicating whether to log1p transform expression values.
#'  Default is TRUE
#'
#' @return SCE object
#' @importFrom Matrix rowMeans
#' @export
normalize_data <- function(x, 
                           genes.use=NULL, 
                           droplets.use=NULL, 
                           scale_factor=1, 
                           logt=TRUE){
    if (is.null(genes.use)) genes.use <- rownames(x@gene_data)[x@gene_data[,"exprsd"]]
    if (is.null(droplets.use)) droplets.use <- x@cluster_set

    if (length(genes.use) == 0) stop("0 genes specified in genes.use.")
    if (length(droplets.use) == 0) stop("0 droplets specified in droplets.use")

    expr <- x@counts[genes.use, droplets.use]
    x@norm <- norm_counts(expr, scale_factor, logt=logt)

    return(x)
}

#' Get PCs from normalized count data
#'
#' @param x An SCE object.
#' @param n_pcs Number of PCs to output.
#' 
#' @return An SCE object.
#' @importFrom irlba prcomp_irlba
#' @export
get_pcs <- function(x, n_pcs=30){
    if (length(x@norm) == 0) stop("Normalize counts before getting PCs.")
    nc <- Matrix::t(x@norm)
    keep <- Matrix::colSums(nc) > 0
    nc <- nc[,keep]
    x@pcs <- irlba::prcomp_irlba(nc, n=n_pcs, scale.=TRUE)$x
    rownames(x@pcs) <- rownames(nc)
    return(x)
}

#' Set droplets for EM testing
#'
#' Sets the count threshold \emph{T} for testing droplets. Droplets with counts below \emph{T} 
#' are fixed as debris, while droplets greater than or equal to \emph{T} are set as unlabeled. The 
#' threshold \emph{T} is found by ranking droplets in order of decreasing size (total read/UMI counts) 
#' and taking the read/UMI count size of the \code{top_n} barcode. If this ranked count is less than 
#' \code{min_counts}, then set \emph{T} to \code{min_counts} instead. In addition to fixing the labels 
#' of debris droplets, it is possible to fix droplets as clean. Droplets with counts above \code{top_thresh} 
#' are fixed to nuclei. By default, no droplets are fixed as clean and this is not recommended.
#'
#' @param x An SCE object.
#' @param top_n Numeric value specifying the  droplets below \code{top_n} (ranked by total counts) 
#'  are fixed to debris. Set to NULL to ignore.
#' @param min_counts A numeric value that specifies that all droplets below this count threshold are 
#'  always fixed to debris.
#' @param top_thresh A numeric value specifying any droplets above this count threshold are 
#' fixed to nuclei. No droplets are fixed to nuclei by default.
#'
#' @return An SCE object.
#' @importFrom Matrix colSums
#' @export
set_test_set <- function(x, 
                         top_n=1e4, 
                         min_counts=150, 
                         min_genes=150, 
                         cluster_n=NULL, 
                         cluster_quantile=0.995, 
                         cluster_divide=3){
    if (is.null(top_n)){
        top_n <- ncol(x@counts)
    }
    if (top_n < 0){
        stop("top_n must be greater than 0")
    }
    if (min_counts < 0){
        stop("min_counts must be greater than 0")
    }

    top_n <- min(top_n, ncol(x@counts))

    barcodes <- colnames(x@counts)
    dc <- Matrix::colSums(x@counts)
    dg <- Matrix::colSums(x@counts > 0)
    keep <- (dc >= min_counts) & (dg >= min_genes)

    # Keep top top_n
    dc <- dc[keep]
    o <- order(dc, decreasing=TRUE)
    dco <- dc[o]
    min_counts <- max(min_counts, dco[top_n], na.rm=TRUE)
    ts <- names(dco)[dco >= min_counts]
    x@test_set <- ts
    x@bg_set <- setdiff(colnames(x@counts), ts)
    x@min_counts <- min_counts
    
    if (is.null(cluster_n)){
        # Get cluster set
        dg <- dg[x@test_set]
        o <- order(dg)
        dgo <- dg[o]
        quant <- dgo[floor(length(dgo)*cluster_quantile)]
        min_count <- quant/cluster_divide
        x@cluster_set <- names(dg[dg >= min_count])
    }

    return(x)
}

#' Filter out lowly expressed features
#'
#' This function removes genes that are lowly expressed in 
#' either the nuclear or background group. Read/UMI counts are summed across 
#' droplets in the nuclear-enriched and background-enriched groups. The 
#' background-enriched droplets are those with labels fixed to debris, 
#' (see \code{\link{set_test_set}}). Counts per million mapped reads (CPM) 
#' are calculated for the two groups, and only genes with a CPM of at least 
#' \code{cpm_thresh} in both groups are kept. This function ensures that 
#' the likelihood of the multinomial mixture does not collapse to 0, since 
#' a gene with 0 counts has a likelihood of 0. This creates a logical 
#' vector in the gene_info slot that indicates whether a gene passes the 
#' threshold or not.  #'
#' @param x An SCE object.
#' @param cpm_thresh The minimum CPM expression threshold in both groups.
#'
#' @return An SCE object
#' @importFrom Matrix rowSums
#' @export
filter_genes <- function(x, cpm_thresh=10){
    groups <- list(x@test_set, x@bg_set)
    keep_all <- sapply(groups, function(g){
                       expr <- Matrix::rowSums(x@counts[,g])
                       cpm <- 1e6*expr/sum(expr)
                       keep <- cpm >= cpm_thresh
                       return(keep)
                           })
    keep <- apply(keep_all, 1, all)
    if (sum(keep) == 0){
        stop("No genes pass cpm_thresh threshold.")
    }
    x@gene_data[,"exprsd"] <- keep
    return(x)
}

#' Read 10X output
#'
#' Read output from 10X into a sparse matrix. Given a 10X output path, 
#' containing the output files from 10X CellRanger. These should be the 
#' raw output files \code{barcodes.tsv} \code{genes.tsv} \code{matrix.mtx} 
#' if using CellRanger v2, or \code{barcodes.tsv.gz} \code{features.tsv.gz} 
#' \code{matrix.mtx.gz} if using CellRanger v3.
#'
#' @param path Character. File path prefix to 10X output.
#'
#' @return expr Sparse matrix. Expression counts from 10X cell ranger
#'  with cells in the columns and genes in the rows.
#'
#' @import Matrix
#' @export
#'
#' @examples
#' counts <- read_10x("mouse_nuclei_2k/raw_gene_bc_matrices/mm10/")
#'
read_10x <- function(path){
    files <- list.files(path, full.names=TRUE)
    files_names <- list.files(path, full.names=FALSE)
    v3 <- "matrix.mtx.gz" %in% files_names

    if (v3){
        mtx_file <- paste0(path, "/matrix.mtx.gz")
        genes_file <- paste0(path, "/features.tsv.gz")
        barcode_file <- paste0(path, "/barcodes.tsv.gz")
    } else {
        mtx_file <- paste0(path, "/matrix.mtx")
        genes_file <- paste0(path, "/genes.tsv")
        barcode_file <- paste0(path, "/barcodes.tsv")
    }

    if (file.exists(mtx_file)) {
        expr <- readMM(mtx_file)
    } else {
        stop(paste0(mtx_file, " not found"))
    }
    expr <- as(expr, "CsparseMatrix")

    if (file.exists(barcode_file)){
        barcode_names <- readLines(barcode_file)
        barcode_names <- sub("-.*", "", barcode_names)
    } else {
        stop(paste0(barcode_file, " not found"))
    }

    if (file.exists(genes_file)){
        genes <- read.delim(genes_file, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
    } else {
        stop(paste0(genes_file, " not found"))
    }

    colnames(expr) <- barcode_names
    rownames(expr) <- make.unique(genes[,2])

    return(expr)
}

