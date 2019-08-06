#' Divide elements of a column by the column's sum in a sparse matrix
#'
#' @param x Sparse Matrix
#'
#' @return The Sparse Matrix x with columns summing to 1.
#' @importFrom Matrix Diagonal colSums
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
#' Additionally, gene counts can be scaled to mean 0 and variance 1 (default). This is recommended so that 
#' highly expressed genes do not drive the calculation of pi. Note that, when scaling, genes with 0 mean 
#' and zero variance will be removed automatically.
#'
#' @param x SCE.
#' @param genes A character vector to subset count data to. Normalization will only be run for these genes.
#' @param scale_factor A numeric scaling factor to multiply counts after division by column sum
#'  Default is 1.
#' @param logt Logical indicating whether to log1p transform expression values.
#'  Default is TRUE
#'
#' @return SCE object
#' @importFrom Matrix rowMeans
#' @export
normalize <- function(x, 
					  genes=NULL, 
					  scale_factor=1,
					  logt=TRUE){
	if (sum(x@labels) == 0){
		stop("Specify test set with set_test_set function")
	}

	expr <- x@counts[,sum(x@labels) == 0]
	if (!is.null(genes)) expr <- expr[genes,]

	x@norm <- norm_counts(expr, scale_factor, logt=logt)

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
#'  are fixed to debris.
#' @param min_counts A numeric value that specifies that all droplets below this count threshold are 
#'  always fixed to debris.
#' @param top_thresh A numeric value specifying any droplets above this count threshold are 
#' fixed to nuclei. No droplets are fixed to nuclei by default.
#'
#' @return An SCE object.
#' @importFrom Matrix colSums
#' @export
set_test_set <- function(x, top_n=1e4, min_counts=100, top_thresh=NULL){
	if (is.null(top_n)){
		top_n <- ncol(x@counts)
	}
	if (top_n < 0){
		stop("top_n must be greater than 0")
	}
	if (min_counts < 0){
		stop("min_counts must be greater than 0")
	}

	x@labels <- rep(2, ncol(x@counts))
	names(x@labels) <- colnames(x@counts)

	dc <- Matrix::colSums(x@counts)
	top_n_count <- dc[order(dc, decreasing=TRUE)[top_n]]
	min_counts <- max(min_counts, top_n_count)
	test_set <- names(dc)[dc >= min_counts]
	x@labels[test_set] <- 0

	if (!is.null(top_thresh)){
		signal_set <- names(dc)[dc >= top_thresh]
		x@labels[signal_set] <- 1
	} else {
		top_thresh <- Inf
	}

	x@min_counts <- min_counts
	x@top_thresh <- top_thresh

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
#' threshold or not.
#'
#' @param x An SCE object.
#' @param cpm_thresh The minimum CPM expression threshold in both groups.
#'
#' @return An SCE object
#' @importFrom Matrix rowSums
#' @export
filter_genes <- function(x, cpm_thresh=25){
	if (length(x@labels) == 0){
		stop("Fix labels with set_test_set before calling filter_genes.")
	}
	bg_drops <- names(x@labels)[x@labels == 2]
	nbg_drops <- names(x@labels)[x@labels != 2]

	if (length(bg_drops) == 0 | length(nbg_drops) == 0){
		stop("No droplets in the background-enriched or nuclear-enriched groups.")
	}

	bg_expr <- Matrix::rowSums(x@counts[,bg_drops])
	nbg_expr <- Matrix::rowSums(x@counts[,nbg_drops])

	bg_cpm <- 1e6*bg_expr/sum(bg_expr)
	nbg_cpm <- 1e6*nbg_expr/sum(nbg_expr)

	keep <- (bg_cpm >= cpm_thresh) & (nbg_cpm >= cpm_thresh)
	if (sum(keep) == 0){
		stop("No genes pass cpm_thresh threshold.")
	}
	x@gene_info[,"exprsd"] <- keep
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
