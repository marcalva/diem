#' Divide elements of a column in a sparse matrix by the column's sum
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

#' Subset DE S4 object to genes in keep
#'
#' @param x SCE object.
#' @param keep Character vector. Names of genes to keep in DE object
#'
#' @return An SCE object
#' @export
subset_DE <- function(x, keep){
	x@de@low_means <- x@de@low_means[intersect(names(x@de@low_means), keep)]
	x@de@high_means <- x@de@high_means[intersect(names(x@de@high_means), keep)]
	x@de@deg_low <- intersect(x@de@deg_low, keep)
	x@de@deg_high <- intersect(x@de@deg_high, keep)
	x@de@log2fc <- x@de@log2fc[intersect(names(x@de@log2fc), keep)]
	return(x)
}

#' Normalize raw counts for droplets being classified.
#'
#' Counts are normalized by dividing by total counts per droplet. Optionally, droplets can be normalized 
#' by log1p transformation, with a counts being multiplied by a scaling factor \code{scale_factor}.
#' Additionally, gene counts can be scaled to mean 0 and variance 1 (default). This is recommended so that 
#' highly expressed genes do not drive the calculation of pi. Note that, when scaling, genes with 0 mean 
#' and zero variance will be removed automatically.
#'
#' @param x SCE.
#' @param scale_factor Numeric. A scaling factor to multiply values after division by total count.
#' @param logt Boolean. Whether to log1p transform expression values.
#' @return SCE object
#' @importFrom Matrix rowMeans
#' @export
normalize <- function(x, 
					  genes=NULL, 
					  scale_factor=1,
					  logt=TRUE,
					  verbose=FALSE){
	if (length(x@test_IDs) == 0){
		stop("Specify test set with set_test_set function")
	}
	if (verbose) cat("Normalizing\n")

	expr <- x@counts[,x@test_IDs]
	if (!is.null(genes)) expr <- expr[genes,]
	expr <- divide_by_colsum(expr)
	expr <- expr * scale_factor
	if (logt) expr <- log1p(expr)
	x@norm <- expr
	if (verbose) cat("Normalized\n")

	return(x)
}

#' Get pi projection of normalized data
#'
#' This function calculates pi_low and pi_high, defined as the inner products of 
#' normalized gene expression and log fold changes for DE genes enriched in the low (background-enriched) 
#' and high (target-enriched) count droplets, respectively.
#'
#' @param x SCE. SCE object.
#'
#' @return An matrix of pi values
#' @importFrom Matrix Matrix t
#' @export
get_pi <- function(x){

	x <- normalize(x) # Set @norm

	if (length(x@diem@de@prop_diff) == 0){
		x <- set_diff_genes(x)
	}

	de_obj <- x@diem@de

	log2fc <- log2(1e6*de_obj@low_prop) - log2(1e6*de_obj@high_prop)
	log2fc_l <- log2fc[log2fc > 0]
	log2fc_h <- log2fc[log2fc < 0]

	if (any(is.na(log2fc_l)) | any(is.na(log2fc_h))){
		stop("Fold changes cannot be NA.")
	}
	if (any(is.infinite(log2fc_l)) | any(is.infinite(log2fc_h))){
		stop("Fold changes cannot be infinite.")
	}

	genes_low <- intersect(rownames(x@norm), names(log2fc_l))
	genes_high <- intersect(rownames(x@norm), names(log2fc_h))

	if ( (length(genes_low) == 0) | (length(genes_high) == 0) ){
		stop("No differential genes between background and signal.")
	}

	# log2fc_l is positive. log2fc_h is negative (so take negative)
	pi_l <- t(x@norm[genes_low,]) %*% log2fc_l[genes_low]
	pi_h <- t(x@norm[genes_high,]) %*% -log2fc_h[genes_high]

	pi_df <- as.data.frame(cbind(pi_l, pi_h))
	rownames(pi_df) <- colnames(x@norm)
	colnames(pi_df) <- c("pi_l", "pi_h")

	x@diem@pi <- as.matrix(pi_df)

	return(x)
}

#' Get PCs
#'
#' @return SCE object.
#' @export
get_pcs <- function(x, n_pcs=30){
	x <- normalize(x) # Set @norm
	genes <- x@diem@de@diff_genes
	drops <- x@test_IDs
	counts <- Matrix::Matrix(x@norm[genes, drops], sparse=TRUE)
	counts <- Matrix::t(counts)
	sce@pcs <- prcomp(as.matrix(counts), scale.=TRUE)
	return(sce)
}

#' Subset droplets for EM testing
#'
#' Sets the threshold for testing droplets and subsets counts matrix. Includes the top \code{top_n} 
#' droplets for testing (including any ties), and removing any droplets with total counts 
#' less than \code{min_counts}. The count matrix is subset to these droplets for EM testing 
#'
#' @param x SCE. SCE object.
#' @param top_n Numeric. Only the top \code{top_n} droplets ranked by total counts are classified.
#' @param min_counts Integer. Only droplets with at least this number of total counts are classified. 
#'  Can override \code{top_n} if \code{top_n} includes droplets with less than \code{min_counts}.
#'
#' @return An SCE object.
#' @importFrom Matrix colSums
#' @export
set_test_set <- function(x, top_n=1e4, min_counts=30, verbose=FALSE){
	if (is.null(top_n)){
		top_n <- ncol(x@counts)
	}
	if (is.null(min_counts)){
		min_counts <- 1
	}

	if (top_n < 0){
		stop("top_n must be greater than 0")
	}
	if (min_counts < 0){
		stop("min_counts must be greater than 0")
	}

	dc <- Matrix::colSums(x@counts)
	top_n_ix <- order(dc, decreasing=TRUE)[top_n]
	min_counts <- max(dc[top_n_ix], min_counts)

	dc <- dc[dc >= min_counts]
	if (length(dc) < 50){
		stop("Less than 50 barcodes pass filtering. Choose less stringent top_n and min_counts parameters.")
	}

	x@test_IDs <- names(dc)
	x@test_thresh <- min_counts

	if (verbose){
		cat(paste0("Testing ", as.character(length(dc)), " barcodes that have counts greater than or equal to ",
				   as.character(min_counts), ".\n"))
	}
	return(x)
}

#' Fix labels of lower and upper count droplets 
#'
#' Fix labels of the upper and lower tails droplets to signal and background, respectively.
#' The top and bottom \code{pct} are then fixed during EM.
#'
#' @param x SCE. SCE object.
#' @param pct Numeric. Fraction of total test barcodes in each tail to fix
#'  as the background or signal.
#'
#' @return An SCE object
#' @importFrom Matrix colSums
#' @export
fix_tails <- function(x, n_top=1e4, pct_low=0, pct_high=0){
	if (length(x@test_IDs) == 0){
		stop("No test IDs.")
	}

	x@labels <- rep(0, length(x@test_IDs))
	names(x@labels) <- x@test_IDs

	dc <- Matrix::colSums(x@counts[,x@test_IDs])
	n_low <- floor(length(dc)*(pct_low))
	n_high <- floor(length(dc)*(pct_high))

	if (pct_low > 0){
		n_low <- floor(length(dc)*(pct_low))
		botn <- names(dc)[order(dc, decreasing=FALSE)[1:n_low]]
		x@labels[botn] <- 2
	}

	if (pct_high > 0){
		n_high <- floor(length(dc)*(pct_high))
		topn <- names(dc)[order(dc, decreasing=TRUE)[1:n_high]]
		x@labels[topn] <- 1
	}
	return(x)
}

#' Read 10X output
#'
#' Read output from 10X into a sparse matrix. Given path, reads 
#' files path/matrix.MM path/genes.tsv and path/barcodes.tsv if v2, 
#' or path/matrix.mtx.gz path/features.tsv.gz path/barcodes.tsv.gz
#'
#' @param path Character. File path prefix to 10X output.
#'
#' @return expr Sparse matrix. Expression counts from 10X cell ranger
#'  with cells in the columns and genes in the rows.
#'
#' @import Matrix
#' @export
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
