#' Divide elements of a column in a sparse matrix by its column sum
#'
#' @param x Sparse Matrix
#'
#' @return A Sparse Matrix with columns that sum to 1.
#' @importFrom Matrix Diagonal colSums
#' @export
divide_by_colsum <- function(x){
	cs <- Matrix::colSums(x = x)
	d <- Matrix::Diagonal(x = 1/cs)
	x <- x %*% d
	return(x)
}

#' Subset DE to genes in keep
#'
#' @param x SCE object.
#' @param keep Character vector. Names of genes to keep in DE object
#'
#' @return An SCE object
#' @export
subset_DE <- function(x, keep){
	x@de@low_means <- x@de@low_means[intersect(names(x@de@low_means), keep)]
	x@de@high_means <- x@de@high_means[intersect(names(x@de@high_means), keep)]
	x@de@deg <- intersect(x@de@deg, keep)
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
#' @param scale. Boolean. Whether to scale genes to mean 0 and variance 1.
#' @param de_genes Boolean. Whether to normalize only the DE genes.
#' @param scale_factor Numeric. A scaling factor to multiply values after division by total count.
#' @param logt Boolean. Whether to log1p transform expression values.
#'
#' @return SCE object
#' @importFrom Matrix rowMeans
#' @export
normalize <- function(x, 
					  scale.=TRUE,
					  genes=NULL, 
					  scale_factor=1,
					  logt=TRUE,
					  verbose=FALSE){
	if (length(x@diem@counts) == 0){
		stop("Specify test set with set_test_set function")
	}
	if (verbose) cat("Normalizing\n")

	expr <- x@diem@counts
	if (!is.null(genes)) expr <- expr[genes,]
	expr <- divide_by_colsum(expr)
	expr <- expr * scale_factor
	if (logt) expr <- log1p(expr)

	if (scale.) {
		genes <- rownames(expr); droplets <- colnames(expr)
		means <- Matrix::rowMeans(expr); names(means) <- genes
		vars <- fast_varCPP(x=expr, mu=means); names(vars) <- genes

		# Remove zero variance and zero mean genes
		keep <- genes[vars > 0 & means > 0]
		means <- means[keep]
		vars <- vars[keep]
		expr <- expr[keep,]
		x <- subset_DE(x, keep)

		# Scale
		expr <- fast_row_scaleCPP(x=expr, mu=means, sigma_sq=vars)
		rownames(expr) <- keep
		colnames(expr) <- droplets
	}
	x@diem@norm <- expr
	if (verbose) cat("Normalized\n")

	return(x)
}

#' Subset droplets for EM testing
#'
#' Sets the threshold for testing droplets and subsets counts matrix. Only droplets with a total number 
#' of counts of at least \code{min_counts} are included in the semi-supervised EM and assigned as either 
#' background or target. By default, this number is calculated by taking 2 orders of magnitude lower than
#' the maximum droplet count in the data. The order of magnitude is set by \code{log_base}. In other words, 
#' this divides the max droplet count by the square of \code{log_base} and sets this as the minimum count. 
#' The minimum number of counts can be set using the \code{top_n} parameter instead, so that only the 
#' \code{top_n} droplets ranked by total counts are used in testing. Another way is to directly set the 
#' minimum counts through \code{min_counts}. The count matrix is subset to droplets having total counts greater 
#' than this threshold.
#'
#' @param x SCE. SCE object.
#' @param log_base Numeric. Base of log to take.
#' @param top_n Numeric. Only the top \code{top_n} droplets ranked by total counts are classified.
#' @param min_counts Integer. Only droplets with at least this number of total counts are classified.
#'
#' @return An SCE object.
#' @importFrom Matrix colSums
#' @export
set_test_set <- function(x, top_n=1e4, log_base=NULL, min_counts=NULL){
	dc <- Matrix::colSums(x@counts)
	dc_max <- max(dc)
	if (is.null(top_n) & is.null(min_counts)){
		min_counts <- dc_max/(log_base^2)
	} else {
		# If user set parameters
		if ( ( !is.null(top_n) ) & ( !is.null(min_counts) ) ){
			stop("Set only one of top_n or min_counts")
		}
		if ( is.null(min_counts) ){
			top_n_ix <- order(dc, decreasing=TRUE)[top_n]
			min_counts <- as.numeric(dc[top_n_ix])
		}
	}
	# Subset counts matrix and remove 0 mean genes
	expr <- x@counts[, dc >= min_counts]
	x@test_droplets <- colnames(expr)
	x@diem <- DIEM(counts = expr)
	return(x)
}

#' Fix labels of lower and upper count droplets 
#'
#' Fix labels of the upper and lower tails droplets to nuclei and background, respectively.
#' The top and bottom \code{pct} are then fixed during EM.
#'
#' @param x SCE. SCE object.
#' @param log_base Numeric. The base of log to take
#'
#' @return An SCE object
#' @importFrom Matrix colSums
#' @export
fix_tails <- function(x, pct=0.05){
	x@diem@labels <- rep(0, ncol(x@diem@counts))
	names(x@diem@labels) <- colnames(x@diem@counts)
	dc <- Matrix::colSums(x@diem@counts)
	names(dc) <- colnames(x@diem@counts)
	n <- floor(length(dc)*pct)
	topn <- names(dc)[order(dc, decreasing=TRUE)[1:n]]
	botn <- names(dc)[order(dc, decreasing=FALSE)[1:n]]
	x@diem@labels[topn] <- 1
	x@diem@labels[botn] <- 2
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
