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

#' Normalize raw counts for droplets being classified.
#'
#' Counts are normalized by dividing by total counts per droplet, and log1p-transformed. 
#' Optionally, counts can be multiplied by a scaling factor \code{scale_factor}.
#'
#' @param x SCE.
#' @param scale_factor Numeric. A scaling factor to multiply values after division by total count.
#'
#' @return SCE object
#' @importFrom Matrix Diagonal colSums
#' @export
normalize <- function(x, 
					  scale_factor=1e4,
					  logt=FALSE){
	expr <- x@diem@counts
	expr <- divide_by_colsum(expr)
	expr <- expr * scale_factor
	if (logt) expr <- log1p(expr)
	x@diem@norm <- expr
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
subset_for_em <- function(x, log_base=10, top_n=NULL, min_counts=NULL){
	dc <- Matrix::colSums(x@counts)
	dc_max <- max(dc)
	if (is.null(top_n) & is.null(min_counts)){
		x@min_counts <- dc_max/(log_base^2)
	} else {
		# If user set parameters
		if ( ( !is.null(top_n) ) & ( !is.null(min_counts) ) ){
			stop("Set only one of top_n or min_counts")
		}
		if ( is.null(top_n) ){
			x@min_counts <- min_counts
		} else {
			top_n_ix <- order(dc, decreasing=TRUE)[top_n]
			x@min_counts <- dc[top_n_ix]
		}
	}
	# Subset counts matrix
	x@diem <- DIEM(counts = x@counts[,dc > x@min_counts])
	return(x)
}

#' Label upper count droplets 
#'
#' Fix labels of upper count droplets to target according to total counts in a droplet. The threshold 
#' is set by taking the log of the maximum droplet count, substracting by 1, and exponentiating by 
#' base of the log. In other words, this divides the max count by \code{log_base} and sets labels of 
#' droplets with counts greater than this to target.
#'
#' @param x SCE. SCE object.
#' @param log_base Numeric. The base of log to take
#'
#' @return An SCE object
#' @importFrom Matrix colSums
#' @export
set_labels <- function(x, log_base=2){
	x@diem@labels <- rep(0, ncol(x@diem@counts))
	dc <- Matrix::colSums(x@diem@counts)
	dc_max <- max(dc)
	x@fix_counts <- dc_max/log_base
	x@diem@labels[ dc > x@fix_counts ] <- 1
	return(x)
}

#' Read 10X output
#'
#' Read output from 10X into a sparse matrix. Given path, reads 
#' files path/matrix.MM path/genes.tsv, and path/barcodes.tsv
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
