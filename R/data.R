#' Divide elements of a column in a sparse matrix by its column sum
#'
#' @param x Sparse Matrix
#'
#' @return A Sparse Matrix with columns sum to 1.
#' @importFrom Matrix Diagonal colSums
#' @export
divide_by_colsum <- function(x){
	cs <- Matrix::colSums(x = x)
	d <- Matrix::Diagonal(x = 1/cs)
	x <- x %*% d
	return(x)
}

#' Normalize simulated and candidate raw count data
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
					  scale_factor=1){
	expr <- x@sim@counts
	expr <- divide_by_colsum(expr)
	expr <- expr * scale_factor
	# expr <- log1p(expr)
	x@sim@norm <- expr
	return(x)
}

#' Get expression matrix of candidate counts
#'
#' @param x SCE.
#'
#' @return \code{\link[Matrix]{Matrix}} object
#' @export
get_candd_counts <- function(x){
	dc <- Matrix::colSums(x@counts)
	high_droplets <- dc > x@min_counts
	return(x@counts[, high_droplets])
}

#' Set counts threshold for testing droplets
#'
#' Sets the threshold for testing droplets. Only droplets with a total number of counts of at least
#' \code{min_counts} are included in the semi-supervised EM and assigned as either background or target. 
#' By default, this number is calculated by dividing the maximum droplet count by 100 in order to test 
#' the top 2 orders of magnitude in the data. 
#' The minimum number of counts can be set using the \code{top_n} parameter instead, so that only the 
#' \code{top_n} droplets ranked by total counts are used in testing. Another way is to directly set the 
#' minimum counts through \code{min_counts}.
#'
#' @param x SCE. SCE object.
#' @param top_n Numeric. Only the top \code{top_n} droplets ranked by total counts are used for testing.
#' @param min_counts Integer. Only droplets with at least this number of total counts are used for testing.
#'
#' @return An SCE object with \code{min_counts} slot set
#' @importFrom Matrix colSums
#' @export
set_min_counts <- function(x, top_n=NULL, min_counts=NULL){
	dc <- Matrix::colSums(x@counts)
	if (is.null(top_n) & is.null(min_counts)){
		x@min_counts <- max(dc)/100
		return(x)
	}

	# If user set parameters
	if ( ( !is.null(top_n) ) & ( !is.null(min_counts) ) ){
		stop("Set only one of top_n or min_counts")
	}
	if ( !is.null(top_n) ){
		x@min_counts <- min_counts
	} else {
		top_n_ix <- order(dc, decreasing=TRUE)[top_n]
		x@min_counts <- max(min_counts, dc[top_n_ix])
	}
	return(x)
}


#' Subset droplets in SCE object
#'
#' Given thresholds min_c and max_c, subset droplets by read counts
#' Keep droplets where the count is greater than or 
#' equal to min_c and less than or equal to max_c.
#'
#' @param x SCE object. Must contain data in counts slot
#' @param min_c Numeric. Minimum row sum
#' @param max_c Numeric. Maximum row sum
#'
#' @return x SCE object subsetted
#' @export
subset_dropls <- function(x, min_c=-Inf, max_c=Inf){

	keep = x@dropl_info[,"total_counts"] > min_c & x@dropl_info[,"total_counts"] <= max_c
	for (i in c("counts", "norm")){
		m <- slot(x, i)
		if (length(m) == 0) next
		slot(x, i) <- slot(x, i)[,keep]
	}
	x@dropl_info <- x@dropl_info[keep,,drop=FALSE]
	if ( ! all(is.na(x@pcs)) )  x@pcs$x <- x@pcs$x[keep,]
	return(x)
}

#' Subset to top n droplets in SCE object
#'
#' @param x SCE object. Must contain data in counts slot
#' @param n Integer. Number of droplets to subset
#'
#' @return x SCE object subsetted
#' @export
subset_n_dropls <- function(x, n){
	ranked <- order(x@dropl_info[,"total_counts"], decreasing=TRUE)
	if (n > length(ranked)) n <- length(ranked)
	ix <- ranked[1:n]
	for (i in c("counts", "norm")){
		m <- slot(x, i)
		if (length(m) == 0) next
		slot(x, i) <- slot(x, i)[,ix]
	}
	x@dropl_info <- x@dropl_info[ix,,drop=FALSE]
	if ( ! all(is.na(x@pcs$x)) )  x@pcs$x <- x@pcs$x[ix,]
	return(x)
}

#' Subset genes in SCE object
#'
#' Given thresholds min_c and max_c, subset genes by mean counts.
#' Keep genes where the count is greater than or 
#' equal to min_m and less than or equal to max_m.
#'
#' @param x SCE object. Must contain data in counts slot
#' @param min_c Numeric. Minimum gene counts
#' @param max_c Numeric. Maximum gene counts
#'
#' @return x SCE object subsetted
#' @export
subset_genes <- function(x, min_c=-Inf, max_c=Inf, genes=NULL){
	keep <- (x@gene_info[,"total_counts"] > min_c) & 
			(x@gene_info[,"total_counts"] <= max_c)
	keep <- keep & (rownames(x@gene_info) %in% genes)
	for (i in c("counts", "norm")){
		m <- slot(x, i)
		if (length(m) == 0) next
		slot(x, i) <- slot(x, i)[keep,]
	}
	# x@gene_info <- x@gene_info[keep,,drop=FALSE]
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
