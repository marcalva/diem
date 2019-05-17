#' Normalize raw count data
#'
#' Raw expression counts are divided by total counts per droplet, 
#' multiplied by scaling factor, and (optionally) log1p transformed (default). 
#' This normalizes the counts in 
#' \code{raw} slot and places them in the @norm slot. The scaling factor 
#' defaults to the median total read count across all droplets.
#'
#' @param x SCE. SCE 
#' @param logt Boolean. Log transform
#' @param sf Numeric. A scaling factor to multiply values after division by total count. If NULL, use the median of counts
#' @param verbose Boolean.
#'
#' @return An SCE object with normalized expression data
#' @importFrom Matrix Diagonal
#' @export
normalize <- function(x, logt=TRUE, sf=1, verbose=FALSE){
	if (verbose) cat("Normalizing gene expression\n")
	expr <- x@raw
	if (is.null(sf)) sf <- median(Matrix::colSums(expr))
	expr <- expr * sf
	d <- Matrix::Diagonal(x = 1/x@dropl_info[,"total_counts"])
	expr <- expr %*% d
	if (logt) expr <- log1p(expr)
	x@norm <- expr
	if (verbose) cat("Normalized gene expression\n")
	return(x)
}

#' Get variable genes
#'
#' Finds top \code{nf} variable genes by running a loess regression 
#' on log(variance) and log(mean). Residuals are used to rank 
#' features, and the top \code{nf} are output.
#'
#' @param x SCE. SCE object
#' @param min_mean Numeric. Minimum mean of counts to consider
#' @param nf Numeric. Output top nf variable features
#' @param lss Numeric. The span paramater for loess function.
#' @param verbose Boolean.
#'
#' @return SCE object with \code{vg} and \code{vg_info} slots filled for variable genes
#' and variable gene info. \code{vg_info} is a data frame with log mean, log var, 
#' fitted log var values, and residuals. \code{vg} is a vector of genes
#' @export
get_var_genes <- function(x, min_mean = 0.0001, nf = 2000, lss=0.3, verbose=FALSE){
	if (verbose) cat("Getting variable genes\n")
	if (min_mean <= 0){
		keep <- x@gene_info[,"mean"] > 0
	} else {
		keep <- x@gene_info[,"mean"] > min_mean
	}
	if (table(keep)["TRUE"] < nf){
		stop("Error: Specifying more variable features than there are features passing min_mean\n")
	}
	gene_means <- x@gene_info[,"mean"][keep]
	gene_names <- rownames(x@gene_info)[keep]

	log_mean <- log10(gene_means + 1)
	log_var <- log10(fast_varCPP(x@raw[keep,], gene_means) + 1)
	fit <- loess(log_var ~ log_mean, span=lss)
	rsd <- log_var - fit$fitted
	topi <- order(rsd, decreasing=TRUE)[1:nf]
	vg <- gene_names[topi]
	df <- data.frame(mean=log_mean, var=log_var,
					 fit=fit$fitted, rsd=rsd)
	rownames(df) <- gene_names
	x@vg <- vg
	x@vg_info <- df
	if (verbose) cat("Found variable genes\n")
	return(x)
}

#' Get DE genes between high and low read count droplets
#'
#' Find differentially expressed genes between low and high read count droplets. Sums
#' gene expression across all cells within the limits specified in the SCE object, 
#' divides by the total number of reads in each, and calculates the difference between 
#' the 2. The top \code{n_genes} genes sorted by absolute difference in proportion are 
#' designated as DE.
#'
#' @param x SCE. SCE object
#' @param n_genes Numeric. Number of genes to output as differentially expressed between background and candidate
#' @param verbose Boolean.
#'
#' @return SCE object with \code{gene_info[,"diff_prop"]} and \code{deg} slots filled. 
#' The \code{diff_prop} column contains the background - candidate difference in proportions.
#' The \code{deg} slot contains the names of the genes output as DE.
#' @export
get_de_genes <- function(x, n_genes=2000, verbose=FALSE){
	if (verbose) cat("Getting genes differentially expressed between high and low count droplets\n")
	dc <- x@dropl_info[,"total_counts"]
	low_expr <- x@raw[, (dc > x@limits$min_bg_count) & (dc <= x@limits$max_bg_count)]
	high_expr <- x@raw[, (dc > x@limits$min_tg_count) & (dc <= x@limits$max_tg_count)]

	low_prop <- Matrix::rowSums(low_expr); low_prop <- low_prop/sum(low_prop)
	high_prop <- Matrix::rowSums(high_expr); high_prop <- high_prop/sum(high_prop)
	diff_prop <- low_prop - high_prop
	
	x@gene_info[,"diff_prop"] <- diff_prop
	de_genes <- names(sort(abs(diff_prop), decreasing=TRUE))[1:n_genes]
	x@deg <- de_genes
	if (verbose) cat("Found DE genes\n")
	return(x)
}

#' Subset droplets in SCE object
#'
#' Given thresholds min_c and max_c, subset droplets by read counts
#' Keep droplets where the count is greater than or 
#' equal to min_c and less than or equal to max_c.
#'
#' @param x SCE object. Must contain data in raw slot
#' @param min_c Numeric. Minimum row sum
#' @param max_c Numeric. Maximum row sum
#'
#' @return x SCE object subsetted
#' @export
subset_dropls <- function(x, min_c=-Inf, max_c=Inf){

	keep = x@dropl_info[,"total_counts"] > min_c & x@dropl_info[,"total_counts"] <= max_c
	for (i in c("raw", "norm")){
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
#' @param x SCE object. Must contain data in raw slot
#' @param n Integer. Number of droplets to subset
#'
#' @return x SCE object subsetted
#' @export
subset_n_dropls <- function(x, n){
	ranked <- order(x@dropl_info[,"total_counts"], decreasing=TRUE)
	if (n > length(ranked)) n <- length(ranked)
	ix <- ranked[1:n]
	for (i in c("raw", "norm")){
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
#' @param x SCE object. Must contain data in raw slot
#' @param min_c Numeric. Minimum gene counts
#' @param max_c Numeric. Maximum gene counts
#'
#' @return x SCE object subsetted
#' @export
subset_genes <- function(x, min_c=-Inf, max_c=Inf, genes=NULL){
	keep <- (x@gene_info[,"total_counts"] > min_c) & 
			(x@gene_info[,"total_counts"] <= max_c)
	keep <- keep & (rownames(x@gene_info) %in% genes)
	for (i in c("raw", "norm")){
		m <- slot(x, i)
		if (length(m) == 0) next
		slot(x, i) <- slot(x, i)[keep,]
	}
	x@gene_info <- x@gene_info[keep,,drop=FALSE]
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
