
#' Normalize (total count and log) raw count data
#'
#' Raw expression counts are divided by total counts per droplet, 
#' multiplied by scaling factor, and log1p transformed. Same 
#' normalization style as Seurat. This normalizes the counts in 
#' \code{raw} slot and places them in the @norm slot.
#'
#' @param x SCE. SCE 
#' @param sf Numeric. A scaling factor to multiply values after division by total count
#'
#' @return A matrix with normalized expression data
#' @importFrom Matrix Diagonal
#' @export
normalize <- function(x, sf=1e4){
	expr <- x@raw
	expr <- expr * sf
	d <- Matrix::Diagonal(x = 1/x@dropl_counts)
	expr <- expr %*% d
	expr <- log1p(expr)
	x@norm <- expr
	return(x)
}

#' Get variable genes
#'
#' Finds top \code{nf} variable genes by running a loess regression 
#' on log(variance) and log(mean). Residuals are used to rank 
#' features, and the top \code{nf} are output.
#'
#' @param x SCE. SCE object
#' @param nf Numeric. Output top nf variable features
#' @param lss Numeric. The span paramater for loess function.
#'
#' @return SCE object with \code{vg} and \code{vg_info} slots filled for variable genes
#' and variable gene info. \code{vg_info} is a data frame with log mean, log var, 
#' fitted log var values, and residuals. \code{vg} is a vector of genes
#' @export
get_var_genes <- function(x, min_mean = 0.0001, nf = 2000, lss=0.3){
	if (min_mean <= 0){
		keep <- x@gene_means > 0
	} else {
		keep <- x@gene_means >= min_mean
	}
	if (table(keep)["TRUE"] < nf){
		stop("Error: Specifying more variable features than there are features passing min_mean\n")
	}
	gene_means <- x@gene_means[keep]
	gene_names <- rownames(x@raw)[keep]
	# gene_means <- Matrix::rowMeans(x@raw)

	log_mean <- log10(gene_means + 1)
	log_var <- log10(fast_varCPP(x@raw[keep,], gene_means) + 1)
	# log_varR <- log10(apply(expr, 1, var) + 1)
	# print(table(log_var == log_varR))
	fit <- loess(log_var ~ log_mean, span=lss)
	rsd <- log_var - fit$fitted
	topi <- order(rsd, decreasing=TRUE)[1:nf]
	vg <- gene_names[topi]
	df <- data.frame(mean=log_mean, var=log_var,
					 fit=fit$fitted, rsd=rsd)
	rownames(df) <- gene_names
	x@vg <- vg
	x@vg_info <- df
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
#' @param slot_name Character. Slot of SCE subset, either raw or norm, or both
#'
#' @return x SCE object subsetted
#' @export
subset_dropls <- function(x, min_c=-Inf, max_c=Inf, slot_name="raw"){
	if ( is.null(x@dropl_counts) )
		x <- fill_counts(x)

	keep = x@dropl_counts >= min_c & x@dropl_counts <= max_c
	for (i in slot_name){
		slot(x, slot_name) <- slot(x, slot_name)[,keep]
	}
	x@dropl_counts <- x@dropl_counts[keep]
	if ( ! all(is.na(x@pcs)) )  x@pcs$x <- x@pcs$x[ix,]
	return(x)
}

#' Subset to top n droplets in SCE object
#'
#' @param x SCE object. Must contain data in raw slot
#' @param n Integer. Number of droplets to subset
#' @param slot_name Character. Slot of SCE to subset, either raw or norm, or both
#'
#' @return x SCE object subsetted
#' @export
subset_n_dropls <- function(x, n, slot_name="raw"){
	if ( is.null(x@dropl_counts) )
		x <- fill_counts(x)

	ranked <- order(x@dropl_counts, decreasing=TRUE)
	ix <- ranked[1:n]
	for (i in slot_name){
		slot(x, i) <- slot(x, i)[,ix]
	}
	x@dropl_counts <- x@dropl_counts[ix]
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
#' @param. slot_name Character. Data to subset, either raw or norm.
#'
#' @return x SCE object subsetted
#' @export
subset_genes <- function(x, min_c=-Inf, max_c=Inf, slot_name="raw"){
	if ( is.null(x@gene_counts) )
		x <- fill_counts(x)

	keep = x@gene_counts >= min_c & x@gene_counts <= max_c
	for (i in slot_name){
		slot(x, i) <- slot(x, i)[keep,]
	}
	x@gene_means <- x@gene_means[keep]
	x@gene_counts <- x@gene_counts[keep]
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
#' @importFrom Matrix readMM
#' @export
read_10x <- function(path){
	files <- list.files(path, full.names=TRUE)
	mtx_file <- grep("matrix", files, value=TRUE)
	genes_file <- grep("genes|features", files, value=TRUE)
	barcode_file <- grep("barcodes", files, value=TRUE)

	expr <- readMM(mtx_file)
	expr <- as(expr, "CsparseMatrix")

	barcode_names <- readLines(barcode_file)
	barcode_names <- sub("-.*", "", barcode_names)

	genes <- read.delim(genes_file, header = FALSE, stringsAsFactors = FALSE, sep = "\t")

	colnames(expr) <- barcode_names
	rownames(expr) <- make.unique(genes[,2])

	return(expr)
}
	

