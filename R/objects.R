
#' EMO
#' @name EMO-class
#' @rdname EMO-class
#' @exportClass EMO
EMO <- setClass(Class = "EMO", 
				slots = c(Z = "matrix", 
						  mu = "list", 
						  sgma = "list", 
						  tau = "vector", 
						  llks = "vector", 
						  n_iter = "numeric"))

#' DE
#' @name DE-class
#' @rdname DE-class
#' @exportClass DE
DE <- setClass(Class = "DE",
			   slots = c(low_means = "numeric",
			   			 high_means = "numeric", 
			   			 low_n = "numeric", 
			   			 high_n = "numeric", 
			   			 deg = "character", 
			   			 deg_low = "character", 
			   			 deg_high = "character", 
			   			 log2fc = "numeric"))

#' DIEM
#' @name DIEM-class
#' @rdname DIEM-class
#' @exportClass DIEM
DIEM <- setClass(Class = "DIEM", 
				slots = c(counts = "dgCMatrix", 
						  norm = "ANY", 
						  pi = "data.frame", 
						  labels = "vector", 
						  emo = "EMO"))

#' SCE
#'
#' Single Cell Expression object
#' @name SCE-class
#' @rdname SCE-class
#' @exportClass SCE
SCE <- setClass(Class = "SCE", 
				slots = c(counts = "dgCMatrix", 
						  norm = "ANY",
						  de_cutpoint = "numeric", 
						  min_counts = "numeric", 
						  fix_counts = "numeric", 
						  de = "DE", 
						  diem = "DIEM", 
						  gene_info = "data.frame",
						  dropl_info = "data.frame", 
                          name = "character"))

#' @method dim SCE
#' @export
dim.SCE <- function(x){
	return(dim(x@counts))
}

#' Fill information from raw counts
#'
#' @param x SCE. SCE object with x@counts data present
#'
#' @importFrom Matrix Matrix rowMeans rowSums colSums
#' @export
fill_counts <- function(x){
	x@gene_info <- x@gene_info[rownames(x@counts), , drop=FALSE]; rownames(x@gene_info) <- rownames(x@counts)
	x@dropl_info <- x@dropl_info[colnames(x@counts), , drop=FALSE]; rownames(x@dropl_info) <- colnames(x@counts)
	x@gene_info[,"mean"] <- Matrix::rowMeans(x@counts)
	x@gene_info[,"total_counts"] <- Matrix::rowSums(x@counts)
	x@gene_info[,"n_cells"] <- Matrix::rowSums(x@counts > 0)
	x@dropl_info[,"total_counts"] <- Matrix::colSums(x@counts)
	x@dropl_info[,"n_genes"] <- Matrix::colSums(x@counts > 0)
	return(x)
}

#' Create SCE object from sparse matrix
#'
#' @param x Sparse matrix. Raw expression counts from single-cell experiment, genes in rows, droplets in columns
#' @param name Character. Name for SCE object.
#'
#' @importFrom Matrix Matrix
#' @return SCE object
#' @export
create_SCE <- function(x, name="SCE"){
	if (is.null(rownames(x))) rownames(x) <- paste0("Gene", seq_len(nrow(x)))
	if (is.null(colnames(x))) colnames(x) <- paste0("Drop", seq_len(ncol(x)))
	# Get sparse matrix
	if (!inherits(x, what="dgCMatrix")){
		sce <- SCE(counts = Matrix(x, sparse=TRUE))
	}else{
		sce <- SCE(counts = x)
	}
	sce@gene_info <- data.frame(row.names=rownames(sce@counts))
	sce@dropl_info <- data.frame(row.names=colnames(sce@counts))
	# rownames(sce@gene_info) <- rownames(sce@counts)
	sce <- fill_counts(sce)
	keep <- sce@dropl_info[,"total_counts"] > 0
	sce@counts <- sce@counts[,keep]
	sce@dropl_info <- sce@dropl_info[keep,]
	sce <- fill_counts(sce)
	sce@name <- name
	return(sce)
}

#' Convert SCE object to Seurat
#'
#' Convert SCE object to Seurat. if \code{targets} is true (default), outputs only droplets that are 
#' called targets. In addition, outputs meta information calculated, such as bg_score. Additional arguments to 
#' this function are passed onto \code{\link[Seurat]{CreateSeuratObject}}. Common arguments include 
#' "min.cells = 3" and "min.features = 200".
#'
#' @param x SCE object
#' @param targets Boolean. If \code{TRUE}, then output droplets that are specified as targets
#' @param meta Boolean. Place information in dropl_info slot into meta.data in Seurat object
#' @param ... Arguments to \code{\link[Seurat]{CreateSeuratObject}}, such as \code{project} for project name
#'
#' @importFrom Seurat CreateSeuratObject
#' @return A Seurat object
#' @export
convert_to_seurat <- function(x, targets=TRUE, meta=TRUE, ...){
	if (!requireNamespace("Seurat", quietly = TRUE)) {
		stop("Package \"Seurat\" needed for this function to work. Please install it.",
		    	       call. = FALSE)
	}
	if (targets) keep <- rownames(x@dropl_info)[ grep("Nucleus", x@dropl_info[,"Call"]) ]
	else keep <- rownames(x@dropl_info)
	if (meta) meta.data <- x@dropl_info[keep,,drop=FALSE]
	else meta.data <- NULL
	seur <- Seurat::CreateSeuratObject(counts=x@counts[,keep], meta.data=meta.data, ...)
	return(seur)
}

