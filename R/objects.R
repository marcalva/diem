
#' SCE
#'
#' Single Cell Expression object
#' @name SCE-class
#' @rdname SCE-class
#' @exportClass SCE
SCE <- setClass(Class = "SCE", 
				slots = c(raw = "dgCMatrix", 
						  norm = "ANY", 
						  bg_info = "list", 
						  emo = "list", 
						  targets = "vector", 
						  other_data = "data.frame", 
                          gene_means = "vector", 
                          n_genes = "vector", 
                          gene_counts = "vector", 
                          dropl_counts = "vector", 
                          vg = "vector", 
                          vg_info = "data.frame", 
                          pcs = "ANY"))

#' @method dim SCE
#' @export
dim.SCE <- function(x){
	return(dim(x@raw))
}

#' Fill information from raw counts
#'
#' @param x SCE. SCE object with x@raw data present
#'
#' @importFrom Matrix Matrix rowSums colSums
#' @export
fill_counts <- function(x){
	x@gene_means <- Matrix::rowMeans(x@raw)
	x@gene_counts <- Matrix::rowSums(x@raw)
	x@dropl_counts <- Matrix::colSums(x@raw)
	x@n_genes <- Matrix::colSums(x@raw > 0)
	return(x)
}

#' Create SCE object from sparse matrix
#'
#' @param x Sparse matrix. Raw expression counts from single-cell experiment
#'
#' @return SCE object
#' @export
create_SCE <- function(x){
	if (is.null(rownames(x))) rownames(x) <- paste0("Gene", seq_len(nrow(x)))
	if (is.null(colnames(x))) colnames(x) <- paste0("Drop", seq_len(ncol(x)))
	# Get sparse matrix
	if (!inherits(x, what="dgCMatrix")){
		sce <- SCE(raw = Matrix(x, sparse=TRUE))
	}else{
		sce <- SCE(raw = x)
	}
	rownames(sce@raw) <- as.character(rownames(sce@raw))
	sce <- fill_counts(sce)
	keep <- sce@dropl_counts > 0
	sce@raw <- sce@raw[,keep]
	sce@dropl_counts <- sce@dropl_counts[keep]
	sce@n_genes <- sce@n_genes[keep]
	return(sce)
}

	
#' Convert SCE object to Seurat
#'
#' @param x SCE object
#' @param targets Boolean. If \code{TRUE}, then output droplets that are specified as targets
#' @param ... Arguments to \code{\link[Seurat]{CreateSeuratObject}}, such as \code{project} for project name
#'
#' @return A Seurat object
#' @importFrom Seurat CreateSeuratObject
#' @export
convert_to_seurat <- function(x, targets=TRUE, ...){
	if(targets) {
		seur <- CreateSeuratObject(counts=x@raw[,x@targets], ...)
	} else {
		seur <- CreateSeuratObject(counts=x@raw, ...)
	}
	return(seur)
}

