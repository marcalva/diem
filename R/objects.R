
#' SCE
#'
#' Single Cell Expression object
#' @name SCE-class
#' @rdname SCE-class
#' @exportClass SCE
SCE <- setClass(Class = "SCE", 
				slots = c(raw = "dgCMatrix", 
						  norm = "ANY", 
						  gene_info = "data.frame",
						  dropl_info = "data.frame", 
						  emo = "list", 
                          vg = "vector", 
                          vg_info = "data.frame", 
                          de_genes = "vector", 
                          bins = "matrix", 
                          bins_top = "vector", 
                          pcs = "ANY",
                          name = "character"))

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
	x@gene_info <- x@gene_info[rownames(x@raw), , drop=FALSE]; rownames(x@gene_info) <- rownames(x@raw)
	x@dropl_info <- x@dropl_info[colnames(x@raw), , drop=FALSE]; rownames(x@dropl_info) <- colnames(x@raw)
	x@gene_info[,"mean"] <- Matrix::rowMeans(x@raw)
	x@gene_info[,"total_counts"] <- Matrix::rowSums(x@raw)
	x@gene_info[,"n_cells"] <- Matrix::rowSums(x@raw > 0)
	x@dropl_info[,"total_counts"] <- Matrix::colSums(x@raw)
	x@dropl_info[,"n_genes"] <- Matrix::colSums(x@raw > 0)
	return(x)
}

#' Create SCE object from sparse matrix
#'
#' @param x Sparse matrix. Raw expression counts from single-cell experiment, genes in rows, droplets in columns
#'
#' @importFrom Matrix Matrix
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
	sce@gene_info <- data.frame(row.names=rownames(sce@raw))
	sce@dropl_info <- data.frame(row.names=colnames(sce@raw))
	# rownames(sce@gene_info) <- rownames(sce@raw)
	sce <- fill_counts(sce)
	keep <- sce@dropl_info[,"total_counts"] > 0
	sce@raw <- sce@raw[,keep]
	sce@dropl_info <- sce@dropl_info[keep,]
	sce <- fill_counts(sce)
	return(sce)
}

	
#' Convert SCE object to Seurat
#'
#' @param x SCE object
#' @param targets Boolean. If \code{TRUE}, then output droplets that are specified as targets
#' @param meta Boolean. Place information in dropl_info slot into meta.data in Seurat object
#' @param ... Arguments to \code{\link[Seurat]{CreateSeuratObject}}, such as \code{project} for project name
#'
#' @return A Seurat object
#' @importFrom Seurat CreateSeuratObject
#' @export
convert_to_seurat <- function(x, targets=TRUE, meta=TRUE, ...){
	if (targets) keep <- x@dropl_info[,"Target"]
	else keep <- rep(TRUE, nrow(x@dropl_info))
	if (meta) meta.data <- x@dropl_info[keep,,drop=FALSE]
	else meta.data <- NULL
	seur <- CreateSeuratObject(counts=x@raw[,keep], meta.data=meta.data, ...)
	return(seur)
}

