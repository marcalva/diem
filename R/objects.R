
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
						  limits = "list",
						  emo = "list", 
                          vg = "vector", 
                          vg_info = "data.frame", 
                          deg = "vector", 
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
	sce@name <- name
	return(sce)
}

#' Set count and gene limits for background and target
#'
#' Set the minimum and maximum number of counts for a droplet to be considered
#' either background or droplet. Also sets the limit of minimum or maximum number 
#' of detected genes that a droplet can be considered a target. On top of setting 
#' the number explicitly, the min_tg_count and max_bg_count can be set by ranking 
#' of droplet total counts. For example, you can set droplets above the 10,000th 
#' ranked droplet as candidate targets and the all droplets under the 20,000th 
#' ranked droplet as background. The parameters top_n_tg and top_n_bg override 
#' the numeric specifications of max_bg_count and min_tg_count.
#'
#' @param x SCE object
#' @param min_bg_count Numeric. Minimum number of counts for a droplet to be considered background
#' @param max_bg_count Numeric. Maximum number of counts for a droplet to be considered background
#' @param min_tg_count Numeric. Minimum number of counts for a droplet to be considered target
#' @param max_tg_count Numeric. Maximum number of counts for a droplet to be considered target
#' @param min_tg_gene Numeric. Minimum number of genes detected for a droplet to be considered target
#' @param max_tg_gene Numeric. Maximum number of genes detected for a droplet to be considered target
#' @param top_n_tg Numeric. Set the top n ranked by counts as the threshold for considering targets
#' @param top_n_bg Numeric. Set the top n ranked by counts as the threshold for considering background
#'
#' @return SCE object
#' @export
set_limits <- function(x, 
					   min_bg_count=0,
					   max_bg_count=150,
					   min_tg_count=199,
					   max_tg_count=Inf,
					   min_tg_gene=199,
					   max_tg_gene=Inf,
					   top_n_tg=NULL,
					   top_n_bg=NULL){
	limits <- list(min_bg_count=min_bg_count,
				   max_bg_count=max_bg_count,
				   min_tg_count=min_tg_count,
				   max_tg_count=max_tg_count,
				   min_tg_gene=min_tg_gene,
				   max_tg_gene=max_tg_gene)
	if (!is.null(top_n_tg)){
		ranks <- order(x@dropl_info[,"total_counts"], decreasing=TRUE)
		limits$min_tg_count <- x@dropl_info[ranks[top_n_tg],"total_counts"]
	}
	if (!is.null(top_n_bg)){
		ranks <- order(x@dropl_info[,"total_counts"], decreasing=TRUE)
		limits$max_bg_count <- x@dropl_info[ranks[top_n_bg],"total_counts"]
	}
	if (limits$min_tg_count < limits$max_bg_count){
		stop("Minimum count threshold for targets must be greater than maximum count threshold for background")
	}
	x@limits <- limits
	return(x)
}

#' Convert SCE object to Seurat
#'
#' @param x SCE object
#' @param targets Boolean. If \code{TRUE}, then output droplets that are specified as targets
#' @param meta Boolean. Place information in dropl_info slot into meta.data in Seurat object
#' @param ... Arguments to \code{\link[Seurat]{CreateSeuratObject}}, such as \code{project} for project name
#'
#' @return A Seurat object
#' @export
convert_to_seurat <- function(x, targets=TRUE, meta=TRUE, ...){
	if (!requireNamespace("Seurat", quietly = TRUE)) {
		stop("Package \"Seurat\" needed for this function to work. Please install it.",
		    	       call. = FALSE)
	}
	if (targets) keep <- x@dropl_info[,"Target"]
	else keep <- rep(TRUE, nrow(x@dropl_info))
	if (meta) meta.data <- x@dropl_info[keep,,drop=FALSE]
	else meta.data <- NULL
	seur <- Seurat::CreateSeuratObject(counts=x@raw[,keep], meta.data=meta.data, ...)
	return(seur)
}

