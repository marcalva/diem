
#' @importClassesFrom Matrix dgCMatrix
setClassUnion("any_matrix", c("matrix", "dgCMatrix"))

#' DIEM
#' @name DIEM-class
#' @rdname DIEM-class
#' @exportClass DIEM
DIEM <- setClass(Class = "DIEM", 
				 slots = c(emo = "list", 
						   pi = "matrix", 
						   PP = "matrix", 
						   calls = "character", 
						   converged = "logical"))

# emo contains the EM output of each iteration.
# Each element of this list contains
#    Z matrix  n x k of log likelihoods of each data point
#	 Mu means
#	 Mc Mixture Coefficient
#    loglks is a list of log likelihood value
#    conv.iter Number of iterations to convergene

#' SCE
#'
#' Single Cell Expression object
#' @name SCE-class
#' @rdname SCE-class
#' @exportClass SCE
SCE <- setClass(Class = "SCE", 
				slots = c(counts = "any_matrix", 
						  norm = "any_matrix", 
						  min_counts = "numeric", 
						  top_thresh = "numeric", 						  
						  labels = "numeric", 
						  diem = "DIEM", 
						  pp_thresh = "numeric", 
						  gene_info = "data.frame",
						  dropl_info = "data.frame", 
                          name = "character"))

#' @method dim SCE
#' @export
dim.SCE <- function(x){
	return(dim(x@counts))
}

#' @method rownames SCE
#' @export
rownames.SCE <- function(x){
	return(rownames(x@counts))
}

#' @method colnames SCE
#' @export
colnames.SCE <- function(x){
	return(colnames(x@counts))
}

#' Fill information from raw counts
#'
#' @param x An SCE object with x@counts data present
#'
#' @importFrom Matrix Matrix rowMeans rowSums colSums
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

#' Create an SCE object from a sparse matrix
#'
#' @param x A sparse matrix consisting of raw expression counts from a single-cell experiment, 
#'  with genes in the rows and droplets in the columns.
#' @param name Character name for the SCE object.
#'
#' @importFrom Matrix Matrix
#' @return SCE object
#' @export
#' @examples
#' counts <- read_10x("mouse_nuclei_2k/raw_gene_bc_matrices/mm10/")
#' mb_sce <- create_SCE(x=counts, name="Mouse Brain")
create_SCE <- function(x, name="SCE"){
	if (is.null(rownames(x))) rownames(x) <- paste0("Gene", seq_len(nrow(x)))
	if (is.null(colnames(x))) colnames(x) <- paste0("Drop", seq_len(ncol(x)))
	
	# Get sparse matrix
	if (!inherits(x, what="dgCMatrix")){
		sce <- SCE(counts = Matrix(x, sparse=TRUE))
	}else{
		sce <- SCE(counts = x)
	}

	rownames(sce@counts) <- make.unique(rownames(sce@counts))
	colnames(sce@counts) <- make.unique(colnames(sce@counts))
	
	sce@gene_info <- data.frame(row.names=rownames(sce@counts))
	sce@dropl_info <- data.frame(row.names=colnames(sce@counts))
	sce <- fill_counts(sce)

	keep <- sce@dropl_info[,"total_counts"] > 0
	sce@counts <- sce@counts[,keep]
	sce@dropl_info <- sce@dropl_info[keep,]

	sce@name <- name

	return(sce)
}

#' Convert an SCE object to Seurat
#'
#' Convert an SCE object to a Seurat object. if \code{targets} is true (default), output only droplets that are 
#' called as not debris. If \code{meta} is TRUE, then output meta data from dropl_info to the meta.data 
#' slot in the Seurat object. Additional functions \code{...} to this function are passed onto 
#' \code{\link[Seurat]{CreateSeuratObject}}. Common arguments include \code{min.cells = 3} and 
#' \code{min.features = 200}.
#'
#' @param x An SCE object.
#' @param targets Logical indicating whether to remove droplets called as debris. Default is TRUE.
#' @param meta Logical that indicates whether to place the data from dropl_info into meta.data in the 
#'  resulting Seurat object. Default is TRUE.
#' @param ... Arguments to \code{\link[Seurat]{CreateSeuratObject}}, such as \code{project} for project name.
#'
#' @return A Seurat object
#' @export
#' @examples
#' mm_seur <- convert_to_seurat(x=mb_sce, min.features = 200, min.cells = 3, project=mb_sce@name)
convert_to_seurat <- function(x, targets=TRUE, meta=TRUE, ...){
	if (!requireNamespace("Seurat", quietly = TRUE)) {
		stop("Package \"Seurat\" needed for convert_to_seurat. Please install.",
		    	       call. = FALSE)
	}

	if (targets){
		if (length(x@diem@calls) == 0) stop("Run DIEM before converting to Seurat with targets=TRUE.")
		keep <- get_clean_ids(x)
	}
	else keep <- rownames(x@dropl_info)

	if (meta) meta.data <- x@dropl_info[keep,,drop=FALSE]
	else meta.data <- NULL

	seur <- Seurat::CreateSeuratObject(counts=x@counts[,keep], meta.data=meta.data, ...)
	return(seur)
}
