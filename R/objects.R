#' @useDynLib diem
#' @import Rcpp
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom methods setClassUnion setOldClass
setClassUnion("any_matrix", c("matrix", "dgCMatrix"))

#' SCE
#'
#' The SCE class
#'
#' This class is used to store the expression data from a 
#' single-cell RNA-seq experiment. It is the class that is used by 
#' DIEM. 
#'
#' Single Cell Expression object
#' @name SCE-class
#' @rdname SCE-class
#' @exportClass SCE
SCE <- setClass(Class = "SCE", 
                slots = c(counts = "any_matrix", 
                          norm = "any_matrix", 
                          pcs = "matrix", 
                          test_set = "character", 
                          bg_set = "character", 
                          gene_data = "data.frame", 
                          droplet_data = "data.frame", 
                          kruns = "list", 
                          vg_info = "data.frame", 
                          vg = "character", 
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
    x@gene_data <- x@gene_data[rownames(x@counts), , drop=FALSE]; rownames(x@gene_data) <- rownames(x@counts)
    x@droplet_data <- x@droplet_data[colnames(x@counts), , drop=FALSE]; rownames(x@droplet_data) <- colnames(x@counts)
    x@gene_data[,"mean"] <- Matrix::rowMeans(x@counts)
    x@gene_data[,"total_counts"] <- Matrix::rowSums(x@counts)
    x@gene_data[,"n_cells"] <- Matrix::rowSums(x@counts > 0)
    x@droplet_data[,"total_counts"] <- Matrix::colSums(x@counts)
    x@droplet_data[,"n_genes"] <- Matrix::colSums(x@counts > 0)
    return(x)
}

#' Create an SCE object from a sparse matrix
#'
#' @param x A sparse matrix consisting of raw expression counts from a 
#'  single-cell RNA experiment, with genes in the rows and 
#'  droplets in the columns.
#' @param name An optional character name for the SCE object.
#'
#' @importFrom Matrix Matrix
#' @importFrom methods new
#' @return SCE object
#' @export
#' @examples
#' counts <- matrix(sample(c(0,1,2), 1000, replace=TRUE),  nrow=10, ncol=100)
#' mb_sce <- create_SCE(x=counts, name="Mouse Brain")
create_SCE <- function(x, name = "SCE"){
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

    sce@gene_data <- data.frame(row.names=rownames(sce@counts))
    sce@droplet_data <- data.frame(row.names=colnames(sce@counts))
    sce <- fill_counts(sce)

    keep <- sce@droplet_data[,"total_counts"] > 0
    sce@counts <- sce@counts[,keep]
    sce@droplet_data <- sce@droplet_data[keep,]

    sce@name <- name

    return(sce)
}

#' Return the droplet data from an SCE object
#'
#' Return the data frame stored in the slot \code{droplet_data}. This 
#' contains data such as number of counts and genes in each droplet, as 
#' well as some of the output from the filtering, such as whether the 
#' droplet is classified as debris or cell/nucleus. The parameter 
#' \code{min_counts} filters out droplets (rows) by removing those 
#' with counts below this number
#'
#' @param x An SCE object.
#' @param min_counts Minimum number of read counts a droplet must have to 
#'  be output.
#' @param type One of either 'all' (default), 'test', 'clean', or 'debris' 
#'  specifying how to subset the data frame to include only those droplets.
#'
#' @return A data frame
#'
#' @export
droplet_data <- function(x, min_counts = 1, type = "all"){
    keep <- x@droplet_data$total_counts >= min_counts
    if (type == "clean"){
        keep <- keep & x@droplet_data$Call == "Clean"
    }
    if (type == "test"){
        keep <- keep & (rownames(x@droplet_data) %in% x@test_set)
    }
    if (type == "debris"){
        keep <- keep & x@droplet_data$Call == "Debris"
    }
    return(x@droplet_data[keep,,drop=FALSE])
}

#' Return the gene data from an SCE object
#' 
#' Return the data frame stored in the slot \code{gene_data}. This 
#' contains data such as the number of droplets the gene is detected in. 
#' The parameter \code{min_droplets} removes genes (rows) taht are 
#' detected in less than this number of droplets.
#'
#' @param x An SCE object.
#' @param min_droplets Minimum number of droplets detected per gene
#'
#' @return A data frame
#'
#' @export
gene_data <- function(x, min_droplets = 1){
    keep <- x@gene_data$n_cells >= min_droplets
    return(x@gene_data[keep,,drop=FALSE])
}

#' Return raw counts
#'
#' Return the raw counts from an SCE object.
#' 
#' @param x An SCE object.
#' 
#' @return Sparse matrix
#'
#' @export
raw_counts <- function(x){
    return(x@counts)
}

#' Get cluster distances to background
#'
#' @param x An SCE object.
#' @param k_init The k_init run to extract the values from.
#'
#' @return A numeric vector containing the distances for each cluster.
#' 
#' @export
distances <- function(x, k_init = NULL){
    k_init <- check_k_init(x, k_init, return_all = FALSE)

    if ("Dist" %in% names(x@kruns[[k_init]])){
        return(x@kruns[[k_init]]$Dist)
    } else {
        return(NULL)
    }
}

#' Get Alpha parameters for clusters
#'
#' @param x An SCE object.
#' @param k_init The k_init run to extract the values from.
#'
#' @return A droplet by cluster matrix containing the alpha parameters of 
#'  of the Dirichlet-multinomial cluster in each column.
#'
#' @export
Alpha <- function(x, k_init = NULL){
    k_init <- check_k_init(x, k_init, return_all = FALSE)
    if ("Alpha" %in% names(x@kruns[[k_init]]$params)){
        return(x@kruns[[k_init]]$params$Alpha)
    } else {
        return(NULL)
    }
}

#' Get percent of reads aligned to given gene(s)
#'
#' Add a data column for each droplet giving the percentage of raw reads/UMIs 
#' that align to genes given in \code{genes}. The column name is specified by 
#' \code{name}.
#'
#' @param x An SCE object.
#' @param genes Genes to calculate percentage of in counts.
#' @param name Column name to place in dropl_info.
#'
#' @return An SCE object.
#' 
#' @importFrom Matrix colSums
#' 
#' @export
#' 
#' @examples
#' # Add MT%
#' mt_genes <- grep(pattern="^mt-", x=rownames(mb_small@gene_data), ignore.case=TRUE, value=TRUE)
#' mb_small <- get_gene_pct(x = mb_small, genes=mt_genes, name="pct.mt")
#' # Add MALAT1
#' genes <- grep(pattern="^malat1$", x=rownames(mb_small@gene_data), ignore.case=TRUE, value=TRUE)
#' mb_small <- get_gene_pct(x = mb_small, genes=genes, name="MALAT1")
get_gene_pct <- function(x, genes, name){
    gi <- intersect(genes, rownames(x@counts))
    if (length(gi) != length(genes)){
        stop("at least one of given genes not found.")
    }
    gene_pct <- 100 * colSums(x@counts[genes,,drop=FALSE]) / colSums(x@counts)
    x@droplet_data[names(gene_pct),name] <- gene_pct
    return(x)
}

#' Convert an SCE object to Seurat
#'
#' Convert an SCE object to a Seurat object. if \code{targets} is true (default), output only droplets that are 
#' called as not debris. If \code{meta} is TRUE, then output meta data from droplet_info to the meta.data 
#' slot in the Seurat object. Additional functions \code{...} to this function are passed onto 
#' \code{\link[Seurat]{CreateSeuratObject}}. Common arguments include \code{min.cells = 3} and 
#' \code{min.features = 200}.
#'
#' @param x An SCE object.
#' @param targets Logical indicating whether to remove droplets called as debris. Default is TRUE.
#' @param meta Logical that indicates whether to place the data from droplet_info into meta.data in the 
#'  resulting Seurat object. Default is TRUE.
#' @param ... Arguments to \code{\link[Seurat]{CreateSeuratObject}}, such as \code{project} for project name.
#'
#' @return A Seurat object
#' @export
#' @examples
#' \donttest{
#' mm_seur <- convert_to_seurat(x = mb_small, 
#'                              targets = FALSE, 
#'                              min.features = 500, 
#'                              min.cells = 3, 
#'                              project = mb_small@name)
#' }
convert_to_seurat <- function(x, targets = TRUE, meta = TRUE, ...){
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Package \"Seurat\" needed for convert_to_seurat. Please install.",
             call. = FALSE)
    } else {

        if (targets) drops <- get_clean_ids(x)
        else drops <- rownames(x@droplet_data)

        if (meta) meta.data <- x@droplet_data[drops,,drop=FALSE]
        else meta.data <- NULL

        seur <- Seurat::CreateSeuratObject(counts = x@counts[,drops,drop=FALSE], 
                                           meta.data = meta.data, ...)
        return(seur)
    }
}

#' Read 10X counts data
#'
#' Read the counts output from into a sparse matrix. 
#' Given a 10X output path, 
#' Read in the the output files from 10X CellRanger. Specifically, \code{path} 
#' should be be a directory containing the raw count output from 10X 
#' CellRnager with  output files \code{barcodes.tsv} 
#' \code{genes.tsv} \code{matrix.mtx}, or 
#' if using CellRanger v2, or \code{barcodes.tsv.gz} \code{features.tsv.gz} 
#' \code{matrix.mtx.gz} if using CellRanger v3. When reading in the barcodes 
#' and gene names, duplicates are removed with the make.unique function. 
#' The separator is specified by \code{sep} by default.
#'
#' @param path The file path prefix to the raw 10X output directory
#' @param clip_end A logical indicating whether to remove the suffix "-N" 
#'  from the barcodes, where N is n integer such as 1. The default it TRUE.
#' @param sep A character indicating the separator to use in the call 
#'  to \code{\link[base]{make.unique}}.
#'
#' @return expr A sparse matrix of counts from 10X cell ranger
#'  with cells in the columns and genes in the rows.
#'
#' @importFrom Matrix readMM
#' @importFrom methods as
#' @importFrom utils read.delim
#' @export
#'
#' @examples
#' \donttest{
#' dir <- "path/to/10x" # should show 3 files described above
#' counts <- read_10x(dir)
#' }
read_10x <- function(path, clip_end=TRUE, sep="."){
    files <- list.files(path, full.names=TRUE)
    files_names <- list.files(path, full.names=FALSE)
    v3 <- "matrix.mtx.gz" %in% files_names

    if (v3){
        mtx_file <- file.path(path, "matrix.mtx.gz")
        genes_file <- file.path(path, "features.tsv.gz")
        barcode_file <- file.path(path, "barcodes.tsv.gz")
    } else {
        mtx_file <- file.path(path, "matrix.mtx")
        genes_file <- file.path(path, "genes.tsv")
        barcode_file <- file.path(path, "barcodes.tsv")
    }

    if (file.exists(mtx_file)) {
        expr <- readMM(mtx_file)
    } else {
        stop(paste0(mtx_file, " not found"))
    }
    expr <- as(expr, "CsparseMatrix")

    if (file.exists(barcode_file)){
        barcode_names <- readLines(barcode_file)
        if (clip_end) barcode_names <- sub("-.*", "", barcode_names)
    } else {
        stop(paste0(barcode_file, " not found"))
    }

    if (file.exists(genes_file)){
        genes <- read.delim(genes_file, 
                            header = FALSE, 
                            stringsAsFactors = FALSE, 
                            sep = "\t")
    } else {
        stop(paste0(genes_file, " not found"))
    }
    
    colnames(expr) <- make.unique(barcode_names, sep=sep)
    rownames(expr) <- make.unique(genes[,2], sep=sep)
    
    return(expr)
}   

