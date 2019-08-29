
#' Run DIEM pipeline
#'
#' Run DIEM pipeline on an SCE object, (see \code{\link{read_10x}} and 
#' \code{\link{create_SCE}} functions).
#'
#' @param sce SCE. SCE object with raw counts.
#' @param top_n Numeric value specifying the  droplets below \code{top_n} (ranked by total counts) 
#'  are fixed to debris.
#' @param min_counts A numeric value that specifies that all droplets below this count threshold are 
#'  always fixed to debris.
#' @param top_thresh A numeric value specifying any droplets above this count threshold are 
#' fixed to nuclei. No droplets are fixed to nuclei by default.
#' @param cpm_thresh The minimum CPM expression threshold in both groups.
#' @param eps Numeric threshold. The EM algorithm converges when the percent change in log likihood is 
#'  less than \code{eps}.
#' @param max_iter The maximum number of iterations allowed to run.
#' @param pp_thresh Numeric threshold, where clean droplets must have a 
#'  posterior probability of at least \code{pp_thresh}.
#' @param min_genes Numeric threshold, where clean droplets must have at least 
#'  \code{min_genes} genes detected.
#' @param verbose Logical indicating verbosity.
#' 
#' @return SCE object.
#' @export
#' @examples
#' counts <- read_10x("mouse_nuclei_2k/raw_gene_bc_matrices/mm10/")
#' mb_sce <- create_SCE(x=counts, name="Mouse Brain")
#' mb_sce <- diem(mb_sce)
#' mm_seur <- convert_to_seurat(x=mb_sce, min.features = 200, min.cells = 3, project=mb_sce@name)
diem <- function(sce,
                 top_n=1e4, 
                 min_counts=150, 
                 cpm_thresh=10, 
                 scale_factor=1, 
                 logt=TRUE, 
                 n_var_genes=2000, 
                 pseudocount=1e-4, 
                 eps=1e-8, 
                 max_iter=100, 
                 nn=50, 
                 use_var_genes=FALSE, 
                 pp_thresh=0.95, 
                 gene_thresh=200, 
                 verbose=TRUE){
    if ((pp_thresh < 0) | (pp_thresh > 1)){
        stop("pp_thresh must be between 0 and 1")
    }
    sce@pp_thresh <- pp_thresh

    # Add MT%
    mt_genes <- grep(pattern="^mt-", x=rownames(sce@gene_data), ignore.case=TRUE, value=TRUE)
    sce <- get_gene_pct(x=sce, genes=mt_genes, name="pct.mt")

    sce <- set_test_set(sce, top_n=top_n, min_counts=min_counts)
    sce <- filter_genes(sce, cpm_thresh=cpm_thresh)
    sce <- normalize_data(sce, scale_factor=scale_factor, logt=logt)
    sce <- get_pcs(sce)
    # sce <- get_snn(sce, nn=nn, weighted=FALSE)
    # sce <- initialize_clusters(sce)
    sce <- get_kclust(sce)
    sce <- run_em(x=sce, eps=eps, max_iter=max_iter, verbose=verbose)
    sce <- call_targets(sce, pp_thresh=sce@pp_thresh, min_genes=gene_thresh)

    if (!sce@emo$converged){
        cat(paste0("Warning: DIEM did not converge within ", as.character(max_iter), " iterations.\n"))
    }

    if (verbose){
        cat(paste0("Identified ", as.character(length(get_clean_ids(sce))), " clean droplets.\n"))
    }

    return(sce)
}

