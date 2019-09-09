
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
                 min_counts=100, 
                 min_genes=100, 
                 fix_debris=NULL, 
                 cpm_thresh=10, 
                 cluster_n=1000, 
                 order_by="gene", 
                 nn=30, 
                 use_var=TRUE,
                 n_var=2000, 
                 lss=0.3, 
                 min_size=20, 
                 eps=1e-8, 
                 max_iter=100, 
                 pseudocount=1e-4, 
                 pp_thresh=0.95, 
                 gene_thresh=200, 
                 verbose=TRUE){
    if ((pp_thresh < 0) | (pp_thresh > 1)){
        stop("pp_thresh must be between 0 and 1")
    }
    sce@pp_thresh <- pp_thresh

    sce <- set_debris_test_set(sce, top_n=top_n, min_counts=min_counts, min_genes=min_genes, fix_debris=fix_debris)
    sce <- filter_genes(sce, cpm_thresh=cpm_thresh)
    sce <- initialize_clusters(sce, 
                               cluster_n=cluster_n, 
                               order_by=order_by, 
                               nn=nn, 
                               use_var=use_var, 
                               n_var=n_var, 
                               lss=lss, 
                               min_size=min_size, 
                               verbose=verbose)
    sce <- run_em(x=sce, eps=eps, max_iter=max_iter, psc=psc, verbose=verbose)
    sce <- call_targets(sce, pp_thresh=sce@pp_thresh, min_genes=gene_thresh)

    if (verbose){
        n_clean <- length(get_clean_ids(sce))
        n_rm <- length(sce@test_set) - n_clean
        cat(paste0("Removed ", as.character(n_rm), " debris droplets from the test set.\n"))
        cat(paste0("Kept ", as.character(n_clean), " clean droplets.\n"))
    }

    return(sce)
}

