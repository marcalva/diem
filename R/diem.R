
#' Run DIEM pipeline
#'
#' Run DIEM pipeline on an SCE object. This function takes as input raw data 
#' from a single-nucleus or single-cell experiment and filters out droplets 
#' that show evidence of contamination from debris/ambient RNA. The input 
#' data is an sce object, returned by the function \code{\link{create_SCE}}. 
#' The function outputs the SCE object with meta data specifying the 
#' filtered and clean droplets.
#'
#' @param sce An SCE object.
#' @param top_n A numeric value giving the total number of test droplets 
#'  and take the counts of the \code{top_n} ranked droplet. Set to \code{NULL} 
#'  to disable and set all droplets above \code{min_counts} and 
#'  \code{min_genes} to the test set.
#' @param min_counts A numeric value that specifies the minimum 
#'  total read/UMI counts a test droplet can have.
#' @param min_genes A numeric value that specifies the minimum 
#'  total genes detected a test droplet can have.
#' @param fix_debris A character vector of droplet IDs that will be assigned 
#'  to the debris set, regardless of its total counts or genes detected.
#' @param cpm_thresh The minimum CPM threshold for removing genes.
#' @param cluster_n Numeric value specifying the number of droplets to use 
#'  in the test set.
#' @param order_by Whether to order the droplets by total number of total 
#'  counts or total number of genes detected in calculating the 
#'  cluster set.
#' @param lss Numeric value of the span parameter of the loess regression 
#'  when identifying variable genes.
#' @param eps Numeric threshold. The EM algorithm converges when the 
#'  percent change in log likihood is less than \code{eps}.
#' @param max_iter The maximum number of iterations allowed to run.
#' @param pp_thresh Numeric threshold, where clean droplets must have a 
#'  posterior probability of at least \code{pp_thresh}.
#' @param gene_thresh Numeric threshold, where clean droplets must have at least 
#'  \code{min_genes} genes detected.
#' @param verbose Logical indicating verbosity.
#' 
#' @return SCE object.
#' @export
#' @examples
#' \donttest{
#' counts <- read_10x("mouse_nuclei_2k/raw_gene_bc_matrices/mm10/")
#' mb_sce <- create_SCE(x=counts, name="Mouse Brain")
#' mb_sce <- diem(mb_sce)
#' mm_seur <- convert_to_seurat(x=mb_sce, min.features = 200, min.cells = 3, project=mb_sce@name)
#' }
diem <- function(sce,
                 top_n=NULL, 
                 min_counts=100, 
                 min_genes=0, 
                 fix_debris=NULL, 
                 order_by="gene", 
                 cpm_thresh=0, 
                 n_var=2000, 
                 lss=0.3, 
                 n_pcs = 30, 
                 k_init = 30, 
                 iter.max = 15, 
                 nstart = 30, 
                 min_size = 10, 
                 fltr = 0.1, 
                 seedn = 1, 
                 psc = 1e-10, 
                 eps = 1e-4, 
                 max_iter = 100, 
                 pp_thresh = 0.95, 
                 gene_thresh = 200, 
                 verbose = TRUE){
    if ((pp_thresh < 0) | (pp_thresh > 1)){
        stop("pp_thresh must be between 0 and 1")
    }
    sce@pp_thresh <- pp_thresh

    sce <- set_debris_test_set(sce, 
                               top_n = top_n, 
                               min_counts = min_counts, 
                               min_genes = min_genes, 
                               fix_debris = fix_debris, 
                               verbose = verbose)
    sce <- filter_genes(sce, 
                        cpm_thresh = cpm_thresh, 
                        verbose = verbose)
    sce <- init(sce, 
                n_pcs = n_pcs, 
                k_init = k_init, 
                iter.max = iter.max, 
                fltr = fltr, 
                nstart = nstart, 
                min_size = min_size, 
                seedn = seedn, 
                psc = psc, 
                verbose = verbose)
    sce <- run_em(sce, 
                  k_init = k_init, 
                  eps = eps, 
                  fltr = fltr, 
                  max_iter = max_iter, 
                  verbose = verbose)
    sce <- call_targets(sce, 
                        k_init = k_init, 
                        pp_thresh = sce@pp_thresh, 
                        min_genes = gene_thresh)

    if (verbose){
        n_clean <- length(get_clean_ids(sce))
        n_rm <- length(sce@test_set) - n_clean
        message(paste0("Removed ", 
                       as.character(n_rm), 
                       " debris droplets from the test set."))
        message(paste0("Kept ", 
                       as.character(n_clean), 
                       " clean droplets."))
    }

    # sce <- remove_debris(sce, verbose = verbose)

    return(sce)
}

