
#' Run DIEM pipeline
#'
#' Run DIEM pipeline on an SCE object. This function takes as input raw data 
#' from a single-nucleus or single-cell experiment, estimates the 
#' parameters of the Dirichlet-multinomial mixture model, and calculates 
#' posterior probabilities that each droplet belongs to a cluster.
#' The input is an SCE object, which can be created using the 
#' function \code{\link{create_SCE}}. This returns the SCE object with 
#' the estimated parameters. To classify the droplets, run 
#' \code{call_targets} on the returned SCE object.
#' 
#' This is a wrapper function that calls the following functions:
#' \enumerate{
#'  \item \code{\link{set_debris_test_set}}
#'  \item \code{\link{filter_genes}}
#'  \item \code{\link{get_pcs}}
#'  \item \code{\link{init}}
#'  \item \code{\link{get_dist}}
#'  \item \code{\link{rm_close}}
#'  \item \code{\link{run_em}}
#' }
#'
#' The diem function takes the sce object \code{sce} and first 
#' partitions the droplets into the debris and sets. Then, 
#' expressed genes are only included in the analysis. PCA is run 
#' on the normalized counts data subset to the variable genes in 
#' order to reduce the dimensionality for k-means. Initialization 
#' is performed with k-means on the PCs. The distance to the 
#' background distribution is estimated, and clusters that have 
#' a distance below a threshold are removed. Finally, EM is run 
#' to estimate the parameters of the mixture distribution.
#' 
#' @param sce An SCE object.
#' @param min_counts A numeric value that specifies the minimum 
#'  total read/UMI counts a test droplet can have.
#' @param cpm_thresh CPMs must be greater than this value. This is 0 by 
#'  default to include all expressed genes.
#' @param n_var_genes Number of top variable genes to use for PCA.
#' @param lss The span parameter of the loess regression, the parameter 
#'  for the function \code{\link[stats]{loess}}. The loess regression is
#'  used to regress out the effect of mean expression on variance.
#' @param n_pcs Number of PCs to return.
#' @param k_init The number of clusters to initialize k-means. This can be 
#'  1 or several values so that multiple values can be compared. This value 
#'  will not include the debris cluster.
#' @param iter.max_init The maximum number of k-means interations 
#'  for the initialization.
#' @param nstart_init The number of starts to use in k-means for 
#'  for the initialization.
#' @param min_size_init The minimum number of droplets that must belong 
#'  to an initialized cluster.
#' @param seedn The seed for random k-means initialization. 
#'  It is set to 1 by default. If you desire truly random initializations
#'  across runs, set to NULL or different values for each run.
#' @param fltr The filter threshold between 0 and 1 
#'  that controls the minimum distance to 
#'  the background distribution that a cluster can have. Remove those 
#'  centers with a distance less than this value.
#' @param eps The delta threshold for when to call convergence for 
#'  the EM estimation of the Dirichlet-multinomial mixture model. The EM 
#'  stops when delta falls below this value. We define delta as the 
#'  average change in posterior probability. By default this is set to 
#'  1e4, so that the EM converges when less than 1 in 10,000 labels 
#'  change on average.
#' @param max_iter_dm Maximum number of iterations for the EM estimation 
#'  of the Dirichlet-multinomial mixture model.
#' @param min_genes A numeric value that specifies the minimum 
#'  total genes detected a test droplet can have.
#' @param top_n A numeric value giving the maximum number of test droplets. 
#'  This is turned off by default.
#' @param debris_ids A character vector of droplet IDs that will be assigned 
#'  to the debris set, regardless of its total counts or genes detected.
#' @param threads Number of threads for parallel execution. Default is 1.
#' @param verbose Logical indicating verbosity.
#' 
#' @return SCE object.
#' @export
#' @examples
#' \donttest{
#' counts <- read_10x("mouse_nuclei_2k/raw_gene_bc_matrices/mm10/")
#' mb_sce <- create_SCE(x=counts, name="Mouse Brain")
#' mb_sce <- diem(mb_sce)
#' mb_sce <- call_targets(mb_sce, pp_thresh = 0.95, min_genes = 200)
#' mm_seur <- convert_to_seurat(x=mb_sce, min.features = 200, min.cells = 3, project=mb_sce@name)
#' }
diem <- function(sce,
                 min_counts = 100, 
                 cpm_thresh = 0, 
                 n_var_genes = 2000, 
                 lss=0.3, 
                 n_pcs = 30, 
                 k_init = 30, 
                 iter.max_init = 15, 
                 nstart_init = 30, 
                 min_size_init = 10, 
                 fltr = 0.1, 
                 seedn = 1, 
                 eps = 1e-4, 
                 max_iter_dm = 100, 
                 min_genes = 0, 
                 top_n = NULL, 
                 debris_ids = NULL, 
                 threads = 1, 
                 verbose = TRUE){
    sce <- set_debris_test_set(sce, 
                               min_counts = min_counts, 
                               min_genes = min_genes, 
                               top_n = top_n, 
                               debris_ids = debris_ids, 
                               verbose = verbose)
    sce <- filter_genes(sce, 
                        cpm_thresh = cpm_thresh, 
                        verbose = verbose)
    sce <- get_pcs(sce, 
                   n_var_genes = n_var_genes, 
                   lss = lss, 
                   n_pcs = n_pcs)
    sce <- init(sce, 
                k_init = k_init, 
                iter.max_init = iter.max_init, 
                nstart_init = nstart_init, 
                min_size_init = min_size_init, 
                seedn = seedn, 
                threads = threads, 
                verbose = verbose)
    sce <- get_dist(sce, 
                    k_init = k_init, 
                    verbose = verbose)
    sce <- rm_close(sce, 
                    fltr = fltr, 
                    k_init = k_init, 
                    verbose = verbose)
    sce <- run_em(sce, 
                  eps = eps, 
                  fltr = fltr, 
                  max_iter_dm = max_iter_dm, 
                  k_init = k_init, 
                  threads = threads, 
                  verbose = verbose)

    return(sce)
}

