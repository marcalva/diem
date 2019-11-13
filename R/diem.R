
#' Run DIEM pipeline
#'
#' Run DIEM pipeline on an SCE object. This function takes as input raw data 
#' from a single-nucleus or single-cell experiment and filters out droplets 
#' that show evidence of contamination from debris/ambient RNA. The input 
#' data is an sce object, returned by the function \code{\link{create_SCE}}. 
#' The function outputs the SCE object with meta data specifying the 
#' filtered and clean droplets.
#'
#' The diem function takes the sce object \code{sce} and first 
#' partitions the droplets into the debris, test, and cluster sets. 
#' \code{min_counts} and \code{min_genes} specify the lower limit 
#' that droplets in the test set must have. In addition, 
#' the parameter \code{top_n} specifies that at most the top \code{top_n} 
#' should be included in the test set. The default is 10,000 
#' in accordance with the typical maximum number of droplets 
#' generated during a single-cell experiment, but should be changed 
#' accordingly.
#'
#' If the read/UMI count of the \code{top_n} count-ranked droplet is greater 
#' than \code{min_counts}, then the \code{min_counts} is 
#' increased to this limit. The \code{fix_debris} parameter allows 
#' additional flexibility in specifying those droplets that are 
#' included in the debris set. For example, it may be useful to 
#' set droplets above a pre-determined MT% threshold as debris, 
#' regardless of the total number of read/UMI counts present.
#' In order to remove lowly expressed genes, the Counts 
#' Per Million mapped reads (CPM) are calculated for the 
#' both the test set and debris set. Then, only genes 
#' that have a CPM of at least \code{cpm_thresh} in both the test set and 
#' debris set are kept.
#' 
#' Instead of randomly initializing the EM, cell types are estimated from 
#' droplets that are expected to contain cells/nuclei. The top 
#' \code{cluster_n} droplets are ranked by "gene" or "count", given by the 
#' parameter \code{order_by}. Then, droplets with at least those ranked 
#' counts/genes are included in the cluster set. The data is then normalized 
#' by first calculating the variable genes. A loess regression line 
#' is fit between the log counts and log variance, and the only top genes 
#' ranked by residual are used to initialize the clusters. The number of 
#' genes is specified with \code{n_var}. Optionally, one can use all genes 
#' by setting \code{use_var} to FALSE. The span of the loess regression line 
#' is given by \code{lss} (default is 0.3).
#' The data is normalized by dividing counts by the total counts per droplet. 
#' Then, the counts are multiplied by a scaling factor, given by 
#' \code{sf} (the median of total counts by default). Finally, the data is 
#' log transformed after adding a constant value of 1.
#' After normalization, the k-nearest neighbors are identified in the 
#' cluster set. The number of nearest neighbors is specified by 
#' \code{nn}. Clusters are identified from the KNN graph 
#' using the Louvain algorithm. Finally, only clusters with at least 
#' \code{min_size} (20 by default) droplets are considered cell types.
#'
#' The number of clusters and their initialized multinomial means 
#' are taken from the initial clustering assignments calculated with 
#' \code{\link{initialize_clusters}}. The EM algorithm is run  
#' by repeatedly updating the membership probabilities and then 
#' estimating the MLE parameters of the multinomial mixture model.
#' The algorithm converges when the percent change in 
#' the log likihood is less than \code{eps}. If the 
#' algorithm doesn't converge by \code{max_iter}, it breaks off.
#' The posterior probability is the calculated by taking the 
#' sum of the likelihood fractions across the cell types.
#'
#' Given the posterior probabilities, the droplets with those 
#' less than \code{pp_thresh} and with less than 
#' \code{n_genes} genes detected are marked as debris.
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
#' @param use_var A logical indicating whether to subset the data to include
#'  only variable genes when initializing the clusters.
#'  The default is TRUE as it may better identify cell types.
#' @param n_var Number of variable genes to use when initializing the clusters.
#' @param lss Numeric value of the span parameter of the loess regression 
#'  when identifying variable genes.
#' @param sf The scaling factor for the normalization. 
#'  Either a numeric scaling factor to multiply counts after 
#'  division by column sums, or "median" indicating to multiply by the median 
#'  number of total read/UMI counts in droplets (default).
#' @param nn Number of nearest neighbors to calculate in constructing the 
#'  kNN graph for cluster initialization.
#' @param min_size Numeric value giving the minimum number of droplets in 
#'  cluster for it to be used for initialization as a cell type for EM. 
#' @param eps Numeric threshold. The EM algorithm converges when the 
#'  percent change in log likihood is less than \code{eps}.
#' @param max_iter The maximum number of iterations allowed to run.
#' @param psc Pseudocount to add to multinomial parameters 
#'  to avoid collapsing likelihood to 0.
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
                 top_n=1e4, 
                 min_counts=100, 
                 min_genes=100, 
                 fix_debris=NULL, 
                 cpm_thresh=10, 
                 cluster_n=1000, 
                 order_by="gene", 
                 use_var=TRUE,
                 n_var=2000, 
                 lss=0.3, 
                 sf = "median", 
                 nn=30, 
                 min_size=20, 
                 eps=1e-8, 
                 max_iter=100, 
                 psc=1e-4, 
                 pp_thresh=0.95, 
                 gene_thresh=200, 
                 verbose=TRUE){
    if ((pp_thresh < 0) | (pp_thresh > 1)){
        stop("pp_thresh must be between 0 and 1")
    }
    sce@pp_thresh <- pp_thresh

    sce <- set_debris_test_set(sce, 
                               top_n = top_n, 
                               min_counts = min_counts, 
                               min_genes = min_genes, 
                               fix_debris = fix_debris)
    sce <- filter_genes(sce, 
                        cpm_thresh = cpm_thresh, 
                        verbose = verbose)
    sce <- set_cluster_set(sce, cluster_n = cluster_n, 
                           order_by = order_by, 
                           verbose = verbose)
    sce <- initialize_clusters(sce, 
                               use_var = use_var, 
                               n_var = n_var, 
                               lss = lss,
                               sf = sf, 
                               nn = nn, 
                               min_size = min_size, 
                               verbose = verbose)
    sce <- run_em(x = sce, 
                  eps = eps, 
                  max_iter = max_iter, 
                  psc = psc, 
                  verbose = verbose)
    sce <- call_targets(sce, 
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

    return(sce)
}

