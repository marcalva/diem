
#' Return design matrix of labs vector
#'
#' @param labs Vector of labels. Ordered by factor and corresponds to 
#'  columns of returned matrix.
#'
#' @return A droplet by cluster matrix
z_table <- function(labs){
    labs <- factor(labs)
    nl <- nlevels(labs)
    z <- matrix(0, nrow = length(labs), ncol = nl)
    for (i in 1:nl){
        z[unclass(labs) == i,i] <- 1
    }
    rownames(z) <- names(labs)
    colnames(z) <- levels(labs)
    return(z)
}

#' Run PCA
#'
#' Get PCs from normalized counts.
#' 
#' @param x An SCE object.
#'  many genes detected.
#' @param n_pcs Number of PCs to calculate.
#' @param seedn The seed to set for irlba PCA calculation. It is set to 
#'  1 for reproducibility but can be set to NULL for a random 
#'  initialization.
#'
#' @return An SCE object with PCs
#'
#' @importFrom irlba prcomp_irlba
#' @importFrom stats prcomp
#' @importFrom Matrix t
run_pca <- function(x, n_pcs = 30, seedn = 1){
    set.seed(seedn, kind = "Mersenne-Twister")
    m <- nrow(x@norm[x@vg,])
    n <- ncol(x@norm[x@vg,])
    if (n_pcs >= 0.5 * min(m, n)){
        pcs <- prcomp(t(x@norm[x@vg,]), center = TRUE, scale. = TRUE)
    } else {
        pcs <- prcomp_irlba(t(x@norm[x@vg,]), n = n_pcs, center = TRUE, scale. = TRUE)
    }
    pcs <- pcs$x
    rownames(pcs) <- colnames(x@norm[x@vg,])
    x@pcs <- pcs[,1:n_pcs,drop=FALSE]
    return(x)
}

#' Get PCs
#'
#' Run PCA and get top \code{n_pcs} PCs for the test set. 
#' The PCs are used as the features for the initial k-means clustering. 
#' Only droplets with at least \code{min_genes} are used in the PCA, and 
#' thus used in the initialization.
#' The counts data for the test set are count-normalized to the median 
#' and log transformed. Then the top \code{n_var_genes} 
#' variable genes are calculated using the function 
#' \code{\link{get_var_genes}}. PCA is run on the normalized count data 
#' for these variable genes only. 
#'
#' @param x An SCE object.
#' @param min_genes Calculate PCs from droplets with at least this 
#'  many genes detected.
#' @param n_var_genes Number of top variable genes to use for PCA.
#' @param lss The span parameter of the loess regression, the parameter 
#'  for the function \code{\link[stats]{loess}}. The loess regression is
#'  used to regress out the effect of mean expression on variance.
#' @param threads Number of threads for parallel execution. Default is 1.
#' @param n_pcs Number of PCs to return.
#' @param seedn The seed to set for irlba PCA calculation. It is set to 
#'  1 for reproducibility but can be set to NULL for a random 
#'  initialization.
#'
#' @return An SCE object with PCs
#'
#' @export
#' @examples
#' \donttest{
#'
#' # Get PCs with default parameters
#' sce <- get_pcs(sce)
#'
#' # Run initialization with droplets that have at least 150 genes
#' # detected
#' sce <- get_pcs(sce, min_genes = 150)
#'
#' # Using top 3,000 variable genes
#' sce <- get_pcs(sce, n_var_genes = 3000)
#'
#' # Use top 50 PCs for initialization
#' sce <- get_pcs(sce, n_pcs = 50)
#'
#' # Return PCs from random irlba initializations
#' sce <- get_pcs(sce, seedn = NULL)
#' sce <- get_pcs(sce, seedn = NULL)
#' sce <- get_pcs(sce, seedn = NULL)
#' 
#' }
get_pcs <- function(x, 
                    min_genes = 200,
                    n_var_genes = 2000, 
                    lss = 0.3, 
                    threads = 1, 
                    n_pcs = 30, 
                    seedn = 1){
    keep <- x@test_data[,"n_genes"] >= min_genes
    droplets.use <- rownames(x@test_data)[keep]

    x <- get_var_genes(x, 
                       droplets.use = droplets.use, 
                       n_genes = n_var_genes, 
                       lss = lss, 
                       threads = threads)
    x <- normalize_data(x, droplets.use = droplets.use)
    x <- run_pca(x, n_pcs = n_pcs, seedn = seedn)
    return(x)
}

#' Re-assign labels by condensing to 1 to k
#' 
#' @param labs Vector of labels.
#'
#' @return Vector of labels in which the numeric values are contiguous.
condense_labs <- function(labs){
    cts <- setdiff(sort(unique(labs)), "1")
    to <- c("1", as.character(seq(2, length(cts) + 1)))
    from <- c("1", cts)
    names(to) <- from
    nms <- names(labs)
    labs <- to[as.character(labs)]
    names(labs) <- nms
    return(labs)
}

#' Merge small clusters
#' 
#' Merge clusters given in labels until all clusters have a size of at 
#' least \code{min_size}. The \code{dat} matrix contains the data used 
#' to generate the cluster means.
#' Data points are assigned to the next closest cluster.
#'
#' @param labs A vector of labels, with numeric values corresponsing 
#'  to clusters
#' @param dat A feature by observation data matrix used to get cluster means.
#' @param min_size The minimum number of members that belong to a cluster.
#' @param verbose Verbosity.
#'
#' @importFrom stats dist
#' 
#' @return A vector containing new labels
merge_size <- function(labs, dat, min_size = 10, verbose = TRUE){
    labs <- as.character(labs)
    freq <- table(labs)
    nm <- 0
    while(min(freq) < min_size){
        clusts <- names(freq)
        small_clust <- names(freq)[which.min(freq)]

        # Find closest clusters
        means <- sapply(clusts, function(i) colMeans(dat[labs == i,,drop=FALSE]))
        dists <- as.matrix(dist(t(means), upper = TRUE, diag = TRUE))
        others <- colnames(dists)[colnames(dists) != small_clust]
        dists <- dists[others,,drop=FALSE]
        closest <- rownames(dists)[which.min(dists[,small_clust])]
        labs[labs == small_clust] <- closest
        freq <- table(labs)
        nm <- nm + 1
    }
    if (verbose){
        if (nm > 0){
            message("merged ", nm, 
                    " small ", 
                    ngettext(nm, "cluster", "clusters"))
        }
    }

    labs <- condense_labs(labs)
    return(labs)
}

#' Initialize clusters
#' 
#' Initialize the parameters of the mixture model. First 
#' run k-means on the principal components of the test droplets to 
#' initialize the cluster memberships. Then use these memberships 
#' to estimate alpha and pi, the parameters of the mixture model. 
#' The mixture model can be one either a Dirichlet-multinomial or 
#' a multinomial, specified by the \code{model} parameter. The 
#' multinomial is much faster but may be less exact in the 
#' clustering step.
#' The parameter \code{k_init} specifies 
#' initial number of clusters k to set for k-means. These clusters 
#' are input into the mixture model along with the labeled 
#' debris clusters.
#'
#' @param x An SCE object.
#' @param k_init The number of clusters to initialize k-means. This 
#'  does not include the debris cluster as the only the unlabeled 
#'  test droplets are classified here.
#' @param nstart_init The number of starts to use in k-means for 
#'  for the initialization.
#' @param min_size_init The minimum number of droplets that must belong 
#'  to an initialized cluster.
#' @param model The mixture model to assume. Can be either "DM" for 
#'  a Dirichlet-multinomial or "mltn" for a multinomial.
#' @param seedn The seed for random k-means initialization. 
#'  It is set to 1 by default. If you desire truly random initializations
#'  across runs, set to NULL or different values for each run.
#' @param threads Number of threads for parallel execution. Default is 1.
#' @param verbose Verbosity.
#'
#' @return An SCE object
#'
#' @importFrom Matrix Diagonal colSums rowMeans rowSums t
#' @importFrom stats kmeans
#'
#' @export
#' @examples
#' \donttest{
#'
#' # Initialize parameters with default values and multiple threads
#' sce <- init(sce, threads = 8)
#'
#' # Specify initial k of 10 cell type clusters
#' sce <- init(sce, k_init = 10, threads = 8)
#'
#' # Initialize parameters for Dirichlet-multinomial
#' sce <- init(sce, k_init = 10, model = "DM", threads = 8)
#'
#' # Set seedn to NULL for random starts
#' sce <- init(sce, k_init = 10, seedn = NULL, threads = 8)
#' sce <- init(sce, k_init = 10, seedn = NULL, threads = 8)
#' sce <- init(sce, k_init = 10, seedn = NULL, threads = 8)
#'
#' # Specify initial k of 30 cell type clusters and 
#' # Only allow clusters with at least 30 droplets
#' sce <- init(sce, k_init = 30, min_size_init = 30, threads = 8)
#
#' }
init <- function(x, 
                 k_init = 20, 
                 nstart_init = 30, 
                 min_size_init = 10, 
                 model = "mltn", 
                 seedn = 1, 
                 threads = 1, 
                 verbose = TRUE){ 
    genes.use <- rownames(x@gene_data)[x@gene_data$exprsd]

    sizes <- colSums(x@counts[genes.use,])

    if (length(x@pcs) == 0){
        stop("calculate PCs before k-means initialization")
    }

    pcs_s <- scale(x@pcs)

    labs <- rep(0, ncol(x@counts[genes.use,]))
    names(labs) <- colnames(x@counts[genes.use,])
    labs[x@bg_set] <- 1

    all_ret <- list()

    k <- as.integer(k_init)
    x@k_init <- k
    if (verbose){
        message("initializing parameters for k_init = ", k)
    }

    if (verbose) message("running k-means")
    old_opt <- options(warn = -1)
    set.seed(seedn, kind = "Mersenne-Twister")

    km <- kmeans(x = pcs_s, 
                 centers = k, 
                 nstart = nstart_init)
    options(old_opt)
    kclust <- km$cluster # integer vector of 1 to k values, corresponding to column in pcs_s

    # Merge small clusters from k-means
    kclust <- merge_size(kclust, pcs_s, min_size_init, verbose)
    kclust <- as.numeric(kclust)
    names(kclust) <- rownames(pcs_s)
    kclust <- kclust

    labs_k <- labs
    labs_k[names(kclust)] <- kclust + 1

    if (verbose) message("estimating initial parameters")
    Z <- z_table(kclust)
    Z <- cbind(0, Z)

    if (model == "DM"){
        Alpha <- get_alpha_dm(counts = x@counts[genes.use, rownames(Z)], 
                           Z = Z, 
                           test_set = x@test_set, 
                           bg_set = x@bg_set, 
                           threads = threads)
    } else if (model == "mltn") {
        fixed_c <- rowSums(x@counts[genes.use,x@bg_set])
        Alpha <- get_alpha_mult(counts = x@counts[genes.use, rownames(Z)], 
                                Z = Z, 
                                add = fixed_c, 
                                add_to = 1)
    } else {
        stop("model must be either ", sQuote("DM"), " or ", sQuote("mltn"))
    }
    Pi <- get_pi(Z, add = length(x@bg_set), add_to = 1)
    llk <- get_llk(counts = x@counts[genes.use, x@test_set], 
                   Alpha = Alpha, 
                   sizes = sizes[x@test_set], 
                   model = model, 
                   threads = threads)

    params <- list("Alpha" = Alpha, "Pi" = Pi)
    x@model <- list("params" = params, "llk" = llk)
    return(x)
}

