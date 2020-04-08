
#' Output summary stats from a DIEM run
#'
#' Return summary statistics from a DIEM run to evaluate clustering 
#' and classification. The function calculates helpful statistics to 
#' evaluate each cluster, including the debris and putatitive cell 
#' types. This returns a data frame with the following 
#' summary stats per cluster:  \describe{
#' \item{n_droplets}{Number of droplets.}
#' \item{avg_n_genes}{Average number of genes detected.}
#' \item{avg_n_counts}{Average number of counts.}
#' \item{avg_logFC}{Average log fold change of the top \code{top_n} 
#'  marker genes.}
#' \item{genes}{Set of top \code{top_n} marker genes separated by 
#'  semicolons.}
#' }
#' The marker genes are calculated as the genes with the largest 
#' difference in proportion between the 
#'
#' @param x An SCE object.
#' @param top_n Number of top DE genes to return in column.
#' @param k_init Run EM on the \code{k_init} initializations. 
#'  If NULL (default), run on all \code{k_init} initializations.
#'
#' @return A data frame
#'
#' @export
summarize_clusters <- function(x, top_n = 20, k_init = NULL){
    k_init <- check_k_init(x, k_init)
    testdat <- droplet_data(x, type = "test")
    dropdat <- droplet_data(x, type = "all")
    dbr_clust <- x@debris_clusters


    if ( ! "Cluster" %in% colnames(testdat))
        stop("Run DIEM and ", sQuote("assign_clusters"), " before summarize")

    if ( ! "score.debris" %in% colnames(testdat))
        x <- estimate_dbr_score(x)

    dropdat[,"Cluster"] <- "1"
    dropdat[rownames(testdat),"Cluster"] <- testdat[,"Cluster"]
    clusters_all <- sort(as.numeric(unique(dropdat[,"Cluster"])))
    clusters_all <- as.character(clusters_all)
    cells <- setdiff(clusters_all, "1")

    means <- Alpha(x)
    means <- sweep(means, 2, colSums(means), "/")


    # Measures
    cs <- rep("Clean", ncol(means))
    cs[clusters_all %in% dbr_clust] <- "Debris"
    ndrops <- sapply(clusters_all, function(ct){sum(testdat[,"Cluster"] == ct)})
    ngenes <- tapply(testdat[,"n_genes"], 
                     testdat[,"Cluster"], 
                     mean)
    ncounts <- tapply(testdat[,"total_counts"], 
                      testdat[,"Cluster"], 
                      mean)
    ds <- tapply(testdat[,"score.debris"], 
                 testdat[,"Cluster"],
                 mean)

    ret <- data.frame("Type" = cs, 
                      "Cluster" = clusters_all, 
                      "n_droplets" = ndrops[clusters_all], 
                      "avg_n_genes" = ngenes[clusters_all],
                      "avg_n_counts" = ncounts[clusters_all], 
                      "avg_dbr_score" = ds[clusters_all]) 

    # Get DE top genes
    de <- top_genes(x, top_n = top_n, k_init = k_init)

    avgs <- sapply(clusters_all, function(ct){
                   des <- de[de[,"cluster"] == ct,,drop=FALSE]
                   genes <- paste(des[,"gene"], collapse=";")
                   return(c("avg_logFC" = mean(des[,"logFC"]), 
                          "genes" = genes))
                      })
    avgs <- t(avgs)
    ret <- cbind(ret, avgs[clusters_all,])

    return(ret)
}

#' Run a DE between two groups with Welch's t-test
#'
#' Run DE between two groups. \code{x} is a sparse matrix with 
#' normalized expression data (typically proportions). \code{c1} 
#' and \code{c2} are vectors containing indices or names 
#' of the columns that belong to group 1 and 2, respectively. 
#' A Welch's t-test is performed.
#'
#' @param x A sparse matrix
#' @param c1 Column indexes or names of first group
#' @param c2 Column indexes or names of second group
#' 
#' @return A data frame of DE results with
#' \describe{
#' \item{diff}{difference between mu1 and mu2}
#' \item{logFC}{log2 fold change}
#' \item{mu1}{mean of group 1}
#' \item{mu2}{mean of group 2}
#' \item{t}{t statistic}
#' \item{p}{p-value}
#' \item{gene}{gene name}
#' }
#' @importFrom Matrix rowMeans
#' @export
de_ttest <- function(x, c1, c2){
    if (length(c1) != ncol(x))
        stop("c1 must have length = number of columns in x")
    if (length(c2) != ncol(x))
        stop("c2 must have length = number of columns in x")

    mu1 <- rowMeans(x[,c1,drop=FALSE])
    mu2 <- rowMeans(x[,c2,drop=FALSE])

    keep <- mu1 > 0 & mu2 > 0
    x <- x[keep,,drop=FALSE]
    mu1 <- mu1[keep]
    mu2 <- mu2[keep]
    
    var1 <- fast_varCPP(x = x[,c1,drop=FALSE], 
                        mu = mu1)
    var2 <- fast_varCPP(x = x[,c2,drop=FALSE], 
                        mu = mu1)
    n1 <- ncol(x[,c1,drop=FALSE])
    n2 <- ncol(x[,c2,drop=FALSE])

    t <- sapply(1:nrow(x), function(i) {
        r <- (mu1[i] - mu2[i]) / sqrt(var1[i]/n1 + var2[i]/n2)
        return(r)
    })
    p <- sapply(1:length(t), function(i){
                num <- (var1[i]/n1 + var2[i]/n2)^2
                d1 <- (var1[i]/n1)^2 / (n1-1)
                d2 <- (var2[i]/n2)^2 / (n2-1)
                dof <- num / (d1 + d2)
                2 * pt(q = -abs(t[i]), df = dof)
    })
    datf <- data.frame(diff = mu1 - mu2, 
                       logFC = log2(mu1) - log2(mu2), 
                       mu1 = mu1, 
                       mu2 = mu2, 
                       t = t, 
                       p = p, 
                       gene = rownames(x))
    return(datf)
}

#' Run a basic DE between clusters
#'
#' @param means A gene by cluster sata frame of cluster means, 
#'  of which the columns sum to 1.
#' @param weights weight to give each column of means when calculating the 
#'  average proportions in means across columns.
#' @param top_n Number of top DE genes to return.
#' 
#' @return A data frame of DE results
de_basic <- function(means, weights, top_n = 20){
    clusters <- 1:ncol(means)
    if (length(weights) != ncol(means))
        stop("length of weights must be the same as number of columns in means.")
    weights <- weights / sum(weights)
    names(weights) <- clusters
    de <- lapply(clusters, function(cl){
                 cln <- setdiff(clusters, cl)
                 p1 <- means[,cl]
                 if (!is.null(weights)){
                     wn <- weights[cln] / sum(weights[cln])
                     p2 <- means[,cln,drop=FALSE] %*% wn
                 }
                 else {
                     p2 <- rowMeans(means[,cln,drop=FALSE])
                 }
                 di <- p1 - p2
                 logfc <- log2(p1) - log2(p2)
                 datf <- data.frame("diff" = di, 
                                    "logFC" = logfc,
                                    "p1" = p1, 
                                    "p2" = p2, 
                                    "cluster" = cl, 
                                    "gene" = rownames(means))
                 o <- order(di, decreasing = TRUE)
                 datf <- datf[o,]
                 return(datf[1:top_n,])
    })
    ret <- do.call(rbind, de)
    ret$gene <- as.character(ret$gene)
    rownames(ret) <- NULL
    return(ret)
}

#' Get top up-regulated genes per cluster
#'
#' Extract marker genes for each cluster by ranking them based on the 
#' difference of their proportion means from the Dirichlet 
#' distributions estimated by DIEM. Output the top \code{top_n} 
#' genes for each cluster.
#'
#' @param x An SCE object.
#' @param top_n Number of top DE genes to return.
#' @param k_init Run EM on the \code{k_init} initialization(s). 
#'  If NULL (default), run on all \code{k_init} initializations.
#'
#' @return An SCE object.
#'
#' @export
top_genes <- function(x, top_n = 20, k_init = NULL){
    genes.use <- rownames(x@gene_data)[x@gene_data$exprsd]
    k_init <- check_k_init(x, k_init)
    means <- Alpha(x, k_init)
    means <- means[genes.use,]
    means <- sweep(means, 2, colSums(means), "/")

    clusters <- 1:ncol(means)
    testdat <- droplet_data(x, type = "test")
    dropdat <- droplet_data(x, type = "all")
    dropdat[,"Cluster"] <- "1"
    dropdat[rownames(testdat),"Cluster"] <- testdat[,"Cluster"]
    clusters_all <- as.character(clusters)

    w <- sapply(clusters_all, function(ct){
                sum(dropdat[dropdat$Cluster == ct, "total_counts"])
    })
    de <- de_basic(means, top_n = top_n, weights = w)
    return(de)
}

#' Estimate debris score per droplet
#'
#' Create a debris score based on expression of debris-enriched genes.
#' For each droplet, it sums the normalized read counts for these genes. 
#' Debris-enriched genes are calculated using a t-test between 
#' droplets in the debris and non-debris clusters. Genes with an 
#' adjusted p-value less than \code{thresh_p} ser 
#' set as debris-enriched genes. The function starts with cluster 1 
#' in the debris group, calculates a debris cores, and then 
#' assigns clusters to the debris group if their average is greater than 
#' or equal to \code(thresh_score}. This is repeated until the set of 
#' debris clusters remains the same as in the previous iteration.
#' 
#'
#' @param x An SCE object.
#' @param thresh_score Threshold for cluster mean debris score. Clusters 
#'  with a mean above this value are considered debris clusters.
#' @param thresh_genes Threshold for cluster mean number of genes 
#'  detected. Clusters with a mean number of genes less than this value 
#'  are set to debris clusters.
#' @param thresh_p P-value threshold for calculating differential 
#'  expression between droplets in the debris and non-debris clusters.
#' @param p_method Method for p-value correction using 
#'  \code{\link[stats]{p.adjust}}.
#'  expression between droplets in the debris and non-debris clusters.
#' @param dbr_fixed Specifies additional clusters to set as debris. Cluster 
#'  1 is always included. Use this only if you know what you are doing.
#' @param verbose Verbosity, indicating whether to show a progress
#'  bar.
#' @param k_init Run EM on the \code{k_init} initialization(s). 
#'  If NULL (default), run on all \code{k_init} initializations.
#' 
#' @importFrom Matrix colSums
#'
#' @return An SCE object
#' @export
estimate_dbr_score <- function(x, 
                               thresh_score = 0.5, 
                               thresh_genes = 200, 
                               thresh_p = 0.05, 
                               p_method = "fdr", 
                               dbr_fixed = NULL, 
                               verbose = TRUE, 
                               k_init = NULL){
    name <- "score.debris"

    k_init <- check_k_init(x, k_init)

    gene_mean <- tapply(x@test_data[,"n_genes"], x@test_data[,"Cluster"], mean)
    lc_clust <- which( gene_mean < thresh_genes)

    fixed <- unique(c(1, dbr_fixed, lc_clust))
    dbr_clust <- c()
    new_dbr_clust <- c(fixed, lc_clust)

    # Specify all clusters
    test_cl <-  x@test_data[,"Cluster"]

    # Normalize test counts
    droplets.use <- rownames(x@test_data)
    genes.use <- rownames(x@gene_data)[x@gene_data$exprsd]

    counts_test <- x@counts[genes.use, droplets.use]
    counts_test <- divide_by_colsum(counts_test)
    sf <- min(colSums(counts_test))
    counts_test <- log1p(sf * counts_test)

    while( !identical(dbr_clust, new_dbr_clust) ){
        dbr_clust <- new_dbr_clust

        # Specify debris droplets
        c1 <- test_cl %in% dbr_clust
        c2 <- !c1

        # Get DE
        de <- de_ttest(counts_test, c1, c2)
        de <- de[order(de[,"diff"], decreasing = TRUE),]
        de[,"p_adj"] <- p.adjust(de[,"p"], method = p_method)
        keep <- de[,"logFC"] > 0 & de[,"p_adj"] < thresh_p
        if (sum(keep) == 0)
            stop("No DE genes found: try increasing thresh_p")
        de <- de[keep,,drop=FALSE]

        scores <- colSums(counts_test[de$gene,])
        clust_mean <- tapply(scores, x@test_data[,"Cluster"], mean)
        scores <- scores - min(clust_mean)
        dbr_mean <- mean(scores[x@test_data[,"Cluster"] %in% dbr_clust])
        scores <- scores / dbr_mean

        clust_mean <- tapply(scores, x@test_data[,"Cluster"], mean)
        new_dbr_clust <- which(clust_mean >= thresh_score)
        new_dbr_clust <- unique(c(fixed, new_dbr_clust))
    }
    x@test_data[,name] <- scores
    x@debris_clusters <- dbr_clust

    if (verbose){
        message("calculated debris scores using ", length(dbr_clust), 
                " debris clusters and ", nrow(de), " debris genes")
    }

    return(x)
}

