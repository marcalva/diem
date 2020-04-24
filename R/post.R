
#' Output summary stats from a DIEM run
#'
#' Return summary statistics from a DIEM run to evaluate clustering 
#' and classification. The function calculates helpful statistics to 
#' evaluate each cluster, including the debris and putatitive cell 
#' types. This returns a data frame with the following 
#' summary stats per cluster:
#' \describe{
#' \item{Type}{Whether the cluster was considered as debris for 
#'  calculating the debris score.}
#' \item{Cluster}{Cluster for which summary stats are shown.}
#' \item{n_droplets}{Number of droplets.}
#' \item{avg_n_genes}{Average number of genes detected.}
#' \item{avg_n_counts}{Average number of counts.}
#' \item{avg_dbr_score}{Average debris score, if available.}
#' \item{genes}{Set of top \code{top_n} marker genes separated by 
#'  semicolons.}
#' }
#' The marker genes are calculated as the genes with the largest 
#' difference in proportion between the 
#'
#' @param x An SCE object.
#' @param top_n Number of top DE genes to return in column.
#'
#' @return A data frame
#'
#' @export
summarize_clusters <- function(x, top_n = 20){
    testdat <- droplet_data(x, type = "test")
    dropdat <- droplet_data(x, type = "all")
    dbr_clust <- x@debris_clusters


    if ( ! "Cluster" %in% colnames(testdat))
        stop("Run DIEM and ", sQuote("assign_clusters"), " before summarize")

    dropdat[,"Cluster"] <- "1"
    dropdat[rownames(testdat),"Cluster"] <- testdat[,"Cluster"]
    clusters_all <- sort(as.numeric(unique(dropdat[,"Cluster"])))
    clusters_all <- as.character(clusters_all)
    cells <- setdiff(clusters_all, "1")

    means <- Alpha(x)
    if (is.null(means))
        stop("Estimate parameters before ", sQuote("summarize_clusters"))

    means <- sweep(means, 2, colSums(means), "/")


    # Measures
    cs <- rep("Clean", ncol(means))
    cs[1] <- "Debris"
    ndrops <- sapply(clusters_all, function(ct){sum(testdat[,"Cluster"] == ct)})
    ngenes <- tapply(testdat[,"n_genes"], 
                     testdat[,"Cluster"], 
                     mean)
    ncounts <- tapply(testdat[,"total_counts"], 
                      testdat[,"Cluster"], 
                      mean)
    if ("score.debris" %in% colnames(testdat)){
        cs[clusters_all %in% dbr_clust] <- "Debris"
        ds <- tapply(testdat[,"score.debris"], 
                     testdat[,"Cluster"],
                     mean)
    } else {
        ds <- rep(NA, length(clusters_all))
        names(ds) <- clusters_all
    }

    ret <- data.frame("Type" = cs, 
                      "Cluster" = clusters_all, 
                      "n_droplets" = ndrops[clusters_all], 
                      "avg_n_genes" = ngenes[clusters_all],
                      "avg_n_counts" = ncounts[clusters_all], 
                      "avg_dbr_score" = ds[clusters_all]) 

    # Get DE top geneu
    de <- top_genes(x, top_n = top_n)

    gene_ids <- sapply(clusters_all, function(ct){
                       des <- de[de[,"Cluster"] == ct,,drop=FALSE]
                       genes <- paste(des[,"gene"], collapse=";")
                       return(genes)
                      })

    ret[,"genes"] <- gene_ids

    return(ret)
}

#' Run a DE between two groups with Welch's t-test
#'
#' Run DE between two groups. \code{x} is a sparse matrix with 
#' count data. These counts can be normalized, or raw counts 
#' can be given. If providiing raw counts, set the \code{normalize} 
#' parameter to TRUE. \code{c1} 
#' and \code{c2} are vectors containing indices or names 
#' of the columns, or boolean indicies for the droplets 
#" that belong to group 1 and 2, respectively. 
#' A Welch's t-test is performed.
#'
#' @param x A sparse matrix
#' @param c1 Column indexes or names of first group
#' @param c2 Column indexes or names of second group
#' @param normalize Whether to scale the droplets to sum to \code{sf} and 
#'  and then log transform.
#' @param sf Scale the counts in each droplet to sum to this value.
#' 
#' @return A data frame of DE results with
#' \describe{
#' \item{gene}{gene name}
#' \item{diff}{difference between mu1 and mu2}
#' \item{logFC}{log2 fold change}
#' \item{mu1}{mean of group 1}
#' \item{mu2}{mean of group 2}
#' \item{t}{t statistic}
#' \item{p}{p-value}
#' }
#'
#' @importFrom Matrix rowMeans
#' @importFrom stats pt
#'
#' @export
de_ttest <- function(x, c1, c2, normalize = FALSE, sf = 1){

    if (normalize){
        x <- divide_by_colsum(x)
        x <- log1p(sf * x)
    }

    mu1 <- rowMeans(x[,c1,drop=FALSE])
    mu2 <- rowMeans(x[,c2,drop=FALSE])

    keep <- mu1 > 0 & mu2 > 0
    x <- x[keep,,drop=FALSE]
    mu1 <- mu1[keep]
    mu2 <- mu2[keep]
    
    var1 <- fast_varCPP(x = x[,c1,drop=FALSE], 
                        mu = mu1)
    var2 <- fast_varCPP(x = x[,c2,drop=FALSE], 
                        mu = mu2)
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
    datf <- data.frame(gene = rownames(x), 
                       diff = mu1 - mu2, 
                       logFC = log2(mu1) - log2(mu2), 
                       mu1 = mu1, 
                       mu2 = mu2, 
                       t = t, 
                       p = p) 
    datf[,"gene"] <- as.character(datf[,"gene"])
    rownames(datf) <- rownames(x)
    return(datf)
}

#' Run a DE across clusters with Welch's t-test
#'
#' Given expression data in \code{counts}, run differential 
#' expression between droplets in one group against all others, for 
#' each group in \code{labels}.
#' These counts can be normalized, or raw counts 
#' can be given. If providiing raw counts, set the \code{normalize} 
#' parameter to TRUE. Returns only genes with an adjust p-value 
#' less than \code{p_thresh} and a log fold change greater than 
#' \code{logFC}.
#'
#' @param counts Gene by droplet sparse matrix of expression data.
#' @param labels Vector of droplet labels giving the cluster that the 
#'  droplet belongs to.
#' @param normalize Whether to scale the droplets to sum to \code{sf} and 
#'  and then log transform.
#' @param sf Scale the counts in each droplet to sum to this value.
#' @param p_thresh Return only genes that have an adjusted p-value less than 
#'  this value.
#' @param logfc_thresh Return genes that have a log fold change of at least 
#'  this value.
#' @param p_method Adjust p-values for multiple testing using this method.
#' 
#' @return A data frame of DE results with
#' \describe{
#' \item{gene}{gene name}
#' \item{diff}{difference between mu1 and mu2}
#' \item{logFC}{log2 fold change}
#' \item{mu1}{mean of group 1}
#' \item{mu2}{mean of group 2}
#' \item{t}{t statistic}
#' \item{p}{p-value}
#' \item{p_adj}{p-value corrected for multiple testing}
#' \item{cluster}{cluster the DE genes belong to}
#' }
#'
#' @importFrom Matrix rowMeans
#' @importFrom stats pt
#'
#' @export
#'
#' @examples
#' \donttest{
#' counts <- raw_counts(sce, type = "test")
#' clusts <- clusters(sce)
#' markers <- de_ttest_all(counts, clusts)
#' }
de_ttest_all <- function(counts, 
                         labels, 
                         normalize = TRUE,
                         sf = 1, 
                         p_thresh = 0.05, 
                         logfc_thresh = 0, 
                         p_method = "bonferroni"){ 

    if (normalize){
        counts_p <- divide_by_colsum(counts)
        counts_s <- log1p(sf * counts_p)
    } else {
        counts_s <- counts
    }

    clusters <- sort(unique(labels))

    de_all <- list()
    for (c1 in clusters){
        c2 <- setdiff(clusters, c1)
        l1 <- labels %in% c1
        l2 <- labels %in% c2

        datf <- de_ttest(counts_s, l1, l2)
        datf$cluster <- c1
        datf$p_adj <- p.adjust(datf$p, method = p_method)
        datf <- datf[order(datf$t, decreasing = TRUE),]
        datf <- datf[(datf$p_adj < p_thresh) & (datf$logFC > logfc_thresh),,drop=FALSE]
        de_all[[as.character(c1)]] <- datf
    }
    de <- do.call(rbind, de_all)
    de[,"gene"] <- as.character(de[,"gene"])
    rownames(de) <- NULL
    return(de)
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
                                    "Cluster" = cl, 
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
#' difference of their proportion means from the multinomial
#' distributions estimated by DIEM. Output the top \code{top_n} 
#' genes for each cluster.
#'
#' @param x An SCE object.
#' @param top_n Number of top DE genes to return.
#'
#' @return An SCE object.
#'
#' @export
top_genes <- function(x, top_n = 20){
    genes.use <- rownames(x@gene_data)[x@gene_data$exprsd]
    means <- Alpha(x)
    if (is.null(means))
        stop("Estimate parameters before getting top DE genes")

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
#' droplets in the debris and non-debris clusters. Only droplets in the 
#' test set are used. Genes with an 
#' adjusted p-value less than \code{thresh_p}, and a log fold change 
#' greater than \code{thresh_logfc} are included in the debris 
#' set. The number of genes included in the debris set can be capped 
#' by setting the \code{max_genes} parameter. This may help ensure 
#' that only the most specific genes are used in the case the 
#' many genes are significant in the DE test.
#' 
#' Differential expression is run by scaling the droplet counts to 
#' sum to \code{sf} and log transforming. Then, a t-test is run 
#' between droplets in the debris set and the cluster being tested.
#'
#' @param x An SCE object.
#' @param thresh_genes Threshold for cluster mean number of genes 
#'  detected. Clusters with an average number of genes detected 
#'  less than this value are set to debris clusters.
#' @param thresh_p P-value threshold for calculating differential 
#'  expression between droplets in the debris and non-debris clusters.
#' @param thresh_logfc Debris genes must have a log2 fold change greater 
#'  than this value.
#' @param max_genes The maximum number of debris genes to include, ranked 
#'  by difference in proportion. This helps if a large number of genes 
#'  are found significantly DE, and may improve the specificity of the 
#'  debris score.
#' @param p_method Method for p-value correction using 
#'  \code{\link[stats]{p.adjust}}.
#'  expression between droplets in the debris and non-debris clusters.
#' @param sf Scaling factor to multiple normalized counts by before 
#'  log transformation. The normalized counts are used for the 
#'  DE analysis.
#' @param verbose Verbosity, indicating whether to show a progress
#'  bar.
#' 
#' @importFrom Matrix colSums
#' @importFrom stats p.adjust
#'
#' @return An SCE object
#' @export
estimate_dbr_score <- function(x, 
                               thresh_genes = 200, 
                               thresh_p = 0.05, 
                               thresh_logfc = 0, 
                               max_genes = 5000, 
                               p_method = "fdr", 
                               sf = 1,
                               verbose = TRUE){ 
    name <- "score.debris"

    if (thresh_logfc < 0)
        stop(sQuote("thresh_logfc"), " must be at least 0, not ", thresh_logfc)

    gene_mean <- tapply(x@test_data[,"n_genes"], x@test_data[,"Cluster"], mean)
    lc_clust <- which(gene_mean < thresh_genes)
    dbr_clust <- unique(c(1, lc_clust))

    test_dat <- x@test_data
    test_drop <- rownames(test_dat)
    test_genes <- rownames(x@gene_data)[x@gene_data$exprsd]
    test_counts <- x@counts[test_genes, test_drop, drop = FALSE]
    test_clusters <- test_dat[,"Cluster"]

    c1 <- test_clusters %in% dbr_clust
    c2 <- ! c1

    if (sum(c1) == 0)
        stop("No test droplets in the debris cluster(s)")

    counts_s <- divide_by_colsum(test_counts)
    counts_s <- log1p(sf * counts_s)

    de <- de_ttest(counts_s, c1, c2)

    # Specify all clusters
    de <- de[order(de[,"diff"], decreasing = TRUE), , drop = FALSE]
    de[,"p_adj"] <- p.adjust(de[,"p"], method = p_method)
    keep <- (de[,"p_adj"] < thresh_p) & (de[,"logFC"] > thresh_logfc)
    
    de_genes <- de[keep,"gene"]
    nr <- min(length(de_genes), max_genes)
    de_genes <- de_genes[1:nr]
    de[,"Debris"] <- FALSE
    de[de_genes,"Debris"] <- TRUE

    if (sum(keep) == 0){
        x@debris_genes <- de
        x@debris_clusters <- dbr_clust
        message("No DE genes found: cannot calculate debris scores")
        return(x)
    }

    scores <- colSums(counts_s[de_genes, test_drop, drop = FALSE])
    clust_mean <- tapply(scores, x@test_data[,"Cluster"], mean)

    scores <- scores - min(clust_mean)
    dbr_mean <- mean(scores[x@test_data[,"Cluster"] %in% dbr_clust])
    scores <- scores / dbr_mean

    x@test_data[,name] <- as.numeric(scores)
    x@debris_genes <- de
    x@debris_clusters <- dbr_clust

    if (verbose){
        message("calculated debris scores using ", length(dbr_clust), 
                " debris clusters and ", length(de_genes), " debris genes")
    }

    return(x)
}

