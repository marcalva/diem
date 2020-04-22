#' @importFrom utils globalVariables

set_breaks_10 <- function(x){
	xmax <- x[2]
	bk <- 10
	brks <- c(bk)
	while (bk < xmax){
		bk <- bk * 10
		brks <- c(brks, bk)
	}
	return(brks)
}

#' Create a barcode rank plot
#'
#' @param sce An SCE object
#' @param title Title of the plot
#' @param ret Logical indicating whether to return a ggplot object
#'
#' @return Nothing. If return=TRUE, then return a ggplot object
#' @import ggplot2
#' @importFrom scales comma
#' @export
barcode_rank_plot <- function(sce, title="", ret=FALSE){
	counts <- sce@droplet_data[order(sce@droplet_data[,"total_counts"], decreasing=TRUE), "total_counts"]
	counts <- sort(counts, decreasing=TRUE)
	ranks <- seq(length(counts))
	counts[duplicated(counts)] <- NA
	ranks[duplicated(counts)] <- NA
	datf <- data.frame("Rank"=ranks, "Count"=counts)
	datf <- datf[!is.na(datf[,2]),,drop=FALSE]
	p <- ggplot(datf, aes_string(x = "Rank", y = "Count")) + 
	geom_point() + 
	scale_x_continuous(trans='log10', breaks=set_breaks_10, labels=comma) + 
	scale_y_continuous(name="Droplet size", trans='log10', labels=comma)  + 
	theme_minimal() + theme(plot.title=element_text(hjust=0.5),
							axis.text.x=element_text(angle=45, hjust=1)) + 
	ggtitle(title)
	if (ret) return(p)
	else print(p)
}

#' Scatterplot of features from meta data or gene expression
#'
#' Plot variables from meta data or from gene expression. \code{feat_x} 
#' will be plotted on the x-axis, while \code{feat_y} will be plotted on 
#' on the y-axis. These variables can come from test_data or the 
#' normalized gene expression values of a features.
#'
#' @param sce An SCE object.
#' @param feat_x Variable of x axis to plot.
#' @param feat_y Varibale of y axis to plot.
#' @param log_x Boolean indicating whether to log transform the y-axis.
#' @param log_y Boolean indicating whether to log transform the y-axis.
#' @param alpha A numeric value controlling the transparency of the points. 
#'  From 0 (transparent) to 1 (no transparency).
#' @param color Color points by values in this column of test_data, such 
#'  as Call.
#' @param data_type Either test or all, specifying the droplets to plot
#' @param ret A logical specifying whether to return the ggplot object 
#'  or just print it out.
#'
#' @return Nothing. If return=TRUE, then return a ggplot object
#' @import ggplot2
#' @export
plot_data <-  function(sce, 
                       feat_x = "total_counts", 
                       feat_y = "n_genes", 
                       log_x = FALSE, 
                       log_y = FALSE, 
                       color = NULL, 
                       alpha = 0.1, 
                       data_type = "test", 
                       ret = FALSE){

    if (data_type == "all"){
        testdat <- sce@droplet_data
    } else {
        testdat <- sce@test_data
    }

    xv <- yv <- NA
    if (feat_x %in% colnames(testdat)){
        xv <- testdat[,feat_x]
        names(xv) <- rownames(testdat)
    } else if (feat_x %in% rownames(sce@norm)){
        di <- intersect(rownames(testdat), colnames(sce@norm))
        xv <- sce@norm[feat_x, di]
        names(xv) <- di
    } else {
        stop(feat_x, " not found in test_data or genes")
    }

    if (feat_y %in% colnames(testdat)){
        yv <- testdat[,feat_y]
        names(yv) <- rownames(testdat)
    } else if (feat_y %in% rownames(sce@norm)){
        di <- intersect(rownames(testdat), colnames(sce@norm))
        yv <- sce@norm[feat_y, di]
        names(yv) <- di
    } else {
        stop(feat_y, " not found in test_data or genes")
    }

    di <- intersect(names(xv),names(yv))
    cols <- rep(NA, length(di))
    names(cols) <- di

    if (!is.null(color)){
        if (color %in% colnames(testdat)){
            cols <- testdat[,feat_y]
        } else if (color %in% rownames(sce@norm)){
            di <- intersect(rownames(testdat), colnames(sce@norm))
            cols <- sce@norm[feat_y, di]
        } else {
            if (!is.na(color)){
                stop(color, " not found in test_data or genes")
            }
        }
        di <- intersect(di, names(cols))
    }

    datf <- do.call(cbind, list(xv[di], yv[di], cols[di]))
    datf <- as.data.frame(datf)
    if (is.null(color)){
        colnames(datf) <- c(feat_x, feat_y, "color")
    } else {
        colnames(datf) <- c(feat_x, feat_y, color)
    }

	p <- ggplot(datf, aes_string(x = feat_x, y = feat_y)) + 
	theme_minimal() + 
    theme(text=element_text(size=16))

    if (log_x)
        p <- p + scale_x_log10()

    if (log_y)
        p <- p + scale_y_log10()


    if (!is.null(color)){
        p <- p + geom_point(alpha=alpha, aes_string(color = color))
    } else {
        p <- p + geom_point(alpha=alpha)
    }

	if (ret) return(p)
	else print(p)
}

#' Scatterplot of genes vs. UMI count
#'
#' @param sce An SCE object.
#' @param alpha A numeric value controlling the transparency of the points. 
#'  From 0 (transparent) to 1 (no transparency).
#' @param color Color points by values in this column of test_data, such 
#'  as Call.
#' @param ret A logical specifying whether to return the ggplot object 
#'  or just print it out.
#'
#' @return Nothing. If return=TRUE, then return a ggplot object
#' @import ggplot2
#' @export
plot_umi_gene <- function(sce, 
                          color = NULL, 
                          alpha = 0.1, 
                          ret = FALSE){

	datf <- sce@test_data

    p <- plot_data(sce, feat_x = "total_counts", feat_y = "n_genes", ret = TRUE, 
                   log_x = TRUE, log_y = TRUE, color = color, alpha = alpha)

	if (ret) return(p)
	else print(p)
}

#' Scatterplot of debris score against feature
#'
#' @param sce An SCE object.
#' @param feature Plot debris score against this feature. This can 
#'  be a column name in the test_data slot or a gene. Default is 
#'  number of genes detected
#' @param log_y Boolean indicating whether to log transform the y-axis.
#' @param alpha A numeric value controlling the transparency of the points. 
#'  From 0 (transparent) to 1 (no transparency).
#' @param color Color points by values in this column of test_data, such 
#'  as Call.
#' @param ret A logical specifying whether to return the ggplot object 
#'  or just print it out.
#'
#' @return Nothing. If return=TRUE, then return a ggplot object
#' @import ggplot2
#' @export
plot_dbr_score <-  function(sce, 
                            feature = "n_genes", 
                            log_y = TRUE, 
                            color = NULL, 
                            alpha = 0.1, 
                            ret = FALSE){

	datf <- sce@test_data

    if ( ! "score.debris" %in% colnames(datf))
        stop("Calculate debris score before plotting")

	datf <- sce@test_data

    p <- plot_data(sce, feat_x = "score.debris", feat_y = feature, ret = TRUE, 
                   log_x = FALSE, log_y = log_y, color = color, alpha = alpha)

	if (ret) return(p)
    else print(p)
}

#' Scatterplot of cluster means for indicated variables
#' 
#' Features should be a column name from droplet_data
#'
#' @param sce An SCE object.
#' @param feat_x Feature to plot on the x-axis.
#' @param feat_y Feature to plot on the y-axis.
#' @param log_x Boolean indicating whether to log transform the x-axis.
#' @param log_y Boolean indicating whether to log transform the y-axis.
#' @param ret A logical specifying whether to return the ggplot object 
#'  or just print it out.
#'
#' @return Nothing. If return=TRUE, then return a ggplot object
#' @import ggplot2
#' @export
plot_clust <- function(sce, 
                       feat_x = "score.debris", 
                       feat_y = "n_genes", 
                       log_x = TRUE, 
                       log_y = FALSE, 
                       ret = FALSE){

	testdat <- sce@test_data

    if ( ! "Cluster" %in% colnames(testdat))
        stop("Assign clusters before plotting")

    if (! feat_x %in% colnames(testdat)){
        stop(feat_x, " not found in droplet data")
    }

    if (! feat_y %in% colnames(testdat)){
        stop(feat_y, " not found in droplet data")
    }

    avg_x <- tapply(testdat[,feat_x], testdat[,"Cluster"], mean, na.rm = TRUE)
    avg_y <- tapply(testdat[,feat_y], testdat[,"Cluster"], mean, na.rm = TRUE)
    clusts <- names(avg_x)
    avg_y <- avg_y[names(avg_x)]
    avgs <- data.frame(x = avg_x, y = avg_y, cl = clusts)
    colnames(avgs) <- c(feat_x, feat_y, "Cluster")

    p <- ggplot(avgs, aes_string(x = feat_x, y = feat_y)) + 
    geom_text(aes_string(label = "Cluster")) +
    theme_minimal()
    
    if (log_x)
        p <- p + scale_x_log10()

    if (log_y)
        p <- p + scale_y_log10()

	if (ret) return(p)
    else print(p)
}

