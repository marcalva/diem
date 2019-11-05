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
#' @param x An SCE object
#' @param title Title of the plot
#' @param ret Logical indicating whether to return a ggplot object
#'
#' @return Nothing. If return=TRUE, then return a ggplot object
#' @import ggplot2
#' @importFrom scales comma
#' @export
barcode_rank_plot <- function(x, title="", ret=FALSE){
	counts <- x@droplet_data[order(x@droplet_data[,"total_counts"], decreasing=TRUE), "total_counts"]
	counts <- sort(counts, decreasing=TRUE)
	ranks <- seq(length(counts))
	counts[duplicated(counts)] <- NA
	ranks[duplicated(counts)] <- NA
	df <- data.frame("Rank"=ranks, "Count"=counts)
	df <- df[!is.na(df[,2]),,drop=FALSE]
	p <- ggplot(df, aes_string(x = "Rank", y = "Count")) + 
	geom_point() + 
	scale_x_continuous(trans='log10', breaks=set_breaks_10, labels=comma) + 
	scale_y_continuous(name="Droplet size", trans='log10', labels=comma)  + 
	theme_minimal() + theme(plot.title=element_text(hjust=0.5),
							axis.text.x=element_text(angle=45, hjust=1)) + 
	ggtitle(title)
	if (ret) return(p)
	else print(p)
}

#' Scatterplot of genes vs. UMI counts, colored by posterior probability
#'
#' @param x An SCE object.
#' @param color Column name of droplet_data to color points by.
#' @param color_name Title of the color legend.
#' @param palette Name of the palette to use from RColorBrewer;
#'  see \code{\link[ggplot2]{scale_color_distiller}}.
#' @param alpha A numeric value controlling the transparency of the points. 
#'  From 0 (transparent) to 1 (no transparency).
#' @param ret A logical specifying whether to return the ggplot object 
#'  or just print it out.
#'
#' @return Nothing. If return=TRUE, then return a ggplot object
#' @import ggplot2
#' @export
plot_umi_gene <- function(x, 
                          color="CleanProb", 
                          color_name="Probability\nClean", 
                          palette="PuBuGn", 
                          alpha=0.1, 
                          ret=FALSE){

	df <- x@droplet_data[x@test_set,]

	p <- ggplot(df, aes_string(x = "total_counts", y = "n_genes")) + 
    geom_point(alpha=alpha, aes_string(colour=color)) + 
	xlab("UMI Counts") +
	ylab("Genes Detected") + 
	scale_x_log10() + scale_y_log10() + 
	theme_minimal() + theme(text=element_text(size=22)) + 
	scale_color_distiller(name=color_name, palette=palette, direction=1) 
	if (ret) return(p)
	else print(p)
}

#' Scatterplot of genes vs. UMI counts, colored by call
#'
#' @param x An SCE object.
#' @param alpha A numeric value controlling the transparency of the points. 
#'  From 0 (transparent) to 1 (no transparency).
#' @param ret A logical specifying whether to return the ggplot object 
#'  or just print it out.
#'
#' @return Nothing, unless return=TRUE then a ggplot
#' @import ggplot2
#' @export
plot_umi_gene_call <- function(x, 
                               alpha=0.1, 
                               ret=FALSE){

    color <- "Call"
    color_name <- "Call"

	df <- x@droplet_data[x@test_set,]
    df <- df[order(df$Call),]
    df$Call <- factor(df$Call, levels = c("Debris", "Clean"))

	p <- ggplot(df, aes_string(x = "total_counts", y = "n_genes")) + 
    geom_point(alpha=alpha, aes_string(colour=color)) + 
	xlab("UMI Counts") +
	ylab("Genes Detected") + 
	scale_x_log10() + scale_y_log10() + 
	theme_minimal() + theme(text=element_text(size=22)) + 
	scale_colour_brewer(name=color_name, palette = "Set1") 
	if (ret) return(p)
	else print(p)
}

