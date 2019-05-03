
#' Create a barcode rank plot
#'
#' @param x SCE. SCE object
#' @param return Boolean. Return a ggplot object
#'
#' @return Nothing, unless return=TRUE then a ggplot
#' @import ggplot2
#' @export
barcode_rank_plot <- function(x, return=FALSE){
	if ( is.null(x@gene_counts) )
		x <- fill_counts(x)
	counts <- x@dropl_counts[ order(x@dropl_counts, decreasing=TRUE) ]
	ranks <- seq(length(counts))
	df <- data.frame(Rank=ranks, Count=counts)
	df <- df[order(df[,"Rank"]),]
	p <- ggplot(df, aes(x=Rank, y=Count)) + 
	geom_point() + 
	scale_x_continuous(trans='log10') + 
	scale_y_continuous(trans='log10')  + 
	theme_bw()
	if (return) return(p)
	else p
}

#' Distribution of likelihood fraction
#'
#' @param x SCE. SCE object
#' @param return Boolean. Return a ggplot object
#'
#' @return Nothing, unless return=TRUE then a ggplot
#' @import ggplot2
#' @export
llk_fraction_plot <- function(x, return=FALSE){
	llfs <- sort(sce@emo$Z[,"Target"], decreasing=TRUE)
	df <- data.frame(Rank=1:length(llfs), llf=llfs)
	p <- ggplot(df, aes(x=Rank, y=llf)) + 
	geom_point() + 
	theme_bw()
	if (return) return(p)
	else p
}
