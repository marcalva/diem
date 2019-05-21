
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
#' @param x SCE. SCE object
#' @param title Character
#' @param return Boolean. Return a ggplot object
#'
#' @return Nothing, unless return=TRUE then a ggplot
#' @import ggplot2
#' @importFrom scales comma
#' @export
barcode_rank_plot <- function(x, title="", return=FALSE){
	counts <- x@dropl_info[order(x@dropl_info[,"total_counts"], decreasing=TRUE), "total_counts"]
	ranks <- seq(length(counts))
	df <- data.frame(Rank=ranks, Count=counts)
	df <- df[order(df[,"Rank"]),]
	p <- ggplot(df, aes(x=Rank, y=Count)) + 
	geom_line() + 
	scale_x_continuous(trans='log10', breaks=set_breaks_10, labels=scales::comma) + 
	scale_y_continuous(name="Droplet size", trans='log10', labels=scales::comma)  + 
	theme_minimal() + theme(plot.title=element_text(hjust=0.5),
							axis.text.x=element_text(angle=45, hjust=1)) + 
	ggtitle(title)
	if (return) return(p)
	else p
}

#' Background score density plot
#'
#' This plot is helpful to see the separation of RNA content between low and high count droplets. 
#' Plot the distribution of background scores for the top 100 count-ranked droplets. 
#' Plot also contains the background score of the low-count gene probabilities 
#' used to simulate background droplets.
#'
#' @param x SCE. SCE object
#' @param title Character
#' @param return Boolean. Return a ggplot object
#'
#' @return Nothing, unless return=TRUE then a ggplot
#' @import ggplot2
#' @importFrom scales comma
#' @export
top100_score_plot <- function(x, title="", return=FALSE){
	top100 <- rownames(x@dropl_info)[order(x@dropl_info$total_counts, decreasing=TRUE)][1:100]
	hp <- sweep_total(x@raw[,top100], x@dropl_info[top100,"total_counts"])
	top100scores <- Matrix::colSums(hp[x@deg,])
	bgscore <- sum(x@gene_prob[x@deg, 1])
	df <- data.frame(Top100 = top100scores)
	p <- ggplot(df, aes(x=Top100)) + geom_density() + theme_minimal() + 
	geom_vline(xintercept=bgscore, color="darkred") + 
	xlab("Background Score")
	ggtitle(title)
	if (return) return(p)
	else p
}

#' Background score-PC plot
#'
#' Scatter plot of background scores vs PCs.
#'
#' @param x SCE. SCE object
#' @param title Character
#' @param return Boolean. Return a ggplot object
#'
#' @return Nothing, unless return=TRUE then a ggplot
#' @import ggplot2
#' @importFrom scales comma
#' @export
pc_bgscore_plot <- function(x, title="", pc=1, return=FALSE){
	df <- cbind(x@dropl_info[,c("bg_score", "Target")], x@pcs$x)
	df[df$Target,"Target"] = "Target"
	df[df$Target == "FALSE","Target"] = "Background"
	df[x@dropl_info$background == 1,"Target"] = "Simulated"
	yname <- colnames(df)[2+pc]
	p <- ggplot(df, aes_string(x="bg_score", y="PC1")) + 
	geom_point(aes(color=factor(Target))) + theme_minimal() +
	xlab("Background Score") + 
	scale_color_hue(name="Droplet") + 
	ggtitle(title)
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
