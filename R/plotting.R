
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

#' Volcano plot of differentially expressed genes between nuclear and background droplets
#'
#' @param x SCE. SCE object
#' @param ret Boolean. Return a ggplot object
#'
#' @return Nothing, unless return=TRUE then a ggplot
#' @import ggplot2
#' @export
volcano_de_plot <- function(x, ret=FALSE){
	df <- x@de@table[!is.na(x@de@table$log2fc),]
	df <- df[!is.infinite(df$log2fc),]
	df$p <- -10*log10(df$p)
	df$DE <- rownames(df) %in% x@de@deg

	p <- ggplot(data=df, aes(x=log2fc, y=p, color=DE)) + 
	geom_point() + 
	theme_minimal() + 
	scale_color_manual(values=c("black", "red")) + 
	xlab(expression(paste("log"[2], "(fold change) background/nuclear"))) + 
	ylab(expression(paste("-10*log"[10], "(p)"))) + 
	theme(legend.position="None")
	if (ret) return(p)
	else p
}

#' Correlation plot of low vs high expression
#'
#' @param x SCE. SCE object
#' @param ret Boolean. Return a ggplot object
#'
#' @return Nothing, unless return=TRUE then a ggplot
#' @import ggplot2
#' @importFrom Matrix colSums rowSums
#' @export
low_high_cor_plot <- function(x, scale_factor=1, ret=FALSE){
	dc <- Matrix::colSums(x@counts)
	low_p <- Matrix::rowSums(x@counts[,dc < x@de_cutpoint])
	high_p <- Matrix::rowSums(x@counts[,dc >= x@de_cutpoint])
	low_p <- log1p(scale_factor*low_p/sum(low_p))
	high_p <- log1p(scale_factor*high_p/sum(high_p))

	df <- data.frame(low=low_p,high=high_p)
	df <- df[rowSums(df) > 0,]

	df$DEG <- ""
	if ( length(x@de@deg) > 0){
		df[x@de@deg,"DEG"] <- "DEG"
	}

	df <- df[order(df$DEG),]
	limits <- c(min(min(low_p), min(high_p)), max(max(low_p), max(high_p)))

	p <- ggplot(data=df, aes(x=low, y=high, color=DEG)) + 
	geom_point() + 
	theme_minimal() + 
	xlim(limits) + ylim(limits) + 
	xlab("Gene expression below DE cutpoint") + 
	ylab("Gene expression above DE cutpoint")
	if (ret) return(p)
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
	llfs <- sort(sce@emo@Z[,"Target"], decreasing=TRUE)
	df <- data.frame(Rank=1:length(llfs), llf=llfs)
	p <- ggplot(df, aes(x=Rank, y=llf)) + 
	geom_point() + 
	theme_bw()
	if (return) return(p)
	else p
}

#' Scatterplot of pi
#'
#' @param x SCE. SCE object
#' @param return Boolean. Return a ggplot object
#'
#' @return Nothing, unless return=TRUE then a ggplot
#' @import ggplot2
#' @export
plot_pi <- function(x, color=NULL, ret=FALSE){
	df <- subset(x@dropl_info, !is.na(x@dropl_info$pi_l))

	p <- ggplot(df, aes_string(x="pi_l", y="pi_h")) + geom_point(alpha=0.5, aes_string(color=Call)) + 
	theme_minimal() + 
	xlab(expression(paste(pi[low],''))) + 
	ylab(expression(paste(pi[high],''))) + 
	theme(text=element_text(size=22))
	if (ret) return(p)
	else p
}

#' Contour plot of pi values
#'
#' @param x SCE. SCE object
#' @param return Boolean. Return a ggplot object
#'
#' @return Nothing, unless return=TRUE then a ggplot
#' @import ggplot2
#' @import reshape2
#' @export
contour_plot_pi <- function(x, ret=FALSE){
	df <- subset(x@dropl_info, !is.na(x@dropl_info$pi_l))

	p <- ggplot(df, aes(x=pi_l, y=pi_h)) + geom_point(alpha=0.1, aes(color=Call)) + 
	geom_density2d(aes(color=factor(Call))) + 
	xlab(expression(paste(pi[low],''))) +
	ylab(expression(paste(pi[high],''))) + 
	theme_minimal() + theme(text=element_text(size=22)) + 
	scale_color_discrete(name="Droplet")
	if (ret) ret(p)
	else p
}

