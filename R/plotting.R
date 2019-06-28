
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

#' Correlation plot of low vs high expression
#'
#' @param x SCE. SCE object
#' @param ret Boolean. Return a ggplot object
#'
#' @return Nothing, unless return=TRUE then a ggplot
#' @import ggplot2
#' @importFrom Matrix colSums rowSums
#' @export
de_cor_plot <- function(x, scale_factor=1e4, ret=FALSE){
	dc <- Matrix::colSums(x@counts)
	low_p <- Matrix::rowSums(x@counts[, x@low_droplets])
	high_p <- Matrix::rowSums(x@counts[, x@high_droplets])
	low_p <- log1p(scale_factor*low_p/sum(low_p))
	high_p <- log1p(scale_factor*high_p/sum(high_p))

	df <- data.frame(low=low_p,high=high_p)
	dfk <- sweep(df, 2, colSums(df), "/")
	dfk <- dfk[rowSums(dfk) > 0,]
	non_zero <- apply(dfk, 1, function(y) all(y > 0))
	kld1 <- sum(apply(dfk[non_zero,], 1, function(y) -y[1]*log(y[2]/y[1])))
	kld2 <- sum(apply(dfk[non_zero,], 1, function(y) -y[2]*log(y[1]/y[2])))
	klds <- kld1 + kld2

	pr <- round(cor(dfk[,1], dfk[,2]), 2)

	limits <- c(min(min(low_p), min(high_p)), max(max(low_p), max(high_p)))

	p <- ggplot(data=df, aes(x=low, y=high)) + 
	geom_point() + 
	theme_minimal() + 
	xlim(limits) + ylim(limits) + 
	xlab("Gene expression below DE cutpoint") + 
	ylab("Gene expression above DE cutpoint") + 
	ggtitle(paste0("R: ", as.character(pr), "; KLD: ", as.character(klds)))
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

#' Scatterplot of pi with density plots on the margin
#'
#' @param x SCE. SCE object
#' @param return Boolean. Return a ggplot object
#'
#' @return Nothing, unless return=TRUE then a ggplot
#' @import ggplot2
#' @export
plot_pi <- function(x, color=NULL, ret=FALSE){
	df <- subset(x@dropl_info, !is.na(x@dropl_info$pi_l))

	p <- ggplot(df, aes_string(x="pi_l", y="pi_h")) + geom_point(alpha=0.5, aes(color=Call)) + 
	theme_minimal() + 
	xlab(expression(paste(pi[low],''))) + 
	ylab(expression(paste(pi[high],''))) + 
	theme(text=element_text(size=22))

	if (ret) return(p)
	else p
}

#' Scatterplot of pi with density plots on the margin
#'
#' @param x SCE. SCE object
#' @param return Boolean. Return a ggplot object
#'
#' @return Nothing, unless return=TRUE then a ggplot
#' @import ggplot2
#' @importFrom grid grid.draw
#' @importFrom gridExtra grid.arrange
#' @export
plot_pi_marginal <- function(x, alpha=0.05, pdfname=NULL, w=7, h=7){
	df <- subset(x@dropl_info, !is.na(x@dropl_info$pi_l))

	common_theme <- theme_minimal() + 
	theme(axis.text=element_blank(), text=element_text(size=22))

	p1 <- ggplot(df, aes(x=pi_l)) + 
	geom_density() + 
	xlab(expression(paste(pi[low],''))) + 
	common_theme + theme(axis.title.x=element_blank(), 
		panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank())

	blank <- ggplot()+geom_point(aes(1,1), colour="white")+
         theme(axis.ticks=element_blank(), 
               panel.background=element_blank(), 
               axis.text.x=element_blank(), axis.text.y=element_blank(),           
               axis.title.x=element_blank(), axis.title.y=element_blank())

	p2 <- ggplot(df, aes_string(x="pi_l", y="pi_h")) + geom_point(alpha=alpha) + 
	xlab(expression(paste(pi[low],''))) + 
	ylab(expression(paste(pi[high],''))) + 
	common_theme

	p3 <- ggplot(df, aes(x=pi_h)) + 
	geom_density() + 
	xlab(expression(paste(pi[low],''))) + 
	coord_flip() + 
	common_theme + theme(axis.title.y=element_blank(), 
				panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank())

	if (!is.null(pdfname)){
		pdf(pdfname, width=w, height=h)
		p <- gridExtra::grid.arrange(p1, blank, p2, p3, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
		grid::grid.draw(p)
		dev.off()
	} else {
		gtable <- gridExtra::arrangeGrob(p1, blank, p2, p3, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
		return(gtable)
	}
}

#' Scatterplot of pi coloring top expressed droplets
#'
#' @param x SCE. SCE object
#' @param return Boolean. Return a ggplot object
#'
#' @return Nothing, unless return=TRUE then a ggplot
#' @import ggplot2
#' @export
plot_pi_top <- function(x, top_n=1000, color=NULL, ret=FALSE){
	df <- subset(x@dropl_info, !is.na(x@dropl_info$pi_l))
	top_drop <- rownames(df)[order(df$total_counts, decreasing=TRUE)[1:top_n]]
	df$Rank <- rep(paste0("Bottom ", as.character(nrow(df) - top_n)), nrow(df))
	df[top_drop, "Rank"] <- paste0("Top ", as.character(top_n))
	# df <- df[order(df$Rank),]

	p <- ggplot(df, aes_string(x="pi_l", y="pi_h")) + geom_point(alpha=0.5, aes(color=Rank)) + 
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

