
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
#' @param ret Boolean. Return a ggplot object
#'
#' @return Nothing, unless return=TRUE then a ggplot
#' @import ggplot2
#' @importFrom scales comma
#' @export
barcode_rank_plot <- function(x, title="", ret=FALSE){
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
	if (ret) return(p)
	else print(p)
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
de_cor_plot <- function(x, scale_factor=1e4, ret=FALSE, main=NULL){
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
	klds <- round(kld1 + kld2, 2)

	pr <- round(cor(dfk[,1], dfk[,2]), 2)

	anno <- paste0("R: ", as.character(pr), "\nKLD: ", as.character(klds))

	limits <- c(min(min(low_p), min(high_p)), max(max(low_p), max(high_p)))

    if (is.null(main)){
        if (is.null(x@name)){
            main = ""
        } else {
            main <- x@name
        }
    }

	p <- ggplot(data=df, aes(x=low, y=high)) + 
	geom_point() + 
	theme_minimal() + 
	xlim(limits) + ylim(limits) + 
	xlab("Debris-enriched") + 
	ylab("Nuclear-enriched") + 
	annotate("text", x=limits[1], y=limits[2], label=anno, hjust=0, vjust=1) + 
	ggtitle(main)
	if (ret) return(p)
	else print(p)
}

#' Distribution of likelihood fraction
#'
#' @param x SCE. SCE object
#' @param ret Boolean. Return a ggplot object
#'
#' @return Nothing, unless return=TRUE then a ggplot
#' @import ggplot2
#' @export
llk_fraction_plot <- function(x, ret=FALSE){
	llfs <- sort(sce@emo@Z[,x@diem@emo@assign["Target"]], decreasing=TRUE)
	df <- data.frame(Rank=1:length(llfs), llf=llfs)
	p <- ggplot(df, aes(x=Rank, y=llf)) + 
	geom_point() + 
	theme_bw()
	if (ret) return(p)
	else print(p)
}

#' Scatterplot of pi with density plots on the margin
#'
#' @param x SCE. SCE object
#' @param ret Boolean. Return a ggplot object
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
	theme(text=element_text(size=22), axis.text=element_blank())

	if (ret) return(p)
	else print(p)
}

#' Scatterplot of pi with density plots on the margin
#'
#' @param x SCE. SCE object
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

#' Scatterplot of pi colored by tails
#'
#' @param x SCE. SCE object
#' @param ret Boolean. Return a ggplot object
#'
#' @return Nothing, unless return=TRUE then a ggplot
#' @import ggplot2
#' @export
plot_pi_tails <- function(x, pct=5, ret=FALSE){
	df <- x@dropl_info
	df <- df[!is.na(df$Call),]
	top_n <- floor(nrow(df)*pct/100)
	top_drop <- rownames(df)[order(df$total_counts, decreasing=TRUE)[1:top_n]]
	bot_drop <- rownames(df)[order(df$total_counts, decreasing=FALSE)[1:top_n]]
	df$Rank <- rep("", nrow(df))
	df[top_drop, "Rank"] <- paste0("Top ", as.character(pct), "%")
	df[bot_drop, "Rank"] <- paste0("Bottom ", as.character(pct), "%")
	df <- df[order(df$Rank),]
	df <- df[df$Rank != "",]

	p <- ggplot(df, aes_string(x="pi_l", y="pi_h")) + 
	geom_point(alpha=0.5, aes(color=Rank)) + 
	theme_minimal() + 
	scale_color_discrete(name="Rank by\ncounts") + 
	xlab(expression(paste(pi[low],''))) + 
	ylab(expression(paste(pi[high],''))) + 
	theme(text=element_text(size=22), axis.text=element_blank())

	if (ret) return(p)
	else print(p)
}

#' Contour plot of pi values
#'
#' @param x SCE. SCE object
#' @param ret Boolean. Return a ggplot object
#'
#' @return Nothing, unless return=TRUE then a ggplot
#' @import ggplot2
#' @import reshape2
#' @export
plot_pi_contour <- function(x, ret=FALSE){
	df <- subset(x@dropl_info, !is.na(x@dropl_info$pi_l))

	p <- ggplot(df, aes(x=pi_l, y=pi_h)) + geom_point(alpha=0.1, aes(color=Call)) + 
	geom_density2d(aes(color=factor(Call))) + 
	xlab(expression(paste(pi[low],''))) +
	ylab(expression(paste(pi[high],''))) + 
	theme_minimal() + theme(text=element_text(size=22), axis.text=element_blank()) + 
	scale_color_discrete(name="Droplet")
	if (ret) return(p)
	else print(p)
}

#' Contour plot of pi values through iterations
#'
#' @param x SCE. SCE object
#' @param ret Boolean. Return a ggplot object
#'
#' @return Nothing, unless return=TRUE then a ggplot
#' @import ggplot2
#' @import reshape2
#' @export
plot_pi_contour_iter <- function(x, ret=FALSE){
	to_ret <- list()
	for (i in 1:length(x@prev_iter)){
		df <- x@prev_iter[[i]]
		df <- df[!is.na(df$pi_l), ]

		p <- ggplot(df, aes(x=pi_l, y=pi_h)) + geom_point(alpha=0.1, aes(color=Call)) + 
		geom_density2d(aes(color=factor(Call))) + 
		xlab(expression(paste(pi[low],''))) +
		ylab(expression(paste(pi[high],''))) + 
		theme_minimal() + theme(text=element_text(size=22), axis.text=element_blank()) + 
		scale_color_discrete(name="Droplet")
		to_ret[[i]] <- p
	}
	if (ret) return(to_ret)
	else{
		for (i in 1:length(to_ret)){
			print(to_ret[[i]])
		}
	}
}

#' Plot of pi values through iterations in a gif
#'
#' @param x SCE. SCE object
#'
#' @return Nothing
#' @export
plot_pi_iter_gif <- function(x, prefix="pi_iter"){
	png(paste0(prefix, ".tmp.%03d.png"))
	plots <- plot_pi_contour_iter(x, ret=TRUE)
	for (p in plots){ print(p) }
	dev.off()
	system(paste0("convert -delay 50 ", prefix, ".tmp.*.png", " ", prefix, ".gif"))
	system(paste0("rm ", prefix, ".tmp.*.png"))
}

