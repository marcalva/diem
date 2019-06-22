
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
	hp <- sweep_cols(x@raw[,top100], x@dropl_info[top100,"total_counts"])
	top100scores <- Matrix::colSums(hp[x@de@deg,])
	bgscore <- sum(x@gene_prob[x@de@deg, 1])
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
#' @param ret Boolean. Return a ggplot object
#'
#' @return Nothing, unless return=TRUE then a ggplot
#' @import ggplot2
#' @importFrom scales comma
#' @export
pc_bgscore_plot <- function(x, title="", pc=1, ret=FALSE){
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
	if (ret) return(p)
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

#' Contour plot of PCs + bg_score, grouped by sim
#'
#' @param x SCE. SCE object
#' @param return Boolean. Return a ggplot object
#'
#' @return Nothing, unless return=TRUE then a ggplot
#' @import ggplot2
#' @export
contour_plot_sim <- function(sce, return=FALSE){
	pcs <- sce@pcs$x
	bg_pc_cor <- cor(sce@dropl_info$bg_score, pcs)
	pci <- which.max(bg_pc_cor)
	pcname <- colnames(pcs)[pci]
	df <- cbind(sce@dropl_info[,c("background", "Target", "bg_score")], sce@pcs$x)
	df[df$background == 1,"background"] <- "Simulated"
	df[df$background == 0,"background"] <- "Candidate"
	df[df$Target == TRUE,"Target"] <- "Nucleus"
	df[df$Target == FALSE,"Target"] <- "Debris"
	df$Border <- NA
	p <- ggplot(df, aes_string(x=pcname, y="bg_score", shape="Target")) + geom_point(alpha=0.5) + 
	geom_density2d(aes(color=background)) + theme_minimal() +
	scale_color_discrete(name="Droplet") + scale_shape(name="Call") + 
	theme(text=element_text(size=22))
	if (return) return(p)
	else p
}

#' Contour plot of PCs + bg_score, grouped by call
#'
#' @param x SCE. SCE object
#' @param return Boolean. Return a ggplot object
#'
#' @return Nothing, unless return=TRUE then a ggplot
#' @import ggplot2
#' @import reshape2
#' @export
contour_plot_call <- function(sce, n_pcs=4, return=FALSE){
	pcs <- sce@sim@pcs$x
	n_pcs <- min(ncol(pcs), n_pcs)
	pcs <- pcs[,1:n_pcs,drop=FALSE]
	df <- as.data.frame(pcs)
	df[,"Droplet"] <- "Debris"; df[targets_ids(sce),"Droplet"] <- "Nucleus"; df[grep("SIM", rownames(df)),"Droplet"] <- "Simulated"
	df[,"bg_score"] <- sce@sim@bg_score

	df <- reshape2::melt(df, id.vars = c("Droplet", "bg_score"), variable.name = "PC", value.name = "PC_score")

	p <- ggplot(df, aes(x=PC_score, y=bg_score)) + geom_point(alpha=0.1, aes(color=factor(Droplet))) + 
	geom_density2d(aes(color=factor(Droplet))) + 
	facet_wrap(~PC, scales="free") + 
	theme_minimal() + theme(text=element_text(size=22)) + 
	scale_color_discrete(name="Droplet")
	if (return) return(p)
	else p
}

#' Density plot of PCs + bg_score, grouped by call
#'
#' @param x SCE. SCE object
#' @param return Boolean. Return a ggplot object
#'
#' @return Nothing, unless return=TRUE then a ggplot
#' @import ggplot2
#' @import reshape2
#' @export
density_plot_call <- function(sce, n_pcs=4, return=FALSE){
	pcs <- sce@sim@pcs$x
	n_pcs <- min(ncol(pcs), n_pcs)
	pcs <- pcs[,1:n_pcs,drop=FALSE]
	df <- as.data.frame(pcs)
	df[,"bg_score"] <- sce@sim@bg_score
	df[,"Droplet"] <- "Debris"; df[targets_ids(sce),"Droplet"] <- "Nucleus"; df[grep("SIM", rownames(df)),"Droplet"] <- "Simulated"

	df <- reshape2::melt(df, id.vars = "Droplet", variable.name = "Feature", value.name = "Value")

	p <- ggplot(df, aes(Value, color=Droplet)) + geom_density() + 
	facet_wrap(~Feature, scales="free") + 
	theme_minimal() + theme(text=element_text(size=22)) + 
	scale_color_discrete(name="Droplet") + 
	scale_x_continuous(limits = c(-10,10))
	if (return) return(p)
	else p
}


