
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
de_cor_plot <- function(x, scale_factor=1e4, ret=FALSE){

	low_p <- x@diem[[iteration]]@de@low_prop
	high_p <- x@diem[[iteration]]@de@high_prop
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

	p <- ggplot(data=df, aes(x=low, y=high)) + 
	geom_point() + 
	theme_minimal() + 
	xlim(limits) + ylim(limits) + 
	xlab("Debris") + 
	ylab("Signal") + 
	annotate("text", x=limits[1], y=limits[2], label=anno, hjust=0, vjust=1)
	if (ret) return(p)
	else print(p)
}

#' Plot posterior probability against pi_low - pi_high
#'
#' @param x SCE. SCE object
#' @param ret Boolean. Return a ggplot object
#'
#' @return Nothing, unless return=TRUE then a ggplot
#' @import ggplot2
#' @export
pp_plot <- function(x, ret=FALSE){

	Z <- x@diem@emo[["Z"]]
	asgn <- x@diem@emo[["assign"]]
	pidf <- as.data.frame(x@diem@pi)
	llfs <- sort(Z[,asgn["Signal"]])
	df <- data.frame(Rank=1:length(llfs), llf=llfs)
	p <- ggplot(df, aes(x=Rank, y=llf)) + 
	geom_point() + 
	xlab("Rank") + 
	ylab("Posterior Probability") + 
	theme_bw()
	if (ret) return(p)
	else print(p)
}

#' Heatmap of pi-gene expression correlations
#'
#' @param x SCE. SCE object
#' @param ret Boolean. Return a ggplot object
#'
#' @return Nothing, unless ret=TRUE then a ggplot
#' @import ggplot2
#' @export
heatmap_pi_genes <- function(x, top_n=15, ret=FALSE){



	genes <- x@diem@DE@diff_genes
	dcor <- cor(t(as.matrix(x@norm[genes,])), x@diem@pi)
	dcor <- as.data.frame(dcor)
	dcormax <- abs(apply(dcor, 1, max))
	dcor[,"Gene"] <- rownames(dcor)

	# Sort by last dcor
	genesl <- rownames(dcor)[order(dcor[,"pi_l"], decreasing=TRUE)][1:top_n]
	genesh <- rownames(dcor)[order(dcor[,"pi_h"], decreasing=TRUE)][1:top_n]
	genes <- c(genesl, rev(genesh))
	genes_ord <- genes

	corsdf <- do.call(rbind, cors)
	corsdf <- corsdf[corsdf[,"Gene"] %in% genes,,drop=FALSE]
	corsdf <- reshape2::melt(corsdf, id.vars="Gene")


    hmcol = colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(100)
    ylab <- c(expression(paste(pi[low],'')),
    		  expression(paste(pi[high],'')))

    p <- ggplot(data=corsdf, aes(x=variable, y=Gene, fill=value)) + 
	geom_tile() + 
    theme_classic() +
    scale_fill_gradientn(colours=hmcol, limits=c(-1, 1), name=expression(italic("R"))) + 
    scale_x_discrete(position = "top", labels=ylab)  + 
    scale_y_discrete(limits=genes) + 
    xlab("") +
    theme(text=element_text(size=18), 
    	  axis.text.x=element_text(angle=90, size=22, vjust=0.5, hjust=0),
          axis.ticks=element_blank(),
          axis.line=element_blank())
    
    if (ret) return(p)
	else print(p)
}

#' Scatterplot of pi
#'
#' @param x SCE. SCE object
#' @param alpha SCE. Transparency of points. 0 (transparent) to 1 (no transparency)
#' @param ret Boolean. Return a ggplot object
#'
#' @return Nothing, unless return=TRUE then a ggplot
#' @import ggplot2
#' @export
plot_pi_call <- function(x, alpha=0.1, ret=FALSE){

	df <- as.data.frame(x@diem@pi)
	df <- df[x@test_IDs,]
	df[names(x@diem@calls),"Call"] <- x@diem@calls
	di <- x@dropl_info[rownames(df),]
	di <- di[, !(colnames(di) %in% colnames(df))]
	df <- cbind(df, di)


	p <- ggplot(df, aes(x=pi_l, y=pi_h)) + geom_point(alpha=alpha, aes(color="Call")) + 
	xlab(expression(paste(pi[low],''))) +
	ylab(expression(paste(pi[high],''))) + 
	theme_minimal() + theme(text=element_text(size=22), axis.text=element_blank()) + 
	scale_color_discrete(name="Call") + 
	guides(color = guide_legend(override.aes = list(alpha = 1)))
	if (ret) return(p)
	else print(p)
}

#' Scatterplot of pi colored by posterior probability
#'
#' @param x SCE. SCE object
#' @param alpha SCE. Transparency of points. 0 (transparent) to 1 (no transparency)
#' @param ret Boolean. Return a ggplot object
#'
#' @return Nothing, unless return=TRUE then a ggplot
#' @import ggplot2
#' @export
plot_pi_pp <- function(x, alpha=0.1, ret=FALSE){

	df <- as.data.frame(x@diem@pi)
	df <- df[x@test_IDs,]
	df[names(x@diem@calls),"Call"] <- x@diem@calls
	di <- x@dropl_info[rownames(df),]
	di <- di[, !(colnames(di) %in% colnames(df))]
	df <- cbind(df, di)

	p <- ggplot(df, aes(x=pi_l, y=pi_h)) + geom_point(alpha=alpha, aes(colour=PP)) + 
	xlab(expression(paste(pi[low],''))) +
	ylab(expression(paste(pi[high],''))) + 
	theme_minimal() + theme(text=element_text(size=22), axis.text=element_blank()) + 
	scale_color_distiller(name="Posterior\nProbability", palette="RdBu", direction=1) 
	if (ret) return(p)
	else print(p)
}

#' Scatterplot of genes vs. UMI counts, colored by posterior probability
#'
#' @param x SCE. SCE object
#' @param alpha SCE. Transparency of points. 0 (transparent) to 1 (no transparency)
#' @param ret Boolean. Return a ggplot object
#'
#' @return Nothing, unless return=TRUE then a ggplot
#' @import ggplot2
#' @export
plot_umi_gene_pp <- function(x, alpha=0.1, ret=FALSE){

	df <- x@dropl_info[x@test_IDs,]

	p <- ggplot(df, aes(x=total_counts, y=n_genes)) + geom_point(alpha=alpha, aes(colour=PP)) + 
	xlab("UMI Counts") +
	ylab("Genes Detected") + 
	scale_x_log10() + scale_y_log10() + 
	theme_minimal() + theme(text=element_text(size=22)) + 
	scale_color_distiller(name="Posterior\nProbability", palette="RdBu", direction=1) 
	if (ret) return(p)
	else print(p)
}

#' Scatterplot of pi colored by posterior probability
#'
#' @param x SCE. SCE object
#' @param alpha SCE. Transparency of points. 0 (transparent) to 1 (no transparency)
#' @param ret Boolean. Return a ggplot object
#'
#' @return Nothing, unless return=TRUE then a ggplot
#' @import ggplot2
#' @export
plot_scatter_pp <- function(sce, x, y, alpha=0.1, ret=FALSE){

	df <- as.data.frame(sce@diem@pi)
	df <- df[sce@test_IDs,]
	df[names(sce@diem@calls),"Call"] <- sce@diem@calls
	di <- sce@dropl_info[rownames(df),]
	di <- di[, !(colnames(di) %in% colnames(df))]
	df <- cbind(df, di)

	p <- ggplot(df, aes_string(x=x, y=y)) + geom_point(alpha=alpha, aes(colour=PP)) + 
	theme_minimal() + theme(text=element_text(size=22)) + 
	scale_color_distiller(name="Posterior\nProbability", palette="RdBu", direction=1) 
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
plot_pi_marginal <- function(x, alpha=0.05, pdfname="pi_marginal.pdf", w=7, h=7){

	df <- as.data.frame(x@diem@pi)

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
	common_theme + theme(axis.title.y=element_blank(), 
				panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank())

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
plot_pi_tails <- function(x, pct=0.05, ret=FALSE){

	df <- as.data.frame(x@diem@pi)
	df[,"Call"] <- x@diem@calls
	di <- x@dropl_info[rownames(df),]
	di <- di[, !(colnames(di) %in% colnames(df))]
	df <- cbind(df, di)

	top_n <- floor(nrow(df)*pct)
	top_drop <- rownames(df)[order(df$total_counts, decreasing=TRUE)[1:top_n]]
	bot_drop <- rownames(df)[order(df$total_counts, decreasing=FALSE)[1:top_n]]
	df$Rank <- rep("", nrow(df))
	df[top_drop, "Rank"] <- paste0("Top ", as.character(100*pct), "%")
	df[bot_drop, "Rank"] <- paste0("Bottom ", as.character(100*pct), "%")
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
