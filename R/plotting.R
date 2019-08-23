
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
	counts <- x@droplet_data[order(x@droplet_data[,"total_counts"], decreasing=TRUE), "total_counts"]
	counts <- sort(counts, decreasing=TRUE)
	ranks <- seq(length(counts))
	counts[duplicated(counts)] <- NA
	ranks[duplicated(counts)] <- NA
	df <- data.frame(Rank=ranks, Count=counts)
	df <- df[!is.na(df[,2]),,drop=FALSE]
	p <- ggplot(df, aes(x=Rank, y=Count)) + 
	geom_point() + 
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

	df <- x@droplet_data[x@test_set,]

	p <- ggplot(df, aes(x=total_counts, y=n_genes)) + geom_point(alpha=alpha, aes(colour=CleanProb)) + 
	xlab("UMI Counts") +
	ylab("Genes Detected") + 
	scale_x_log10() + scale_y_log10() + 
	theme_minimal() + theme(text=element_text(size=22)) + 
	scale_color_distiller(name="Posterior\nProbability", palette="RdBu", direction=1) 
	if (ret) return(p)
	else print(p)
}

#' Scatterplot of genes vs. UMI counts, colored by call
#'
#' @param x SCE. SCE object
#' @param alpha SCE. Transparency of points. 0 (transparent) to 1 (no transparency)
#' @param ret Boolean. Return a ggplot object
#'
#' @return Nothing, unless return=TRUE then a ggplot
#' @import ggplot2
#' @export
plot_umi_gene_call <- function(x, alpha=0.1, ret=FALSE){

	df <- x@droplet_data[x@test_set,]

	p <- ggplot(df, aes(x=total_counts, y=n_genes)) + geom_point(alpha=alpha, aes(colour=Call)) + 
	xlab("UMI Counts") +
	ylab("Genes Detected") + 
	scale_x_log10() + scale_y_log10() + 
	theme_minimal() + theme(text=element_text(size=22)) + 
	scale_colour_discrete(name="Posterior\nProbability") 
	if (ret) return(p)
	else print(p)
}

