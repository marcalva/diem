
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
de_cor_plot <- function(x, scale_factor=1e4, iteration=NULL, ret=FALSE){
	if (length(x@iter_use) == 0) {
		iter_use <- length(x@diem)
	} else {
		iter_use <- x@iter_use
	}
	if (is.null(iteration)){
		iteration <- iter_use
	}
	if ( (iteration > length(x@diem)) | (iteration < 1 ) ){
		stop("Iteration cannot be larger than the number of DIEM iterations run or smaller than 0.")

	}

	low_p <- x@diem[[iteration]]@de@low_means
	high_p <- x@diem[[iteration]]@de@high_means
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
llk_fraction_plot <- function(x, iteration=NULL, ret=FALSE){
	if (length(x@iter_use) == 0) {
		iter_use <- length(x@diem)
	} else {
		iter_use <- x@iter_use
	}
	if (is.null(iteration)){
		iteration <- iter_use
	}
	if ( (iteration > length(x@diem)) | (iteration < 1 ) ){
		stop("Iteration cannot be larger than the number of DIEM iterations run or smaller than 0.")

	}

	Z <- x@diem[[iteration]]@emo[["Z"]]
	asgn <- x@diem[[iteration]]@emo[["assign"]]
	pidf <- as.data.frame(x@diem[[iteration]]@pi)
	pi_all <- pidf[,"pi_l"] - pidf[,"pi_h"]

	ord <- order(pi_all, decreasing=FALSE)

	llfs <- Z[ord,asgn["Signal"]]
	df <- data.frame(pi_all=pi_all[ord], llf=llfs)
	p <- ggplot(df, aes(x=pi_all, y=llf)) + 
	geom_point() + 
	xlab(expression(paste(pi[low],"-",pi[high],''))) + 
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
heatmap_pi_genes <- function(x, top_n=15, iterations=NULL, ret=FALSE){
	if (is.null(iterations)){
		iterations <- 1:length(x@diem)
	}

	cors <- lapply(iterations, function(y){
		genes <- c(x@diem[[y]]@de@deg_low, x@diem[[y]]@de@deg_high)
		dcor <- cor(t(as.matrix(x@norm[genes,])), x@diem[[y]]@pi)
		dcor <- as.data.frame(dcor)
		dcormax <- abs(apply(dcor, 1, max))
		# dcor <- dcor[order(dcormax, decreasing=TRUE),,drop=FALSE]
		dcor[,"Iteration"] <- rep(y, nrow(dcor))
		dcor[,"Gene"] <- rownames(dcor)
		return(dcor)
	})

	# Sort by last dcor
	dcor <- cors[[x@iter_use]]
	genesl <- rownames(dcor)[order(dcor[,"pi_l"], decreasing=TRUE)][1:top_n]
	genesh <- rownames(dcor)[order(dcor[,"pi_h"], decreasing=TRUE)][1:top_n]
	genes <- c(genesl, rev(genesh))
	genes_ord <- genes

	corsdf <- do.call(rbind, cors)
	corsdf <- corsdf[corsdf[,"Gene"] %in% genes,,drop=FALSE]
	corsdf <- reshape2::melt(corsdf, id.vars=c("Iteration", "Gene"))


    hmcol = colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(100)
    ylab <- c(expression(paste(pi[low],'')),
    		  expression(paste(pi[high],'')))

    p <- ggplot(data=corsdf, aes(x=variable, y=Gene, fill=value)) + 
    facet_wrap(~Iteration, nrow=1, strip.position="bottom") + 
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
#' @param ret Boolean. Return a ggplot object
#'
#' @return Nothing, unless return=TRUE then a ggplot
#' @import ggplot2
#' @export
plot_pi <- function(x, color=NULL, alpha=0.1, iteration=NULL, ret=FALSE){
	if (length(x@iter_use) == 0) {
		iter_use <- length(x@diem)
	} else {
		iter_use <- x@iter_use
	}
	if (is.null(iteration)){
		iteration <- iter_use
	}
	if ( (iteration > length(x@diem)) | (iteration < 1 ) ){
		stop("Iteration cannot be larger than the number of DIEM iterations run or smaller than 0.")

	}

	df <- as.data.frame(x@diem[[iteration]]@pi)
	df[,"Call"] <- x@diem[[iteration]]@calls
	di <- x@dropl_info[rownames(df),]
	di <- di[, !(colnames(di) %in% colnames(df))]
	df <- cbind(df, di)


	p <- ggplot(df, aes(x=pi_l, y=pi_h)) + geom_point(alpha=alpha, aes(color=Call)) + 
	xlab(expression(paste(pi[low],''))) +
	ylab(expression(paste(pi[high],''))) + 
	theme_minimal() + theme(text=element_text(size=22), axis.text=element_blank()) + 
	scale_color_discrete(name="Droplet") + 
	guides(color = guide_legend(override.aes = list(alpha = 1)))
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
plot_pi_marginal <- function(x, alpha=0.05, pdfname="pi_marginal.pdf", iteration=NULL, w=7, h=7){
	if (length(x@iter_use) == 0) {
		iter_use <- length(x@diem)
	} else {
		iter_use <- x@iter_use
	}
	if (is.null(iteration)){
		iteration <- iter_use
	}
	if ( (iteration > length(x@diem)) | (iteration < 1 ) ){
		stop("Iteration cannot be larger than the number of DIEM iterations run or smaller than 0.")

	}

	df <- as.data.frame(x@diem[[iteration]]@pi)

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
plot_pi_tails <- function(x, pct=0.05, iteration=NULL, ret=FALSE){
	if (length(x@iter_use) == 0) {
		iter_use <- length(x@diem)
	} else {
		iter_use <- x@iter_use
	}
	if (is.null(iteration)){
		iteration <- iter_use
	}
	if ( (iteration > length(x@diem)) | (iteration < 1 ) ){
		stop("Iteration cannot be larger than the number of DIEM iterations run or smaller than 0.")

	}

	df <- as.data.frame(x@diem[[iteration]]@pi)
	df[,"Call"] <- x@diem[[iteration]]@calls
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

#' Contour plot of MVN densities for the 2 groups
#'
#' @param x SCE. SCE object
#' @param ret Boolean. Return a ggplot object
#'
#' @return Nothing, unless return=TRUE then a ggplot
#' @import ggplot2
#' @import reshape2
#' @export
plot_mvn_contour <- function(x, length.out=1e2, iteration=NULL, ret=FALSE){
	if (length(x@iter_use) == 0) {
		iter_use <- length(x@diem)
	} else {
		iter_use <- x@iter_use
	}
	if (is.null(iteration)){
		iteration <- iter_use
	}
	if ( (iteration > length(x@diem)) | (iteration < 1 ) ){
		stop("Iteration cannot be larger than the number of DIEM iterations run or smaller than 0.")

	}

	mins <- apply(x@diem[[iteration]]@pi, 2, min)
	maxs <- apply(x@diem[[iteration]]@pi, 2, max)
	dl1 <- (maxs[1] - mins[1])/10
	dl2 <- (maxs[2] - mins[2])/10
	mins[1] <- mins[1]-dl1; mins[2] <- mins[2]-dl2
	maxs[1] <- maxs[1]+dl1; maxs[2] <- maxs[2]+dl2

	x1seq <- seq(mins[1], maxs[1], length.out=length.out)
	x2seq <- seq(mins[2], maxs[2], length.out=length.out)
	z1 <- matrix(nrow=length.out, ncol=length.out)
	colnames(z1) <- x2seq; rownames(z1) <- x1seq
	z2 <- matrix(nrow=length.out, ncol=length.out)
	colnames(z2) <- x2seq; rownames(z2) <- x1seq

	emo <- x@diem[[iteration]]@emo

	for (i in 1:length(x1seq)){
		for (j in 1:length(x2seq)){
			z1[i,j] <- EMCluster::dmvn(c(x1seq[i],x2seq[j]), 
					mu=emo[["Mu"]][1,], 
					LTsigma=emo[["LTSigma"]][1,])
			z2[i,j] <- EMCluster::dmvn(c(x1seq[i],x2seq[j]), 
					mu=emo[["Mu"]][2,], 
					LTsigma=emo[["LTSigma"]][2,])

		}
	}
	c1m <- reshape2::melt(z1)
	c1m[,"Droplet"] <- "Signal"
	c2m <- reshape2::melt(z2)
	c2m[,"Droplet"] <- "Background"
	cm <- rbind(c1m, c2m)

	p <- ggplot(cm, aes(x=Var1, y=Var2, z=value, color=Droplet)) + 
	geom_contour() + 
	xlim(mins[1]-dl1, maxs[1]+dl1) + 
	ylim(mins[2]-dl2, maxs[2]+dl2) + 
	xlab(expression(paste(pi[low],''))) + 
	ylab(expression(paste(pi[high],''))) + 
	theme_minimal()

	if (ret) return(p)
	else print(p)
}

#' Plot through iterations, return a list of ggplot2 objects
#'
#' @param x SCE. SCE object
#' @param FUN character string. Name of function for plotting. Can be one of 
#'  of the following:
#'  \itemize{
#'		\item de_cor_plot
#'		\item llk_fraction_plot
#'		\item plot_pi
#'		\item plot_mvn_contour
#'  }
#' @param ret Logical. If TRUE, return the list of ggplots, else, print them
#'
#' @return Nothing
#' @export
plot_iter <- function(x, FUN, ret=FALSE, ...){
	if (length(x@name) == 0) x@name <- ""
	FUN <- match.fun(FUN)
	to_ret <- list()
	for (iteration in 1:length(x@diem)){
		p <- FUN(x, iteration=iteration, ret=TRUE, ...) + ggtitle(paste0(x@name, " Iteration ", as.character(iteration)))
		to_ret[[iteration]] <- p
	}
	if (ret) return(to_ret)
	else{
		for (i in 1:length(to_ret)){
			print(to_ret[[i]])
		}
	}
}

#' Plot through iterations in a gif
#'
#' @param x SCE. SCE object
#' @param FUN character string. Name of function for plotting. See \code{plot_iter} for more details.
#' @param prefix character string. Prefix of GIF output file name.
#'
#' @return Nothing
#' @export
plot_iter_gif <- function(x, FUN, prefix="iter", ...){
	png(paste0(prefix, ".tmp.%03d.png"))
	FUN <- match.fun(FUN)
	plots <- plot_iter(x, FUN, ret=TRUE, ...)
	for (p in plots){ print(p) }
	dev.off()
	system(paste0("convert -delay 50 ", prefix, ".tmp.*.png", " ", prefix, ".gif"))
	system(paste0("rm ", prefix, ".tmp.*.png"))
}
