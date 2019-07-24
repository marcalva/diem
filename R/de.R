
bin_by_counts <- function(counts, min_cell=10, log_base=4){
	tc <- Matrix::colSums(counts)
	min_c <- min(tc)
	max_c <- max(tc)
	min_l <- floor(log(min_c, base=log_base))
	max_l <- ceiling(log(max_c, base=log_base))
	limits <- log_base^(min_l:max_l)
	n_bins <- length(limits) - 1
	bins <- lapply(1:n_bins, function(i) c(limits[i], limits[i+1]))

	# Remove bins with 0 counts
	keep <- sapply(bins, function(bin){
		keep <- (tc > bin[1]) & (tc <= bin[2])
		sum(keep) >= min_cell
		})
	bins <- bins[keep]
	n_bins <- length(bins)

	bin_names <- sapply(bins, function(bin) {paste0(as.character(bin[1]), "-", as.character(bin[2]))})
	bin_counts <- matrix(nrow=nrow(counts), ncol=n_bins)
	rownames(bin_counts) <- rownames(counts); colnames(bin_counts) <- bin_names
	for (i in seq_along(bins)){
		bin <- bins[[i]]
		keep <- (tc > bin[1]) & (tc <= bin[2])
		bin_counts[,i] <- Matrix::rowSums(counts[,keep,drop=FALSE])
	}
	return(list(counts=bin_counts, bins=bins))
}

# get slope
get_slope <- function(X,y){
	betas <- apply(X, 1, function(xi){
		(t(y)%*%xi)/(t(xi)%*%xi)
		})
	return(betas)
}

#' Set count cutpoint for DE testing between low and high count droplets
#'
#' In order to find genes enriched in the background or the signal, we need to 
#' determine how to group droplets for DE. This function initializes a cut point 
#' to decide which droplets are assigned to either the high or low group. This number 
#' is determined by dividing the maximum droplet count by \code{log_base}. 
#' Alternatively, the you can set the cutpoint threshold manually with \code{de_cutpoint}.
#'
#' @param x SCE. SCE object.
#' @param log_base Numeric. Base of log to take.
#' @param diff_thresh Numeric. Manually initialize the DE cutpoint.
#' @param verbose Boolean. Verbosity.
#'
#' @return SCE object.
#' @importFrom DropletUtils barcodeRanks
#' @export
get_diff_thresh <- function(x, 
							top_n=1e4, 
							diff_thresh=NULL, 
							verbose=FALSE){
	if (verbose) cat("Setting DE cutpoint\n")

	dc <- Matrix::colSums(x@counts)
	
	if (is.null(diff_thresh)){
		br.out <- barcodeRanks(x@counts[,x@test_IDs])
		diff_thresh <- br.out$inflection
	}
	#if (is.null(diff_thresh)){
	#	diff_thresh <- (dc[order(dc, decreasing=TRUE)])[top_n]
	#}

	n_over <- sum(dc > diff_thresh)
	if (n_over < 10){
		stop("Less than 10 droplets are over the DE threshold. Set to a lower value.")
	}
	n_less <- sum(dc < diff_thresh)
	if (n_less < 10){
		stop("Less than 10 droplets are under the DE threshold. Set to a higher value.")
	}

	if (verbose) cat(paste0("Set DE cutpoint at ", as.character(diff_thresh), " counts.\n"))

	x@diff_thresh <- diff_thresh
	return(x)
}

#' Set differential genes between high and low read count droplets
#'
#' Find differentially expressed genes between low and high read count droplets. Groups droplets into
#' low or high count groups. The count threshold for separating low and high is determined by taking the 
#' log10 of the maximum droplet count, subtracting 1, and calculating 10 to the power of this number. In other 
#' words, the thresold is found by dividing the max count by 10. 
#'
#' @param x SCE. SCE object.
#' @param n_genes Numeric. Number of genes to set as differential.
#' @param cpm_thresh Numeric. CPM treshold for expressed genes in both low and high groups.
#' @param verbose Boolean.
#'
#' @return DE object.
#' @importFrom Matrix rowMeans rowSums
#' @export
set_diff_genes <- function(x, n_genes=100, cpm_thresh=3, verbose=FALSE){

	if (is.null(n_genes)){
		n_genes <- nrow(x@counts)
	}
	
	if (length(x@diff_thresh) == 0){
		x <- get_diff_thresh(x)
	}

	dc <- Matrix::colSums(x@counts)

	bg_member <- rep(1, ncol(x@counts))
	names(bg_member) <- colnames(x@counts)
	tg_member <- rep(0, ncol(x@counts))
	names(tg_member) <- colnames(x@counts)

	high_droplets <- dc > x@diff_thresh
	bg_member[high_droplets] <- 0
	tg_member[high_droplets] <- 1

	# Get props and cpm in each group
	low_sum <- x@counts %*% bg_member; low_sum <- low_sum[,1]
	low_sum <- low_sum/sum(low_sum)
	low_cpm <- 1e6*low_sum
	high_sum <- x@counts %*% tg_member; high_sum <- high_sum[,1]
	high_sum <- high_sum/sum(high_sum)
	high_cpm <- 1e6*high_sum

	# Subset to genes expressed
	genes <- rownames(x@counts)
	genes_expr <- genes[ low_cpm > cpm_thresh & high_cpm > cpm_thresh ]
	if (n_genes <= 0 ){
		stop(paste0("n_genes must be greater than 0."))
	}
	low_sum <- low_sum[genes_expr]; high_sum <- high_sum[genes_expr]
	low_cpm <- low_cpm[genes_expr]; high_cpm <- high_cpm[genes_expr]

	prop_diff <- low_sum - high_sum

	diff_genes <- names(prop_diff)[order(abs(prop_diff), decreasing=TRUE)][1:n_genes]
	diff_low <- prop_diff[diff_genes][prop_diff[diff_genes] > 0]
	diff_high <- prop_diff[diff_genes][prop_diff[diff_genes] < 0]

	x@diem@de <- DE(low_prop=low_sum,
					high_prop=high_sum,
					diff_genes=diff_genes,
					prop_diff=prop_diff)
	return(x)
}

#' Run differential expression between high and low read count droplets
#'
#' Find differentially expressed genes between low and high read count droplets. Groups droplets into
#' low or high count groups. The count threshold for separating low and high is determined by taking the 
#' log10 of the maximum droplet count, subtracting 1, and calculating 10 to the power of this number. In other 
#' words, the thresold is found by dividing the max count by 10. 
#'
#' @param x SCE. SCE object.
#' @param cpm_thresh Numeric. CPM treshold for expressed genes in both low and high groups.
#' @param log2fc_thresh Numeric. Call DE genes with log2 fold change above or below this threshold.
#' @param verbose Boolean.
#'
#' @return DE object.
#' @importFrom Matrix rowMeans rowSums
#' @export
get_log2fc <- function(x, 
					   cpm_thresh=3, 
					   log2fc_thresh=0, 
					   verbose=FALSE){

	if (length(x@de_cut) == 0){
		stop("Initialize DE cutpoint first")
	}

	dc <- Matrix::colSums(x@counts)

	bg_member <- rep(1, ncol(x@counts))
	names(bg_member) <- colnames(x@counts)
	tg_member <- rep(0, ncol(x@counts))
	names(tg_member) <- colnames(x@counts)

	# Set hard threshold from initialization
	high_droplets <- dc > x@de_cut
	bg_member[high_droplets] <- 0
	tg_member[high_droplets] <- 1

	low_sum <- x@counts %*% bg_member; low_sum <- low_sum[,1]
	low_sum <- low_sum/sum(low_sum)
	low_cpm <- 1e6*low_sum

	high_sum <- x@counts %*% tg_member; high_sum <- high_sum[,1]
	high_sum <- high_sum/sum(high_sum)
	high_cpm <- 1e6*high_sum

	genes <- rownames(x@counts)
	genes_expr <- genes[ low_cpm > cpm_thresh & high_cpm > cpm_thresh ]

	low_sum <- low_sum[genes_expr]; high_sum <- high_sum[genes_expr]
	low_cpm <- low_cpm[genes_expr]; high_cpm <- high_cpm[genes_expr]
	log2fc <- log2(low_cpm) - log2(high_cpm)
	names(log2fc) <- genes_expr

	deg_low <- genes_expr[ log2fc > log2fc_thresh ]
	deg_high <- genes_expr[ log2fc < -log2fc_thresh ]

	# Save in DE
	x@diem@de <- DE(low_means=low_sum[c(deg_low, deg_high)],
					high_means=high_sum[c(deg_low, deg_high)],
					deg_low=deg_low, 
					deg_high=deg_high,
					log2fc=log2fc[c(deg_low, deg_high)])
	return(x)
}
