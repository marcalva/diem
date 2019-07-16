
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
#' @param de_cutpoint Numeric. Manually initialize the DE cutpoint.
#' @param verbose Boolean. Verbosity.
#'
#' @return SCE object.
#' @importFrom DropletUtils barcodeRanks
#' @export
init_de_cutpoint <- function(x, 
							log_base=NULL, 
							de_cutpoint=NULL, 
							verbose=FALSE){
	if (verbose) cat("Setting DE cutpoint\n")

	dc <- Matrix::colSums(x@counts)
	
	# 	max_count <- max(dc)
	# 	de_cutpoint <- max_count/log_base
	#} 
	
	br.out <- barcodeRanks(x@counts[,x@test_IDs])
	de_cutpoint <- br.out$inflection

	n_over <- sum(dc > de_cutpoint)
	n_less <- sum(dc < de_cutpoint)
	if (n_over < 10){
		stop("Less than 10 droplets are over the DE cutpoint. Set to a lower value.")
	}
	if (n_less < 10){
		stop("Less than 10 droplets are under the DE cutpoint. Set to a higher value.")
	}

	if (verbose) cat(paste0("Initialized DE cutpoint at ", as.character(de_cutpoint), " counts.\n"))

	# dc <- dc[order(dc, decreasing=TRUE)]
	
	x@de_cut_init <- de_cutpoint
	
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
	iteration=NULL, 
	verbose=FALSE){
	
	if (length(x@de_cut_init) == 0){
		stop("Initialize DE cutpoint first")
	}

	if (length(x@diem) == 0){
		x@diem[[1]] <- DIEM()
	}

	if (is.null(iteration)){
		iteration <- length(x@diem)
	} else {
		if (iteration > length(x@diem)) stop("Iteration cannot be larger than the number of DIEM iterations run.")
	}

	dc <- Matrix::colSums(x@counts)

	bg_member <- rep(1, ncol(x@counts))
	names(bg_member) <- colnames(x@counts)
	tg_member <- rep(0, ncol(x@counts))
	names(tg_member) <- colnames(x@counts)

	if (iteration == 1){
		# Set hard threshold from initialization
		high_droplets <- dc > x@de_cut_init
		bg_member[high_droplets] <- 0
		tg_member[high_droplets] <- 1
	} else {
		# Use soft threshold of posterior probabilities if available from previous iteration
		Z <- x@diem[[iteration-1]]@emo[["Z"]]
		asgn <- x@diem[[iteration-1]]@emo[["assign"]]

		bg_member[rownames(Z)] <- Z[,asgn["Background"]]
		tg_member[rownames(Z)] <- Z[,asgn["Signal"]]
	}

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
	x@diem[[iteration]]@de <- DE(low_means=low_sum[c(deg_low, deg_high)],
								 high_means=high_sum[c(deg_low, deg_high)],
								 deg_low=deg_low, 
								 deg_high=deg_high,
								 log2fc=log2fc[c(deg_low, deg_high)])
	return(x)
}
