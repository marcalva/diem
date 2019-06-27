
#' Set count cutpoint for DE testing between low and high count droplets
#'
#' In order to find genes enriched in the background or the nucleus, we need to 
#' determine how to group droplets for DE. This function finds a proper cut point 
#' to decide which droplets are assigned to either the high or low group. This number 
#' is determined by dividing the maximum droplet count by \code{log_base}. 
#' Alternatively, the you can set the cutpoint threshold manually with \code{de_cutpoint}.
#'
#' @param x SCE. SCE object.
#' @param log_base Numeric. Base of log to take.
#' @param de_cutpoint Numeric. Manually set the DE cutpoint.
#' @param verbose Boolean. Verbosity.
#'
#' @return SCE object.
#' @export
set_de_cutpoint <- function(x, 
							log_base=100, 
							de_cutpoint=NULL, 
							verbose=FALSE){
	if (verbose) cat("Setting DE cutpoint\n")
	if (!is.null(de_cutpoint)){
		x@de_cutpoint <- de_cutpoint
	} else {
		max_count <- max(x@dropl_info$total_counts)
		x@de_cutpoint <- max_count/log_base
	}
	n_over <- sum(x@dropl_info$total_counts > x@de_cutpoint)
	if (n_over < 10){
		stop("Less than 10 droplets are over the DE cutpoint. Set to a lower value.\n")
	}
	if (verbose) cat(paste0("Set DE cutpoint at ", as.character(x@de_cutpoint), " counts\n"))
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
#' @return SCE object.
#' @importFrom Matrix rowMeans rowSums
#' @export
get_de <- function(x, 
	cpm_thresh=3, 
	log2fc_thresh=0.25, 
	verbose=FALSE){
	
	if (is.null(x@de_cutpoint)){
		stop("Calculate DE cutpoint first")
	}

	dc <- x@dropl_info[,"total_counts"]
	low_droplets <- (dc <= x@de_cutpoint)
	high_droplets <- (dc > x@de_cutpoint)

	genes <- rownames(x@counts)

	low_expr <- x@counts[, low_droplets]
	high_expr <- x@counts[, high_droplets]

	low_sum <- Matrix::rowSums(low_expr)
	high_sum <- Matrix::rowSums(high_expr)

	low_sum <- 1e6*low_sum/sum(low_sum)
	high_sum <- 1e6*high_sum/sum(high_sum)

	genes_expr <- genes[ low_sum > cpm_thresh & high_sum > cpm_thresh ]
	low_sum <- low_sum[genes_expr]; high_sum <- high_sum[genes_expr]
	log2fc <- log2(low_sum) - log2(high_sum)

	deg_low <- genes_expr[ log2fc > log2fc_thresh ]
	deg_high <- genes_expr[ log2fc < -log2fc_thresh ]
	deg <- c(deg_low, deg_high)

	# Save in DE
	de_obj <- DE(low_means=low_sum[deg],
				 high_means=high_sum[deg],
				 low_n=ncol(low_expr),
				 high_n=ncol(high_expr),
				 deg=deg, 
				 deg_low=deg_low, 
				 deg_high=deg_high,
				 log2fc=log2fc[deg])
	x@de <- de_obj
	return(x)
}
