
#' Set low and high droplet labels
#'
#' Set the droplet IDs to be used as the low and high groups for 
#' DE testing. By default, the high group is set as the nuclei called 
#' by the previous iteration of DIEM, while the low group is set as 
#' everything else. If \code{low} and \code{high} are given, as in 
#' during initialization when nuclei haven't been called yet, the low 
#' and high groups are set to these for logFC calculation.
#'
#' @param x SCE. SCE object.
#' @param low Character. Droplet IDs to set as the low group.
#' @param high Character. Droplet IDs to set as the high group.
#'
#' @return SCE object.
#' @export
set_de_groups <- function(x, low=NULL, high=NULL){

	# Set DE from calls in droplet info
	if (is.null(low)){
		low <- rownames(x@dropl_info)[x@dropl_info[,"Call"] %in% "Debris"]
	}
	if (is.null(high)){
		high <- rownames(x@dropl_info)[x@dropl_info[,"Call"] %in% "Nucleus"]
	}

	# Incorporate droplets below test set into background
	dc <- x@dropl_info[,"total_counts"]
	names(dc) <- rownames(x@dropl_info)
	min_count <- min(dc[x@test_droplets])

	below_test <- names(dc)[dc < min_count]

	x@low_droplets <- union(low, below_test)
	x@high_droplets <- high
	return(x)
}

#' Set count cutpoint for DE testing between low and high count droplets
#'
#' In order to find genes enriched in the background or the nucleus, we need to 
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
#' @export
init_de_cutpoint <- function(x, 
							log_base=50, 
							de_cutpoint=NULL, 
							verbose=FALSE){
	if (verbose) cat("Setting DE cutpoint\n")

	dc <- x@dropl_info[,"total_counts"]
	names(dc) <- rownames(x@dropl_info)

	if (is.null(de_cutpoint)){
		max_count <- max(dc)
		de_cutpoint <- max_count/log_base
	}
	n_over <- sum(dc > de_cutpoint)
	if (n_over < 10){
		stop("Less than 10 droplets are over the DE cutpoint. Set to a lower value.\n")
	}

	if (verbose) cat(paste0("Initialized DE cutpoint at ", as.character(de_cutpoint), " counts\n"))

	dc <- dc[order(dc, decreasing=TRUE)]

	x <- set_de_groups(x, low=names(dc)[dc <= de_cutpoint], high=names(dc)[(dc > de_cutpoint)])

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
get_log2fc <- function(x, 
	cpm_thresh=3, 
	log2fc_thresh=0.25, 
	verbose=FALSE){
	
	if (length(x@high_droplets) == 0){
		stop("Initialize DE cutpoint first")
	}

	low_sum <- Matrix::rowSums(x@counts[, x@low_droplets])
	high_sum <- Matrix::rowSums(x@counts[, x@high_droplets])

	low_sum <- 1e6*low_sum/sum(low_sum)
	high_sum <- 1e6*high_sum/sum(high_sum)

	genes <- rownames(x@counts)
	genes_expr <- genes[ low_sum > cpm_thresh & high_sum > cpm_thresh ]
	low_sum <- low_sum[genes_expr]; high_sum <- high_sum[genes_expr]
	log2fc <- log2(low_sum) - log2(high_sum)
	names(log2fc) <- genes_expr

	deg_low <- genes_expr[ log2fc > log2fc_thresh ]
	deg_high <- genes_expr[ log2fc < -log2fc_thresh ]
	deg <- c(deg_low, deg_high)

	# Save in DE
	de_obj <- DE(low_means=low_sum[deg],
				 high_means=high_sum[deg],
				 deg=deg, 
				 deg_low=deg_low, 
				 deg_high=deg_high,
				 log2fc=log2fc[deg])
	x@de <- de_obj
	return(x)
}
