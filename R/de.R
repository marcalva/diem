
##' Get variable genes
##'
##' Finds top \code{nf} variable genes by running a loess regression 
##' on log(variance) and log(mean). Residuals are used to rank 
##' features, and the top \code{nf} are output.
##'
##' @param x SCE. SCE object.
##' @param min_mean Numeric. Minimum mean of counts to consider.
##' @param nf Numeric. Output top nf variable features.
##' @param lss Numeric. The span paramater for loess function.
##' @param verbose Boolean.
##'
##' @return SCE object with \code{vg} and \code{vg_info} slots filled for variable genes
##' and variable gene info. \code{vg_info} is a data frame with log mean, log var, 
##' fitted log var values, and residuals. \code{vg} is a vector of genes
##' @export
#get_var_genes <- function(x, min_mean = 0.0001, nf = 2000, lss=0.3, verbose=FALSE){
#	if (verbose) cat("Getting variable genes\n")
#	if (min_mean <= 0){
#		keep <- x@gene_info[,"mean"] > 0
#	} else {
#		keep <- x@gene_info[,"mean"] > min_mean
#	}
#	if (table(keep)["TRUE"] < nf){
#		stop("Error: Specifying more variable features than there are features passing min_mean\n")
#	}
#	gene_means <- x@gene_info[,"mean"][keep]
#	gene_names <- rownames(x@gene_info)[keep]
#
#	log_mean <- log10(gene_means + 1)
#	log_var <- log10(fast_varCPP(x@counts[keep,], gene_means) + 1)
#	fit <- loess(log_var ~ log_mean, span=lss)
#	rsd <- log_var - fit$fitted
#	topi <- order(rsd, decreasing=TRUE)[1:nf]
#	vg <- gene_names[topi]
#	df <- data.frame(mean=log_mean, var=log_var,
#					 fit=fit$fitted, rsd=rsd)
#	rownames(df) <- gene_names
#	x@vg <- vg
#	x@vg_info <- df
#	if (verbose) cat("Found variable genes\n")
#	return(x)
# }

#' Run t test given parameters of the mean, variance, and size.
#'
#' This t test function rapidly calculates p value and log fold change
#' given the parameters.
#'
#' @param x_mean Numeric. Mean of group x.
#' @param x_var Numeric. Variance of group x.
#' @param x_n Numeric. Size of group x.
#' @param y_mean Numeric. Mean of group y.
#' @param y_var Numeric. Variance of group y.
#' @param y_n Numeric. Size of group y.
#'
#' @return list with 2 components:
#' \describe{
#'   \item{p}{p-value}
#'   \item{log2fc}{Log base 2 fold change (x over y)}
#' }
#' 
#' @export
t_test_sparse <- function(x_mean, x_var, x_n, y_mean, y_var, y_n){
	tv <- (x_mean - y_mean) / sqrt( (x_var/x_n) + (y_var/y_n) )
	dof <- ( (x_var/x_n) + (y_var/y_n) )^2 / ( (x_var/x_n)^2/(x_n-1) + (y_var/y_n)^2/(y_n-1) )
	pv_lt <- pt(tv, df=dof, lower.tail=TRUE)
	pv_ut <- pt(tv, df=dof, lower.tail=FALSE)
	pv <- min(pv_lt, pv_ut) * 2
	log2fc <- log2(x_mean/y_mean); names(log2fc) <- ""
	return(list(p=pv, log2fc=log2fc))
}

#' Set count cutpoint for DE testing between low and high count droplets
#'
#' In order to find genes enriched in the background or the nucleus, we need to 
#' determine how to group droplets for DE. This function finds a proper cut point 
#' to decide which droplets are assigned to either the high or low group. This number 
#' is determined by taking the log10 of the maximum droplet count, subtracting 1, 
#' and calculating 10 to the power of this number (10^(log10(x)-1) where x is the max 
#' count. Alternatively, the you can set the cutpoint threshold manually.
#'
#' @param x SCE. SCE object.
#' @param de_cutpoint Numeric. Manually set the DE cutpoint.
#' @param verbose Boolean. Verbosity.
#'
#' @return SCE object.
#' @export
set_de_cutpoint <- function(x, 
							de_cutpoint=NULL, 
							verbose=FALSE){
	if (verbose) cat("Setting DE cutpoint\n")
	if (!is.null(de_cutpoint)){
		x@de_cutpoint <- de_cutpoint
	} else {
		max_count <- max(x@dropl_info$total_counts)
		max_count_l <- log10(max_count)
		de_cutpoint_l <- max_count_l-1
		x@de_cutpoint <- 10^de_cutpoint_l
	}
	n_over <- sum(x@dropl_info$total_counts > x@de_cutpoint)
	if (n_over < 10){
		stop("Less than 10 droplets are over the DE cutpoint. Set to a lower value.\n")
	}
	if (verbose) cat(paste0("Set DE cutpoint at ", as.character(x@de_cutpoint), " counts\n"))
	return(x)
}


#' Get DE genes between high and low read count droplets using a t-test
#'
#' Find differentially expressed genes between low and high read count droplets. Groups droplets into
#' low or high if they have less than or equal to or greater than \code{de_cutpoint}, respectively. Must run
#' set \code{de_cutpoint} by running \code{\link{set_de_cutpoint}} before.A t-test is run between the 2 groups 
#' for each gene, calculating log fold change and p-value. Statistically significant genes are determined by adjusting 
#' the p-value using the \code{de_correct}. Can be any method used in \code{\link[stats]{p.adjust}}, typically either 
#' "bonferroni" (default) or "fdr."
#'
#' @param x SCE. SCE object.
#' @param log2fc_thresh Numeric. Test only genes that have this log2 fold-change threshold between low and high
#' droplets will be considered as differentially expressed. Fold change is for both comparison low vs high and high vs low.
#' @param de_p_thresh Numeric. Adjusted p-value threshold.
#' @param de_correct Character. p-value correction method.
#' @param verbose Boolean.
#'
#' @return SCE object.
#' @importFrom Matrix Diagonal colSums rowMeans
#' @export
get_de <- function(x, 
				   log2fc_thresh=0.5, 
				   de_p_thresh=0.05, 
				   de_correct="bonferroni", 
				   verbose=FALSE){
	if (verbose) cat("Finding differentially expressed genes\n")

	if (is.null(x@de_cutpoint)){
		stop("Can't find DE genes. Set DE cutpoint first.")
	}
	dc <- x@dropl_info[,"total_counts"]
	low_droplets <- (dc <= x@de_cutpoint)
	high_droplets <- (dc > x@de_cutpoint)

	low_expr <- x@counts[, low_droplets]
	high_expr <- x@counts[, high_droplets]

	# Normalize count size
	low_prop <- divide_by_colsum(low_expr)
	high_prop <- divide_by_colsum(high_expr)

	# low_prop <- log1p(low_prop)
	# high_prop <- log1p(high_prop)

	# Calculate mean, variance, and n
	low_means <- Matrix::rowMeans(low_prop)
	high_means <- Matrix::rowMeans(high_prop)
	low_vars <- fast_varCPP(low_prop, low_means); names(low_vars) <- names(low_means)
	high_vars <- fast_varCPP(high_prop, high_means); names(high_vars) <- names(high_means)
	low_n <- ncol(low_prop)
	high_n <- ncol(high_prop)

	# Run DE
	genes_use <- rownames(high_expr)[ low_means > 0 & high_means > 0 ]
	t_results <- sapply(genes_use, function(i){
						tr <- t_test_sparse(low_means[i], low_vars[i], low_n, high_means[i], high_vars[i], high_n)
						return(c("p"=tr$p, "log2fc"=tr$log2fc))
					 })
	t_results <- as.data.frame(t(t_results))
	t_results[,"p_adj"] <- p.adjust(t_results[,"p"], method=de_correct)

	deg <- rownames(t_results)[ (abs(t_results[,"log2fc"]) > log2fc_thresh) & (t_results[,"p_adj"] < de_p_thresh) ]
	deg_low <- rownames(t_results)[ (t_results[,"log2fc"] > log2fc_thresh) & (t_results[,"p_adj"] < de_p_thresh) ]
	deg_high <- rownames(t_results)[ (t_results[,"log2fc"] < log2fc_thresh) & (t_results[,"p_adj"] < de_p_thresh) ]

	# Save in DE
	de_obj <- DE(low_means=low_means,
				 high_means=high_means,
				 low_vars=low_vars,
				 high_vars=high_vars,
				 low_n=low_n,
				 high_n=high_n,
				 deg=deg, 
				 deg_low=deg_low, 
				 deg_high=deg_high,
				 table=t_results)
	x@de <- de_obj

	# Check if de genes can separate background from top expressed droplets
	if (length(x@de@deg) < 10) stop("Less than 10 genes found differentially expressed between low and high read droplets")
	if ( length(x@de@deg) < 50 ){
		cat("Warning: Less than 50 genes found DE between low and high read droplets")
	} else {
		if (verbose) cat(paste0("Found ", as.character(length(x@de@deg)), " DE genes\n"))
	}

	return(x)
}
