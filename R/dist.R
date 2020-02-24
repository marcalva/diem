
#' Get cluster distances to background disrtibution
#' 
#' Get the likelihood-based distance of the clusters to the background 
#' distribution. The distances take values between 0 and 1, where a 
#' larger value for a cluster indicates it is further away and thus more 
#' likely to reprsent a cell type. The distance values can plotted against 
#' average number of UMIs in a cluster using \code{\link{plot_dist}} to 
#' help decide which threshold to use. It is recommended to use a threshold 
#' so that improperly initialized clusters can be removed. The distances 
#' can be retrieved by calling the \code{\link{distances}} function.
#'
#' @param x An SCE object.
#' @param k_init An integer specifying the k_init run(s) to remove clusters 
#'  from. The default is to remove from all k_init values.
#' @param verbose Verbosity.
#'
#' @return An SCE object
#'
#' @export
get_dist <- function(x, k_init = NULL, verbose = TRUE){
    if (verbose) message("checking distances to background distribution...")

    if (length(x@kruns) == 0) stop("initialize parameters before running get_dist")

    k_init <- check_k_init(x, k_init)

    for (k in k_init){
        kc <- as.character(k)

        if (ncol(x@kruns[[kc]]$params$Alpha) == 1){
            if (verbose){
                warning("warning: only the debris distribution is present, no distances calculated")
            }
            return(x)
        }

        Z <- get_z(x@kruns[[kc]]$llk, x@kruns[[kc]]$params$Pi)

        d <- c(0)
        for (k in 2:ncol(x@kruns[[kc]]$params$Alpha)){
            diffs <- Z[,k] * (x@kruns[[kc]]$llk[,k] - x@kruns[[kc]]$llk[,1])
            wsum <- sum(Z[,k])
            if (wsum == 0){
                d[k] <- 0
            } else {
                d[k] <- sum(diffs) / wsum
            }
        }
        d[d < 0] <- 0
        d[is.na(d)] <- 0
        ds <- c()
        if (max(d) == 0){
            warning("warning: all distances are 0")
            ds <- d
        } else {
            ds <- d / max(d)
        }
        ds[1] <- NA
        x@kruns[[kc]]$Dist <- ds
    }
    return(x)

}

#' Remove clusters close to the background distribution
#'
#' Remove clusters with a distance less than \code{fltr} to the background. 
#' Distances are obtained using \code{\link{get_dist}} from either the 
#' initialized centers or after running EM.
#' The filter value is set to 0.1 by default, but it is recommended to 
#' set this appropriately to the your data set. A helpful way to discrminate 
#' background clusters from cell types is to look at the distribution. The 
#' distance values can be obtained with the function \code{\link{distances}}.
#' It is also helpful to plot them against their 
#' average UMI size using \code{\link{plot_dist}}. Clusters with
#' very low read counts are likely to from droplets with ambient RNA and the 
#' filter threshld should be set to exclude these while including the 
#' cell type clusters.
#' 
#' @param x An SCE object.
#' @param fltr The filter threshold between 0 and 1 
#'  that controls the minimum distance to 
#'  the background distribution that a cluster can have. Remove those 
#'  centers with a distance less than this value.
#' @param k_init An integer specifying the k_init run(s) to remove clusters 
#'  from. The default is to remove from all k_init values.
#' @param verbose Verbosity
#'
#' @return An SCE object.
#'
#' @export
rm_close <- function(x, 
                     fltr = 0.1, 
                     k_init = NULL, 
                     verbose = TRUE){ 

    if (is.na(fltr)){
        stop("fltr value is set to NA, change to a value between 0 and 1")
    }

    if (length(x@kruns) == 0) 
        stop("initialize parameters and get distances before running rm_close")

    k_init <- check_k_init(x, k_init)

    for (k in k_init){
        if (verbose){
            message("removing clusters for k_init = ", k_init, 
                    " using filter = ", fltr)
        }
        kc <- as.character(k)
        i <- length(x@kruns[[kc]])

        if (! "Dist" %in% names(x@kruns[[kc]]))
            stop("calculate distances before running rm_close")

        ds <- x@kruns[[kc]]$Dist[-1]
        skeep <- c(TRUE, (ds >= fltr))
        if (sum(!skeep) > 0){
            x@kruns[[kc]]$params$Alpha <- x@kruns[[kc]]$params$Alpha[,skeep,drop=FALSE]
            x@kruns[[kc]]$params$Pi <- x@kruns[[kc]]$params$Pi[skeep]
            x@kruns[[kc]]$params$Pi <- x@kruns[[kc]]$params$Pi / sum(x@kruns[[kc]]$params$Pi)
            x@kruns[[kc]]$llk <- x@kruns[[kc]]$llk[,skeep,drop=FALSE]
            ds <- x@kruns[[kc]]$Dist[skeep]
            x@kruns[[kc]]$Dist <- ds
        }
        x@kruns[[kc]]$removed <- sum(!skeep)

        if (verbose){
            torm <- x@kruns[[kc]]$removed
            message("removed ", torm, 
                    ngettext(torm, " cluster", " clusters"),
                    " with a distance less than ", 
                    fltr, " to the background distribution", 
                    " for k_init = ", k_init)
            nclusters <- length(x@kruns[[kc]]$params$Pi) - 1
            message("using 1 background and ", 
                    nclusters, " cell type clusters", 
                    " for k_init = ", k_init)
        }
    }

    return(x)
}

