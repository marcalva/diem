
#' Check values of given k_init against kruns
#'
#' @param x An SCE object.
#' @param k_init Run EM on the \code{k_init} initialization(s). 
#'  If NULL (default), run on all \code{k_init} initializations.
#' @param return_all If there are more than one k_init values, return 
#'  all or throw an error.
#'
#' @return integer specifying the k_init value(s).
check_k_init <- function(x, k_init, return_all = TRUE){

    if (length(x@kruns) == 0){
        stop("no runs have been initialized")
    }

    if (is.null(k_init)){
        k_init <- names(x@kruns)
        if (!return_all){
            if (length(x@kruns) > 1){
                stop(sQuote("k_init"), " must be specified to a value from: ",
                     paste(names(x@kruns), collapse = " "))
            }
        }
    } else {
        if (!return_all){
            if (length(k_init) != 1){
                stop(sQuote("k_init"), " must be a single value")
            }
        }
    }
    if ( any(! k_init %in% names(x@kruns) ) ){
        stop(sQuote("k_init"), " must have an initialized k value taken from: ", 
             paste(names(x@kruns), collapse = " "))
    }
    return(k_init)
}

#' fraction of logs
#'
#' @param x numeric vector
fraction_log <- function(x){
    xinf <- is.infinite(x)
    if (any(xinf)){
        frac <- rep(0, length(x))
        frac[which(xinf)] <- 1 
    } else {
        x_c = x - max(x)
        x_c = exp(x_c)
        frac = x_c / sum(x_c);
    }
    return(frac);
}

#' sum of logs
#'
#' @param x numeric vector
sum_log <- function(x){
    max_x <- max(x)
    x_c = x - max_x
    x_sum <- log(sum(exp(x_c))) + max_x
    return(x_sum)
}

