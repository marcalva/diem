
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

