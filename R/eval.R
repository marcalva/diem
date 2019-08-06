
#' Get pi projection of normalized data
#'
#' This function calculates pi_low and pi_high, defined as the inner products of 
#' normalized gene expression and log fold changes between the background and 
#' non-background means, respectively.
#'
#' @param x SCE. SCE object.
#'
#' @return An matrix of pi values
#' @importFrom Matrix t
#' @export
get_pi <- function(x){

    mu <- x@diem@emo[[length(x@diem@emo)]]$Mu
    m_diff <- log2(mu[,2]) - log2(mu[,1]) # background - signal
    names(m_diff) <- rownames(mu)

    genes_b <- names(m_diff)[m_diff > 0]
    genes_s <- names(m_diff)[m_diff < 0]

    genes <- rownames(x@gene_info)[x@gene_info$exprsd]
    droplets <- names(x@labels)[x@labels != 2]

    counts <- x@counts[genes, droplets]
    counts_norm <- norm_counts(counts)

    pi_l <- Matrix::t(counts_norm[genes_b,]) %*% m_diff[genes_b]
    pi_h <- Matrix::t(counts_norm[genes_s,]) %*% -m_diff[genes_s]

    pi_df <- as.data.frame(cbind(as.numeric(pi_l), as.numeric(pi_h)))
    rownames(pi_df) <- colnames(counts_norm)
    colnames(pi_df) <- c("pi_l", "pi_h")

    x@diem@pi <- as.matrix(pi_df)

    return(x)
}

#' KL divergence for probabilities
#' @export
kldiverg <- function(p,q){
    kld <- -sum(p*log(q/p))
    return(kld)
}

#' KL divergence between low count with debris and signal
#' @export
kl_background <- function(x){
    db_mult <- x@diem@emo[[length(x@diem@emo)]]$Mu[,2]
    sg_mult <- x@diem@emo[[length(x@diem@emo)]]$Mu[,1]
    low_drop_mult <- Matrix::rowSums(x@counts[names(db_mult),x@labels == 2])
    low_drop_mult <- low_drop_mult[names(db_mult)]

    keep <- (db_mult > 0) & (sg_mult > 0) & (low_drop_mult > 0)
    db_mult <- db_mult[keep]
    sg_mult <- sg_mult[keep]
    low_drop_mult <- low_drop_mult[keep]

    db_mult <- db_mult/sum(db_mult)
    low_drop_mult <- low_drop_mult/sum(low_drop_mult)

    kld_bg_db <- kldiverg(low_drop_mult, db_mult)
    kld_bg_sg <- kldiverg(low_drop_mult, sg_mult)

    return(c("Debris"=kld_bg_db, "Signal"=kld_bg_sg))
}

#' Get Mu for debris, signal, and low-count fixed droplets
#' @export
get_mu_all <- function(x){
    Mu <- x@diem@emo[[length(x@diem@emo)]]$Mu
    low_drop_mult <- Matrix::rowSums(x@counts[names(db_mult),x@labels == 2])

    keep <- intersect(rownames(Mu), names(low_drop_mult))

    Mu <- Mu[keep,]
    low_drop_mult <- low_drop_mult[keep]

    mu_all <- cbind(Mu, low_drop_mult)
    colnames(mu_all) <- c("Signal", "Debris", "LowCount")

    mu_all <- sweep(mu_all, 2, colSums(mu_all), "/")
    mu_all <- log1p(1e4*mu_all)
    return(mu_all)
}

