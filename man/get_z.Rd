% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/em.R
\name{get_z}
\alias{get_z}
\title{Get posterior probabilities from log likelihoods and mixing coefficients}
\usage{
get_z(llk, Pi)
}
\arguments{
\item{llk}{A droplet by cluster matrix that gives the log liklihood of 
each droplet against each cluster.}

\item{Pi}{A cluster-length numeric vector of mixing coefficients.}
}
\value{
A droplet by cluster matrix of posterior probabilities
}
\description{
Get posterior probabilities from log likelihoods and mixing coefficients
}
