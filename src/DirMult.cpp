
// [[Rcpp::depends(RcppEigen)]]
//
// #include <Rcpp.h>
#include <RcppEigen.h>
#include <math.h>
#include <stdio.h>
#include "digamma.h"
#include "trigamma.h"
using namespace Rcpp;

// Compute log density of Dirichlet-Multinomial PMF with sparse matrix
// x is a gene by droplet sparse matrix
// [[Rcpp::export]]
NumericVector LlkDirMultSparse(Eigen::SparseMatrix<double> x, 
        NumericVector sizes,
        NumericVector alpha, 
        bool debug = false){
    //x = x.transpose(); // x is now gene by droplet
    int n_g = x.rows();
    int n_c = x.cols();
    // int ifault;
    double as = sum(alpha);

    NumericVector llks(n_c);
    for (int i = 0; i < n_c; ++i){
        double s = sizes[i];
        double l = lgamma(sizes[i] + 1) + 
            lgamma(as) - 
            lgamma(sizes[i] + as);
        for (Eigen::SparseMatrix<double>::InnerIterator it(x,i); it; ++it){
            l += lgamma(it.value() + alpha[it.index()]) - 
                lgamma(it.value() + 1) - 
                lgamma(alpha[it.index()]);
        }
        llks[i] = l;
    }
    return(llks);
}

