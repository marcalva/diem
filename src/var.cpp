
// [[Rcpp::depends(RcppEigen)]]
//
// #include <Rcpp.h>
#include <RcppEigen.h>
#include <math.h>
#include <stdio.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector fast_varCPP(Eigen::SparseMatrix<double> x, 
        NumericVector mu){
    x = x.transpose();
    int n_c = x.rows();
    int n_g = x.cols();
    NumericVector var = no_init(n_g);
    for (int i = 0; i < n_g; i++){
        double v = 0;
        int n_zeroes = n_c;
        for (Eigen::SparseMatrix<double>::InnerIterator it(x,i); it; ++it) {
            n_zeroes -= 1;
            v += pow(it.value() - mu[i], 2);
        }
        v += pow(mu[i], 2) * n_zeroes;
        var[i] = v / (n_c - 1);
    }
    return(var);
}
