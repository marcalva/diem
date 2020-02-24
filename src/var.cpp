
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <math.h>
#include <stdio.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericVector fast_varCPP(Eigen::SparseMatrix<double> x, 
        NumericVector mu, 
        int threads = 1, 
        bool display_progress = false){
#ifdef _OPENMP
    if ( threads > 0 ){
        omp_set_num_threads( threads );
    }
    if (display_progress){
        REprintf("Number of threads=%i\n", threads);
    }
#endif

    x = x.transpose();
    int n_c = x.rows();
    int n_g = x.cols();
    NumericVector var(n_g);

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
    for (int i = 0; i < n_g; i++){
        var[i] = 0;
        int n_zeroes = n_c;
        for (Eigen::SparseMatrix<double>::InnerIterator it(x,i); it; ++it) {
            n_zeroes -= 1;
            var[i] += pow(it.value() - mu[i], 2);
        }
        var[i] += pow(mu[i], 2) * n_zeroes;
        var[i] = var[i] / (n_c - 1);
    }
    return(var);
}
