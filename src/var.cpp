
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



//
// [[Rcpp::export]]
NumericVector fast_wvarCPP(Eigen::SparseMatrix<double> x, 
        NumericVector mu, 
        NumericVector weights,
        int threads = 1, 
        bool display_progress = false){
#ifdef _OPENMP
    if ( threads > 0 ){
        omp_set_num_threads( threads );
    }
#endif

    int n_c = x.rows();
    int n_g = x.cols();
    NumericVector var(n_g);

    if (weights.length() != n_c)
        stop("Weights must have the same length as the number of droplets");

    double wsum = 0;
    for (int i = 0; i < n_c; ++i){
        wsum += weights[i];
    }
    
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
    for (int i = 0; i < n_g; ++i){
        var[i] = 0;
        double n_zeroes = wsum;
        for (Eigen::SparseMatrix<double>::InnerIterator it(x,i); it; ++it) {
            n_zeroes -= weights[it.index()];
            var[i] += weights[it.index()] * pow(it.value() - mu[i], 2);
        }
        if (n_zeroes < 0){
            n_zeroes = 0;
        }
        var[i] += pow(mu[i], 2) * n_zeroes;
        var[i] = var[i] / wsum;
    }
    return(var);
}

// [[Rcpp::export]]
NumericVector fast_wmeanCPP(Eigen::SparseMatrix<double> x, 
        NumericVector weights,
        int threads = 1, 
        bool display_progress = false){
#ifdef _OPENMP
    if ( threads > 0 ){
        omp_set_num_threads( threads );
    }
#endif

    int n_c = x.rows();
    int n_g = x.cols();
    NumericVector means(n_g);


    if (weights.length() != n_c)
        stop("Weights must have the same length as the number of droplets");

    double wsum = 0;
    for (int i = 0; i < n_c; ++i){
        wsum += weights[i];
    }
    
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
    for (int i = 0; i < n_g; ++i){
        means[i] = 0;
        for (Eigen::SparseMatrix<double>::InnerIterator it(x,i); it; ++it) {
            means[i] += weights[it.index()] * it.value();
        }
        means[i] = means[i] / wsum;
    }
    return(means);
}

