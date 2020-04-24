
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
NumericMatrix fast_correct(Eigen::SparseMatrix<double> x, 
        NumericMatrix means,
        NumericMatrix props, 
        bool round_count = true,
        bool display_progress = false){

    int n_clust = means.cols();
    if (props.rows() != n_clust){
        stop("Number of rows in props must be same as number of columns in means");
    }
    int n_c = x.cols();
    int n_g = x.rows();
    NumericVector total_counts(n_c);
    for (int i = 0; i < n_c; ++i){
        total_counts(i)  = 0;
        for (Eigen::SparseMatrix<double>::InnerIterator it(x,i); it; ++it) {
            total_counts(i) += it.value();
        }
    }

    NumericMatrix crrctd(n_g, n_c);

    for (int i = 0; i < n_c; i++){
        for (Eigen::SparseMatrix<double>::InnerIterator it(x,i); it; ++it) {
            double subt = 0;
            for (int k = 0; k < n_clust; k++){
                subt += (means(it.index(),k) * props(k,i));
            }
            subt *= total_counts(i);
            if (round_count){
                crrctd(it.index(), i) = it.value() - round(subt);
            } else {
                crrctd(it.index(), i) = it.value() - subt;
            }
            if (crrctd(it.index(), i) < 0)
                crrctd(it.index(), i) = 0;
        }
    }
    return(crrctd);
}

