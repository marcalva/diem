
// [[Rcpp::depends(RcppEigen)]]
//
// #include <Rcpp.h>
#include <RcppEigen.h>
#include <math.h>
#include <stdio.h>
using namespace Rcpp;

// Compute LOO step for Dirichlet-Multinomial
// [[Rcpp::export]]
NumericVector compute_LOO_step(Eigen::SparseMatrix<double> x, 
        NumericVector sizes,
        NumericVector weights,
        NumericVector alpha, 
        int tol = 100, 
        bool debug = false){
    int n_c = x.rows();
    int n_g = x.cols();
    double as = sum(alpha);
    NumericVector alpha_new(n_g);

    double denom = 0;
    for (int i = 0; i < n_c; ++i){
        denom += (weights[i] * sizes[i]) / (sizes[i] - 1 + as);
    }

    if (isnan(denom)){
        Rcout << "Denominator " << denom << "\n";
        stop("NA values encountered. An alpha value is likely 0.");
    }

    // calculate weighted gene sums
    for (int k = 0; k < n_g; ++k){
        double numer = 0;
        for (Eigen::SparseMatrix<double>::InnerIterator it(x,k); it; ++it) {
            double xik = it.value();
            int i = it.index();
            numer += (weights[i] * xik) / (xik - 1 + alpha[k]);
        }
        alpha_new[k] = alpha[k] * (numer / denom);
        if (isnan(alpha_new[k])){
            Rcout << "alpha_k " << alpha[k] << "\n";
            Rcout << "value " << (numer / denom) << "\n";
            Rcout << "Numerator " << numer << "\n";
            Rcout << "Denominator " << denom << "\n";
            stop("NA values encountered. An alpha value is likely 0.");
        }
    }
    return(alpha_new);
}

