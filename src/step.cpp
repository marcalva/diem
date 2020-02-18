// [[Rcpp::depends(RcppEigen)]]
//
// #include <Rcpp.h>
#include <RcppEigen.h>
#include <math.h>
#include <vector>
#include <stdio.h>
#include <algorithm>
using namespace Rcpp;

// Compute LOO step for Dirichlet-Multinomial
// [[Rcpp::export]]
NumericVector compute_LOO_step_all(Eigen::SparseMatrix<double> x, 
        NumericVector sizes,
        NumericVector weights,
        NumericVector alpha, 
        double eps = 1e-4, 
        int max_loo = 500, 
        double psc = 1e-10, 
        bool debug = false){
    int n_c = x.rows();
    int n_g = x.cols();
    double as = sum(alpha);
    NumericVector alpha_new(n_g);
    NumericVector alpha_old = clone(alpha);
    // std::copy( alpha.begin(), alpha.end(), alpha_old.begin() ) ;


    double delt = eps + 1;
    int iter = 1;
    // calculate weighted gene sums
    while (delt > eps && iter <= max_loo){
        double as = sum(alpha_old);

        double denom = 0;
        for (int i = 0; i < n_c; ++i){
            denom += (weights(i) * sizes(i)) / (sizes(i) - 1 + as);
        }
        if (isnan(denom)){
            Rcout << "Denominator " << denom << "\n";
            stop("NA values encountered. An alpha value is likely 0.");
        }

        for (int k = 0; k < n_g; ++k){
            double numer = 0;
            for (Eigen::SparseMatrix<double>::InnerIterator it(x,k); it; ++it) {
                double xik = it.value();
                int i = it.index();
                numer += (weights(i) * xik) / (xik - 1 + alpha_old(k));
            }
            alpha_new(k) = alpha_old(k) * (numer / denom);
            if (isnan(alpha_new(k))){
                Rcout << "alpha_k " << alpha_old(k) << "\n";
                Rcout << "value " << (numer / denom) << "\n";
                Rcout << "Numerator " << numer << "\n";
                Rcout << "Denominator " << denom << "\n";
                stop("NA values encountered. An alpha_old value is likely 0.");
            }
            if (alpha_new(k) < 0) {
                alpha_new(k) = psc;
            } else {
                alpha_new(k) += psc;
            }
        }
        //NumericVector dv(n_g);
        //for (int g = 0; g < n_g; g++){
        //    dv(g) = abs(alpha_new(g) - alpha_old(g));
        //}
        //delt = sum(dv) / sum(alpha_old);
        
        delt = sum(abs(alpha_new - alpha_old)) / sum(alpha_old);
        iter += 1;
        alpha_old = clone(alpha_new);
    }
    return(alpha_new);
}

// Compute LOO step for Dirichlet-Multinomial
// [[Rcpp::export]]
NumericVector compute_LOO_step(Eigen::SparseMatrix<double> x, 
        NumericVector sizes,
        NumericVector weights,
        NumericVector alpha, 
        double psc = 1e-10, 
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

