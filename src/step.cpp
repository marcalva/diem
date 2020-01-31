
// [[Rcpp::depends(RcppEigen)]]
//
// #include <Rcpp.h>
#include <RcppEigen.h>
#include <math.h>
#include <stdio.h>
#include "digamma.h"
#include "trigamma.h"
using namespace Rcpp;

double round_to_digits(double value, int digits)
{
    if (value == 0.0) // otherwise it will return 'nan' due to the log10() of zero
        return 0.0;

    double factr = pow(10, digits);
    double r = round(factr * value) / factr;
    return r;
}

// Compute step for Dirichlet-Multinomial
// [[Rcpp::export]]
NumericVector compute_step(Eigen::SparseMatrix<double> x, 
        NumericVector sizes,
        NumericVector weights,
        NumericVector alpha, 
        int tol = 100, 
        bool debug = false){
    // x = x.transpose(); // x is now droplet by gene
    int n_c = x.rows();
    int n_g = x.cols();
    int ifault;
    double as = sum(alpha);

    // Round weights to tolerance
    for (int i = 0; i < n_c; ++i){
        weights[i] = round_to_digits(weights[i], tol);
    }

    double ws = sum(weights);

    // calculate weighted gene sums
    NumericVector wgs(n_g);
    for (int k = 0; k < n_g; ++k){
        double s = 0;
        for (Eigen::SparseMatrix<double>::InnerIterator it(x,k); it; ++it){
            s += weights[it.index()] * it.value();
        }
        wgs[k] = s;
    }

    
    NumericVector sa(n_c);
    for (int i = 0; i < n_c; ++i){
        sa[i] = sizes[i] + as;
    }
    NumericVector sad = Rcpp::digamma(sa);
    NumericVector sat = Rcpp::trigamma(sa);
    NumericVector sadw(n_c);
    NumericVector satw(n_c);
    for (int i = 0; i < n_c; ++i){
        sadw[i] = weights[i] * sad[i];
        satw[i] = weights[i] * sat[i];
    }
    double ssadw = sum(sadw);
    double ssatw = sum(satw);

    double asd = digamma(as, &ifault);
    double asdw = ws * asd;
    double ast = trigamma(as, &ifault);
    double astw = ws * ast;

    double z = astw - ssatw;

    /* ssadw is weighted sum of digamma(size_i + alpha_sum)
     * ssatw is weighted sum of trigamma(size_i + alpha_sum)
     * asdw is sum of weights times digamma(alpha_sum)
     * astw is sum of weights times trigamma(alpha_sum)
     */

    double t1 = asdw - ssadw;
    NumericVector q(n_g); 
    NumericVector g(n_g);
    double gqfrac = 0;
    double qinv_sum = 0;
    for (int k = 0; k < n_g; ++k){
        if (wgs[k] == 0) continue;
        g[k] = t1;
        int n_zero = n_c;
        for (Eigen::SparseMatrix<double>::InnerIterator it(x,k); it; ++it) {
            if (weights[it.index()] == 0) continue;
            n_zero -= 1;
            double xa = it.value() + alpha[k];
            // g[k] from this term is usually small
            double ds = digamma(xa, &ifault) - digamma(alpha[k], &ifault);
            g[k] += weights[it.index()] * ds;
            // Note very large alphas and small x values produce 0 for qk
            double ts = trigamma(xa, &ifault) - trigamma(alpha[k], &ifault);
            q[k] += weights[it.index()] * ts;
        }
        double gadd = g[k] / q[k];
        // q[k] = 0 tends to occur when alpha is much larger than x
        // g[k] tends to be positive
        if (q[k] == 0){
            // Rcout << "q[k] = 0, and alpha[k] = " << alpha[k] << "\n";
            if (g[k] > 0) gadd = pow(10, tol);
            else gadd = -pow(10, tol);
        }
        //Rcout << "q[k] " << q[k] << "\n";
        if (isinf(gadd)){
            Rcout << "Gene " << k + 1 << " is numerically unstable: " << g[k] << "/" << q[k] << "\n";
            continue;
        }
        gqfrac += gadd;
        qinv_sum += (1 / q[k]);
    }
    double b = gqfrac / ( (1/z) + qinv_sum );

    if (debug){
        Rcout << "z is " << z << "\n";
        Rcout << "Max g is " << max(g) << "\n";
        Rcout << "Min g is " << min(g) << "\n";
        Rcout << "Mean g is " << mean(g) << "\n";
        Rcout << "Median g is " << median(g) << "\n";

        Rcout << "Max q is " << max(q) << "\n";
        Rcout << "Min q is " << min(q) << "\n";
        Rcout << "Mean q is " << mean(q) << "\n";
        Rcout << "Median q is " << median(q) << "\n";
        Rcout << "The value of gqfrac is " << gqfrac << "\n";
        Rcout << "The value of qinv_sum is " << qinv_sum << "\n";
        Rcout << "The value of z is " << z << "\n";
        Rcout << "The value of b is " << b << "\n";
    }

    NumericVector hg(n_g);
    for (int k = 0; k < n_g; ++k){
        if (wgs[k] != 0 && q[k] == 0){
            //Rcout << "Possible mismatch for gene k " << k + 1 << "\n";
            //Rcout << "Alphak " <<  alpha[k] << "\n";
        }
        if (q[k] != 0) hg[k] = (g[k] - b) / q[k];
    }
    return(hg);
}

//
// Compute LOO step for Dirichlet-Multinomial
// [[Rcpp::export]]
NumericVector compute_LOO_step(Eigen::SparseMatrix<double> x, 
        NumericVector sizes,
        NumericVector weights,
        NumericVector alpha, 
        int tol = 100, 
        bool debug = false){
    // x = x.transpose(); // x is now droplet by gene
    int n_c = x.rows();
    int n_g = x.cols();
    double as = sum(alpha);
    NumericVector alpha_new(n_g);

    // Round weights to tolerance
    // for (int i = 0; i < n_c; ++i){
    //    weights[i] = round_to_digits(weights[i], tol);
    // }

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
            if (isnan(numer)){
                stop("NA values encountered. An alpha value is likely 0.");
            }
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

// Compute fixed point iteration for Dirichlet-Multinomial
// [[Rcpp::export]]
NumericVector compute_fp_step(Eigen::SparseMatrix<double> x, 
        NumericVector sizes,
        NumericVector weights,
        NumericVector alpha, 
        int tol = 100, 
        bool debug = false){
    // x = x.transpose(); // x is now droplet by gene
    int ifault;
    int n_c = x.rows();
    int n_g = x.cols();
    double as = sum(alpha);
    NumericVector alpha_new(n_g);

    // Round weights to tolerance
    for (int i = 0; i < n_c; ++i){
        weights[i] = round_to_digits(weights[i], tol);
    }
    double ws = sum(weights);

    double cs = ws * digamma(as, &ifault);
    double denom = -cs;
    for (int i = 0; i < n_c; ++i){
        double s = weights[i] * digamma(sizes[i] + as, &ifault);
        denom += s;
    }

    // calculate weighted gene sums
    for (int k = 0; k < n_g; ++k){
        double numer = 0;
        for (Eigen::SparseMatrix<double>::InnerIterator it(x,k); it; ++it) {
            double xik = it.value();
            int i = it.index();
            numer += weights[i] * (digamma(xik + alpha[k], &ifault) - digamma(alpha[k], &ifault));
        }
        if (isnan(numer)){
            Rcout << "NA numerator " << k << ";";
        }
        alpha_new[k] = alpha[k] * (numer / denom);
    }
    return(alpha_new);
}

