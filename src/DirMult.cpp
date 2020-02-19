
// [[Rcpp::depends(RcppEigen)]]
//
// #include <Rcpp.h>
#include <RcppEigen.h>
#include <math.h>
#include <stdio.h>
using namespace Rcpp;
using namespace std;

// Compute log density of Dirichlet-Multinomial PMF with sparse matrix
// x is a gene by droplet sparse matrix
// [[Rcpp::export]]
NumericVector LlkDirMultSparse(Eigen::SparseMatrix<double> x, 
        NumericVector sizes,
        NumericMatrix alpha, 
        bool debug = false){
    int n_g = x.rows();
    int n_c = x.cols();
    if ( sizes.length() != n_c ){
        stop("Length of sizes must be the same as the number of columns in x");
    }
    if ( alpha.rows() != n_g ){
        stop("The number of rows in alpha must be the same as the number of rows in x");
    }

    int K = alpha.ncol();
    NumericMatrix llks(n_c, K);

    for (int k = 0; k < K; k++){
        NumericVector alpha_k = alpha(_,k);
        double as = sum(alpha_k);

        double l = 0;
        // NumericVector llks(n_c);
        for (int i = 0; i < n_c; ++i){
            double s = sizes(i);
            l = lgamma(sizes(i) + 1) + 
                lgamma(as) - 
                lgamma(sizes(i) + as);
            for (Eigen::SparseMatrix<double>::InnerIterator it(x,i); it; ++it){
                l += lgamma(it.value() + alpha_k(it.index()) ) - 
                    lgamma(it.value() + 1) - 
                    lgamma( alpha_k(it.index()) );
                if (isnan(l)) {
                    Rcout << "alpha_g = " << alpha_k(it.index()) << "\n";
                    Rcout << "Likelihood value is NA for cell " << i << 
                        " after adding gene " << it.index() << "\n";
                    stop("Exiting");
                }
            }
            if (isnan(l)) {
                stop("nan value encountered. There are likely non-positive values in x, sizes, or alpha.");
            }
            llks(i,k) = l;
        }
    }
    return(llks);
}

