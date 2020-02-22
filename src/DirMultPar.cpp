
// [[Rcpp::plugins(openmp)]]

#include <Rcpp.h>
#include <RcppEigen.h>
#include <math.h>
#include <stdio.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <progress.hpp>
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppEigen)]]

// Compute log density of Dirichlet-Multinomial PMF with sparse matrix
// x is a gene by droplet sparse matrix
// [[Rcpp::export]]
NumericVector LlkDirMultSparsePar(Eigen::SparseMatrix<double> x, 
        NumericVector sizes,
        NumericMatrix alpha, 
        int threads = 1, 
        bool display_progress = true,
        bool debug = false){
    int n_g = x.rows();
    int n_c = x.cols();
    if ( sizes.length() != n_c ){
        stop("Length of sizes must be the same as the number of columns in x");
    }
    if ( alpha.rows() != n_g ){
        stop("The number of rows in alpha must be the same as the number of rows in x");
    }

#ifdef _OPENMP
    if ( threads > 0 ){
        int mt = omp_get_max_threads();
        if (threads > mt){
            threads = mt;
        }
        omp_set_num_threads( threads );
    }
    if (display_progress){
        REprintf("Number of threads=%i\n", threads);
    }
#endif

    int K = alpha.ncol();
    double* llk = new double[n_c * K]();
    double* asum = new double[K]();
    for (int i = 0; i < K; ++i){
        asum[i] = sum(alpha(_,i));
    }

    Progress p(n_c * K, display_progress);


#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static) collapse(2) 
#endif
    for (int i = 0; i < n_c; ++i){
        for (int k = 0; k < K; k++){
            if ( ! Progress::check_abort() )
                p.increment();
            llk[i + (k * n_c)] = lgamma(sizes(i) + 1) + 
                lgamma(asum[k]) - 
                lgamma(sizes(i) + asum[k]);
            for (Eigen::SparseMatrix<double>::InnerIterator it(x,i); it; ++it){
                llk[i + (k * n_c)] += lgamma(it.value() + alpha(it.index(), k) ) - 
                    lgamma(it.value() + 1) - 
                    lgamma( alpha(it.index(), k) );
                if (isnan(llk[i + (k * n_c)])) {
                    Rcout << "alpha_g = " << alpha(it.index(), k) << "\n";
                    Rcout << "Likelihood value is NA for cell " << i << 
                        " after adding gene " << it.index() << "\n";
                    stop("Exiting");
                }
            }
            if (isnan(llk[i + (k * n_c)])) {
                stop("nan value encountered. There are likely non-positive values in x, sizes, or alpha.");
            }
        }
    }

    NumericVector llk_rv(llk, llk + sizeof(llk)/sizeof(*llk));
    llk_rv.attr("dim") = Dimension(n_c, K);
    delete [] llk;
    delete [] asum;
    return(llk_rv);
}

