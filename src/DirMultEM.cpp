
// [[Rcpp::plugins(openmp)]]

#include <Rcpp.h>
#include <RcppEigen.h>
#include <math.h>
#include <stdio.h>
#include <chrono>
#include <ctime>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <progress.hpp>
#include "utils.h"
#include "var.h"
#include "MultEM.h"
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppEigen)]]

typedef Eigen::SparseMatrix<double>::InnerIterator cit;

// Compute LOO step for Dirichlet-Multinomial
// [[Rcpp::export]]
Eigen::VectorXd max_loo(const Eigen::SparseMatrix<double> &x, 
        const Eigen::VectorXd &sizes,
        const Eigen::VectorXd &weights,
        const Eigen::VectorXd &alpha, 
        double eps = 1e-4, 
        int max_iter = 1e4, 
        double alpha_prior = 0, 
        double psc = 1e-10, 
        int threads = 1, 
        bool debug = false){

    int n_c = x.rows();
    int n_g = x.cols();
    double as = alpha.sum();
    Eigen::VectorXd alpha_new(n_g);
    Eigen::VectorXd alpha_old(alpha);

    if (weights.size() != n_c)
        stop("weights must be same length as number of rows in x");

    if (sizes.size() != n_c)
        stop("sizes must be same length as number of rows in x");

#ifdef _OPENMP
    if ( threads > 0 ){
        omp_set_num_threads( threads );
    }
#endif

    for (int k = 0; k < n_g; ++k){
        if (alpha_old(k) < 0)
            alpha_old(k) = 0;
        alpha_old(k) += psc;
    }
    double delt = eps + 1;
    int iter = 1;
    Eigen::VectorXd numer(n_g); 
    // calculate weighted gene sums
    while (delt > eps && iter <= max_iter){
        double as = alpha_old.sum();

        double denom = 0;
        for (int i = 0; i < n_c; ++i){
            denom += (weights(i) * sizes(i)) / (sizes(i) - 1 + as);
        }
        if (std::isnan(denom)){
            Rcout << "Denominator " << denom << "\n";
            stop("NA values encountered. An alpha value is likely 0.");
        }
        if (denom == 0){
            Rcout << "Denominator " << denom << " for " << n_c << " cells\n";
            stop("Denominator is 0");
        }

#ifdef _OPENMP
        #pragma omp parallel for num_threads(threads) schedule(dynamic)
#endif
        for (int k = 0; k < n_g; ++k){
            numer(k) = 0;
            for (Eigen::SparseMatrix<double>::InnerIterator it(x,k); it; ++it) {
                numer(k) += (weights(it.index()) * it.value()) / (it.value() - 1 + alpha_old(k));
            }
            alpha_new(k) = alpha_old(k) * (numer(k) / denom);
            if (std::isnan(alpha_new(k))){
                Rcout << n_c << " cells\n";
                Rcout << "alpha_k " << alpha_old(k) << "\n";
                Rcout << "value " << (numer(k) / denom) << "\n";
                Rcout << "Numerator " << numer(k) << "\n";
                Rcout << "Denominator " << denom << "\n";
                stop("NA values encountered. An alpha_old value is likely 0.");
            }
            if (alpha_new(k) < 0) {
                alpha_new(k) = psc + alpha_prior;
            } else {
                alpha_new(k) += psc + alpha_prior;
            }
        }
        
        delt = 0;
        for (int k = 0; k < n_g; ++k){
            delt += std::abs(alpha_new(k) - alpha_old(k)) / alpha_old(k);
            alpha_old(k) = alpha_new(k);
        }

        iter += 1;
    }
    return(alpha_new);
}

// Compute log density of Dirichlet-Multinomial PMF with sparse matrix
// x is a gene by droplet sparse matrix
// [[Rcpp::export]]
Eigen::MatrixXd LlkDirMultSparsePar(const Eigen::SparseMatrix<double> &x, 
        const Eigen::VectorXd &sizes,
        const Eigen::MatrixXd &alpha, 
        Eigen::VectorXi ix, 
        int ix_len, 
        int threads = 1, 
        bool display_progress = true,
        bool debug = false){
    int n_g = x.rows();
    int n_c = x.cols();
    if ( sizes.size() != n_c ){
        stop("Length of sizes must be the same as the number of columns in x");
    }
    if ( alpha.rows() != n_g ){
        stop("The number of rows in alpha must be the same as the number of rows in x");
    }

    // Check for zeroes
    for (int i = 0; i < alpha.rows(); i++){
        for (int j = 0; j < alpha.cols(); j++){
            if (alpha(i,j) <= 0){
                stop("Cannot calculate likelihood with values <= 0 in alpha");
            }
        }
    }
    for (int ixi = 0; ixi < ix_len; ++ixi){
        int i = ix(ixi);
        if (sizes(i) == 0){
            stop("Cannot calculate likelihood with values <= 0 in sizes");
        }
    }

#ifdef _OPENMP
    if ( threads > 0 ){
        omp_set_num_threads( threads );
    }
#endif

    int K = alpha.cols();
    Eigen::VectorXd asum(K);
    for (int i = 0; i < K; ++i){
        asum(i) = alpha.col(i).sum();
    }
    Eigen::MatrixXd llk(n_c, K);
    llk.setZero();

    Progress p(n_c * K, display_progress);

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static) collapse(2) 
#endif
    for (int ixi = 0; ixi < ix_len; ++ixi){
        for (int k = 0; k < K; k++){
            int i = ix(ixi);
            if ( ! Progress::check_abort() )
                p.increment();
            llk(i, k) = lgamma(sizes(i) + 1) + 
                lgamma(asum(k)) - 
                lgamma(sizes(i) + asum(k));
            for (cit it(x,i); it; ++it){
                llk(i, k) += lgamma(it.value() + alpha(it.index(), k) ) - 
                    lgamma(it.value() + 1) - 
                    lgamma( alpha(it.index(), k) );
                if (std::isnan(llk(i, k))) {
                    Rcout << "alpha_g = " << alpha(it.index(), k) << "\n";
                    Rcout << "Likelihood value is NA for cell " << i << 
                        " after adding gene " << it.index() << "\n";
                    stop("Exiting");
                }
            }
            if (std::isnan(llk(i, k))) {
                stop("nan value encountered. There are likely non-positive values in x, sizes, or alpha.");
            }
        }
    }

    return(llk);
}


// Estimate Dirichlet multinomial means with method of moments
// [[Rcpp::export]]
Eigen::MatrixXd dirmult_alpha_mom(const Eigen::SparseMatrix<double> &x, 
        const Eigen::MatrixXd &Z, 
        double alpha_prior = 0, 
        double psc = 1e-10, 
        int threads = 1){

#ifdef _OPENMP
    if ( threads > 0 ){
        omp_set_num_threads( threads );
    }
#endif

    int K = Z.cols();
    int N = x.cols();
    int M = x.rows();

    if (Z.rows() != N){
        stop("Number of rows in Z must be number of columns in x");
    }

    Eigen::SparseMatrix<double> p = col_prop(x);
    Eigen::SparseMatrix<double> pt = p.transpose();
    Eigen::MatrixXd pmean(M,K);
    Eigen::MatrixXd pvar(M,K);
    Eigen::MatrixXd alpha(M,K);
    Eigen::VectorXd gmax(K);
    
    // Get pmean
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
    for (int k = 0; k < K; k++){
        pmean.col(k) = fast_wmeanCPP(pt, Z.col(k), 1, false);
        pvar.col(k) = fast_wvarCPP(pt, pmean.col(k), Z.col(k), 1, false);
    }

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
    for (int k = 0; k < K; k++){
        for (int m = 0; m < M; m++){
            double nm = pmean(m,k) * (1 - pmean(m,k));
            alpha(m,k) = (nm / pvar(m,k)) - 1;
            if (std::isnan(alpha(m,k)))
                alpha(m,k) = 0;
            if (alpha(m,k) == std::numeric_limits<double>::infinity())
                alpha(m,k) = 0;
        }
        gmax(k) = alpha(0,k);
        for (int m = 1; m < M; m++){
            if (alpha(m,k) > gmax(k))
                gmax(k) = alpha(m,k);
        }
    }
    for (int k = 0; k < K; k++){
        for (int m = 0; m < M; m++){
            pmean(m,k) = (pmean(m,k) * gmax(k)) + psc + alpha_prior;
        }
    }

    return(pmean);
}

// Estimate Dirichlet multinomial means with leave-one-out
// [[Rcpp::export]]
Eigen::MatrixXd dirmult_alpha_loo(const Eigen::SparseMatrix<double> &x, 
        const Eigen::MatrixXd &alpha0,
        const Eigen::MatrixXd &Z, 
        double eps = 1e-4,
        int max_iter = 1e4,
        double alpha_prior = 0, 
        double psc = 1e-10, 
        int threads = 1){

    int K = Z.cols();
    int N = x.cols();
    int M = x.rows();

    Eigen::VectorXd xsizes = sparse_colsum(x);
    Eigen::SparseMatrix<double> xt = x.transpose();

    Eigen::MatrixXd alpha(alpha0);

    for (int k = 0; k < K; k++){
        alpha.col(k) = max_loo(xt, xsizes, Z.col(k), alpha0.col(k), 
                eps, max_iter, alpha_prior, psc, threads, false);
    }

    return(alpha);
}

// Estimate Dirichlet multinomial means with 2-step
// [[Rcpp::export]]
Eigen::MatrixXd dirmult_alpha_ml(const Eigen::SparseMatrix<double> &x,
        const Eigen::MatrixXd &Z,
        double eps = 1e-4,
        int max_iter_loo = 1e4, 
        double alpha_prior = 0, 
        double psc = 1e-10, 
        int threads = 1){

    if ( Z.rows() != x.cols()){
        stop("The number of rows in Z must be the same as the number of colmns in x");
    }

    Eigen::MatrixXd Alpha0 = dirmult_alpha_mom(x, Z, alpha_prior, psc, threads);
    Eigen::MatrixXd Alpha = dirmult_alpha_loo(x, Alpha0, Z, eps, 
            max_iter_loo, alpha_prior, psc, threads);

    return(Alpha);
}

// MLE of p(z) parameters for Dirichlet-multinomial
// @param X m by n matrix of counts
// @param Alpha m by k matrix of initial probabilities
// @param Pi k vector of probabilities
// @param Z n by k matrix to update in place
// @param fixed 1-based vector of Z row indices to keep fixed
// @param alpha_prior Prior count to add to MLE of Alpha
// @param pi_prior Prior count to add to MLE of Pi
// [[Rcpp::export]]
List dirmult_em(Eigen::SparseMatrix<double> X,
             Eigen::MatrixXd Alpha, 
             Eigen::VectorXd Pi,
             Eigen::MatrixXd Z,
             Eigen::VectorXi fixed,
             double alpha_prior = 0,
             double pi_prior = 0, 
             int threads = 1,
             int max_iter = 100, 
             int max_iter_loo = 1e4, 
             double eps = 1e-4, 
             double eps_loo = 1e-4, 
             double psc = 1e-10, 
             bool display_progress = true){

#ifdef _OPENMP
    if ( threads > 0 ){
        omp_set_num_threads( threads );
    }
    if (display_progress){
        REprintf("Number of threads=%i\n", threads);
    }
#endif

    int N = X.cols();
    int M = X.rows();
    int K = Alpha.cols();

    if (Alpha.rows() != M){
        stop("Number of rows in Alpha must match X");
    }
    if (Pi.size() != K){
        stop("Length of Pi must match number of columns in Alpha");
    }
    if (Z.rows() != N){
        stop("Rows in Z must match columns in X");
    }
    if (Z.cols() != K){
        stop("Columns in Z must match columns in Alpha");
    }

    // Create copies
    Eigen::MatrixXd Ac(Alpha);
    Eigen::MatrixXd Pc(Pi);
    Eigen::MatrixXd Zc(Z);

    // Droplet sizes
    // Eigen::VectorXd sizes(N);
    Eigen::VectorXd sizes = sparse_colsum(X);

    // Get indices of free
    int fixed_n = fixed.size();
    for (int i = 0; i < fixed_n; i++){
        fixed(i) -= 1;
    }

    int free_n = 0;
    Eigen::VectorXi free_o = vec_cmplmnt(fixed, N, free_n);

    if (display_progress){
        Rcout << "\n";
        Rcout << "Number of clusters:\t" << K << "\n";
        Rcout << "Number of features:\t" << M << "\n";
        Rcout << "Total observations:\t" << N << "\n";
        Rcout << "Free observations:\t" << free_n << "\n";
        Rcout << "Fixed observations:\t" << fixed_n << "\n";
        Rcout << "\n";
    }

    // Initial Z update
    Eigen::MatrixXd llk(N,K);
    llk.setZero();

    llk = LlkDirMultSparsePar(X, sizes, Ac, free_o, free_n, 
            threads = threads, false);
    double delta = mult_z_ml(llk, Pc, free_o, free_n, Zc);

    if (display_progress){
        time_t now = std::time(0);
        std::string s(30, '\0');
        std::strftime(&s[0], s.size(), "%Y-%m-%d %H:%M:%S", std::localtime(&now));

        char* dt = ctime(&now);
        Rcout << s << "\tRunning EM\n";
    }

    delta = eps;
    int iter = 1;
    while ((iter <= max_iter) && (delta >= eps)){

        R_CheckUserInterrupt();

        Pc = mult_pi_ml(Zc, pi_prior);
        Ac = dirmult_alpha_ml(X, Zc, eps_loo, max_iter_loo, alpha_prior, psc, threads);
        llk = LlkDirMultSparsePar(X, sizes, Ac, free_o, free_n, 
                threads = threads, false);
        delta = mult_z_ml(llk, Pc, free_o, free_n, Zc);

        if (display_progress){
            time_t now = std::time(0);
            std::string s(30, '\0');
            std::strftime(&s[0], s.size(), "%Y-%m-%d %H:%M:%S", std::localtime(&now));

            char* dt = ctime(&now);
            Rcout << s << "\titeration " << iter << "; delta=" << delta << "\n";
        }

        iter++;
    }

    if (iter > max_iter)
        Rcout << "warning: failed to converge after " << max_iter << " iterations\n";
    else {
        if (display_progress)
            Rcout << "converged after " << iter-1 << " iterations\n";
    }

    Rcout << "Done\n";
    List l = List::create(Ac, Pc, Zc, llk);

    return l;
}

