
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
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppEigen)]]

typedef Eigen::SparseMatrix<double>::InnerIterator cit;

// Compute log density of Dirichlet-Multinomial PMF with sparse matrix
// x is a gene by droplet sparse matrix
// [[Rcpp::export]]
Eigen::MatrixXd LlkMultSparsePar(const Eigen::SparseMatrix<double> &x, 
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

#ifdef _OPENMP
    if ( threads > 0 ){
        omp_set_num_threads( threads );
    }
#endif

    // Check for zeroes
    for (int i = 0; i < alpha.rows(); i++){
        for (int j = 0; j < alpha.cols(); j++){
            if (alpha(i,j) <= 0){
                stop("Cannot calculate likelihood with values <= 0 in alpha");
            }
        }
    }
    // testing
    if (ix.size() != ix_len){
        stop("ix.size != ix_len");
    }
    for (int ixi = 0; ixi < ix_len; ++ixi){
        int i = ix(ixi);
        if (sizes(i) == 0){
            stop("Cannot calculate likelihood with values <= 0 in sizes");
        }
    }

    int K = alpha.cols();
    Eigen::MatrixXd llk(n_c, K);
    llk.setZero();

    Progress p(n_c * K, display_progress);

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 100) collapse(2) 
#endif
    for (int ixi = 0; ixi < ix_len; ++ixi){
        for (int k = 0; k < K; k++){
            int i = ix(ixi);
            if (i >= n_c){
                stop("Indices must be between 0 and number of columns of x - 1");
            }
            if ( ! Progress::check_abort() )
                p.increment();
            llk(i,k) = std::lgamma(sizes(i) + 1);
            for (Eigen::SparseMatrix<double>::InnerIterator it(x,i); it; ++it){
                llk(i,k) -= std::lgamma(it.value() + 1);
                llk(i,k) += it.value() * log(alpha(it.index(), k));
            }
            if (std::isnan(llk(i,k))) {
                stop("nan value encountered. There are likely non-positive values in x, sizes, or alpha.");
            }
        }
    }

    return(llk);
}

// MLE of alpha parameters for multinomial
// [[Rcpp::export]]
Eigen::MatrixXd mult_alpha_ml(const Eigen::SparseMatrix<double> &x, 
        const Eigen::MatrixXd &Z, 
        double prior = 0, 
        int threads = 1, 
        double psc = 1e-10){ 
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

    Eigen::MatrixXd A(M, K);
    A.setZero();

    for (int k = 0; k < K; k++){
        for (int m = 0; m < M; m++){
            A(m,k) = prior + psc;
        }
    }

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 100)
#endif
    for (int i = 0; i < N; i++){
        for (cit it(x, i); it; ++it){
            int m = it.index();
            for (int k = 0; k < K; k++){
                A(m,k) += (Z(i,k) * it.value());
            }
        }
    }

    Eigen::VectorXd Asum(K);
    Asum.setZero();
    for (int k = 0; k < K; k++){
        for (int m = 0; m < M; m++){
            Asum(k) += A(m,k);
        }
    }

    for (int k = 0; k < K; k++){
        for (int m = 0; m < M; m++){
            A(m,k) /= Asum(k);
        }
    }

    return A;
}

// MLE of p(z) parameters for multinomial
// @param llk n by k matrix of log likelihoods
// @param pi k vector of pi probabilities
// @param ix vector of row indices to update for llk
// @param ix_len length of ix
// @param Z n by k matrix to update in place
// [[Rcpp::export]]
double mult_z_ml(const Eigen::MatrixXd &llk, 
        const Eigen::VectorXd &pi, 
        const Eigen::VectorXi &ix, 
        const int &ix_len, 
        Eigen::MatrixXd &Z, 
        int threads = 1){ 
#ifdef _OPENMP
    if ( threads > 0 ){
        omp_set_num_threads( threads );
    }
#endif

    int K = Z.cols();
    double delta = 0;

    if (pi.size() != K){
        stop("Length of pi must be number of columns in llk");
    }

    if (llk.rows() != Z.rows()){
        stop("Dimensions of llk and Z must match");
    }

    if (llk.cols() != Z.cols()){
        stop("Dimensions of llk and Z must match");
    }

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 100) reduction(+:delta)
#endif
    for (int ixi = 0; ixi < ix_len; ixi++){
        Eigen::VectorXd Zk(K);
        int i = ix(ixi);
        for (int k = 0; k < K; k++){
            Zk(k) = llk(i,k) + std::log(pi(k));
        }
        Zk = fraction_log(Zk);
        double d1 = 0;
        for (int k = 0; k < K; k++){
            d1 += std::fabs(Zk(k) - Z(i,k));
            Z(i,k) = Zk(k);
        }
        delta += d1;
    }

    delta /= ix_len;

    return delta;
}

// MLE of p(z) parameters for multinomial
// @param llk n by k matrix of log likelihoods
// @param pi k vector of pi probabilities
// @param ix vector of row indices to update for llk
// @param ix_len length of ix
// @param Z n by k matrix to update in place
// [[Rcpp::export]]
Eigen::VectorXd mult_pi_ml(const Eigen::MatrixXd &Z, 
        double prior = 0,
        int threads = 1){
#ifdef _OPENMP
    if ( threads > 0 ){
        omp_set_num_threads( threads );
    }
#endif

    int N = Z.rows();
    int K = Z.cols();
    Eigen::VectorXd Pi(K);

    for (int k = 0; k < K; k++){
        Pi(k) = prior;
    }

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static) 
#endif
    for (int k = 0; k < K; k++){
        for (int i = 0; i < N; i++){
            Pi(k) += Z(i,k);
        }
    }

    Pi = fraction_vec(Pi);

    return Pi;
}

// MLE of p(z) parameters for multinomial
// @param X m by n matrix of counts
// @param Alpha m by k matrix of initial probabilities
// @param Pi k vector of probabilities
// @param Z n by k matrix to update in place
// @param fixed 1-based vector of Z row indices to keep fixed
// @param alpha_prior Prior count to add to MLE of Alpha
// @param pi_prior Prior count to add to MLE of Pi
// [[Rcpp::export]]
List mult_em(Eigen::SparseMatrix<double> X,
             Eigen::MatrixXd Alpha, 
             Eigen::VectorXd Pi,
             Eigen::MatrixXd Z,
             Eigen::VectorXi fixed,
             double alpha_prior = 0,
             double pi_prior = 0, 
             int threads = 1,
             int max_iter = 100, 
             double eps = 1e-4, 
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

    int N = X.cols(); // M droplets
    int M = X.rows(); // M genes
    int K = Alpha.cols(); // K groups

    if (Alpha.rows() != M || Alpha.cols() != K){
        stop("Dimensions of Alpha must be (M,K)");
    }
    if (Pi.size() != K){
        stop("Length of Pi must match number of columns in Alpha");
    }
    if (Z.rows() != N || Z.cols() != K){
        stop("Dimensions of Z must match (N,K)");
    }

    // Create copies
    Eigen::MatrixXd Ac(Alpha);
    Eigen::MatrixXd Pc(Pi);
    Eigen::MatrixXd Zc(Z);

    // Droplet sizes
    // Eigen::VectorXd sizes(N);
    Eigen::VectorXd sizes = sparse_colsum(X);

    // Check fixed indices
    int fixed_n = fixed.size();
    for (int i = 0; i < fixed_n; i++){
        fixed(i) -= 1;
        if (fixed(i) < 0 || fixed(i) >= N)
            stop("fixed indices must be >= 0 and < N");
    }

    // get free indices
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

    llk = LlkMultSparsePar(X, sizes, Ac, free_o, free_n, 
            threads = threads, false);
    double delta = mult_z_ml(llk, Pc, free_o, free_n, Zc, threads);

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

        Pc = mult_pi_ml(Zc, pi_prior, threads);
        Ac = mult_alpha_ml(X, Zc, alpha_prior, threads, psc);
        llk = LlkMultSparsePar(X, sizes, Ac, free_o, free_n, threads, false);
        delta = mult_z_ml(llk, Pc, free_o, free_n, Zc, threads);

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

    if (display_progress)
        Rcout << "Done\n";

    List l = List::create(Ac, Pc, Zc, llk, iter);

    return l;
}

